module FerriteIGA

using Reexport: @reexport
@reexport using Tensors
@reexport using Ferrite

using Ferrite: 
    AbstractRefShape, RefHypercube, RefLine, RefQuadrilateral, RefHexahedron, getnbasefunctions,
    VectorInterpolation, VectorizedInterpolation,
    FunctionValues, GeometryMapping, MappingValues, ValuesUpdateFlags
using OrderedCollections:
    OrderedSet
using WriteVTK:
    WriteVTK, VTKCellTypes
using LinearAlgebra: 
    LinearAlgebra, Diagonal, norm, qr, ⋅
using SparseArrays: 
    SparseArrays, SparseMatrixCSC, SparseVector, sparsevec, sparse

const Optional{T} = Union{T, Nothing}
const CoordsAndWeight{sdim,T} = Tuple{ <: AbstractVector{Vec{sdim,T}}, <: AbstractVector{T}}

"""
    BezierExtractionOperator{T}

Represent the Bézier extraction operator for an element, used to map between the local Bernstein polynomial basis and the global B-spline or NURBS basis.
"""
struct BezierExtractionOperator{T}
    C::Vector{SparseVector{T,Int}} #TODO: Consider changing this to SparseMatrixCSR instead
end

include("exports.jl")

"""
    BezierCoords

Data for one cell, from `getcoordinates`.
`x`, `w` are the NURBS points and weights.
`xb`, `wb` are the Bézier points and weights.
`beo` is the extraction operator.
"""
mutable struct BezierCoords{dim_s,T} 
    const xb::Vector{Vec{dim_s,T}}
    const wb::Vector{T}
    const x::Vector{Vec{dim_s,T}}
    const w::Vector{T}
    beo::Optional{BezierExtractionOperator{T}}
end

#`nb` is the number of Bézier points (Bernstein functions), and `n` is the number of 
#control points of the cell. These are equal for NURBS cells, but not for T-spline cells.
function resize_bezier_coord!(X::BezierCoords, nb::Int, n::Int = nb)
    (; xb, wb, x, w) = X
    resize!(xb, nb)
    resize!(wb, nb)
    resize!(x,  n)
    resize!(w,  n)
end

zero_bezier_coord(dim, T, nb, n = nb) = BezierCoords{dim,T}(zeros(Vec{dim,T}, nb), zeros(T, nb), zeros(Vec{dim,T}, n), zeros(T, n), nothing)

#Base.zero(Type{BezierCoords{dim,T}}) where {dim,T} = BezierCoords

"""
    IGAInterpolation{shape, order}()

Bernstein basis of degree `order` on `shape` (`RefLine`, `RefQuadrilateral`, or `RefHexahedron`).
Use with `BezierCellValues`. For a 2D vector unknown write `ip^2`.
"""
struct IGAInterpolation{shape, order} <: Ferrite.ScalarInterpolation{shape, order}
    function IGAInterpolation{shape,order}() where {rdim, shape<:RefHypercube{rdim}, order} 
        #Check if you can construct a Bernstein basis
        Bernstein{shape,order}()
        return new{shape,order}()
    end
end

Ferrite.adjust_dofs_during_distribution(::IGAInterpolation) = false
Ferrite.reference_shape_value(::IGAInterpolation{shape, order}, ξ::Vec, i::Int) where {shape, order} = Ferrite.reference_shape_value(Bernstein{shape, order}(), ξ, i)
Ferrite.reference_coordinates(::IGAInterpolation{shape, order}) where {shape, order} = Ferrite.reference_coordinates(Bernstein{shape, order}())
Ferrite.getnbasefunctions(::IGAInterpolation{shape, order}) where {shape, order} = getnbasefunctions(Bernstein{shape, order}())

Ferrite.nvertices(ip::IGAInterpolation) = getnbasefunctions(ip)
Ferrite.vertexdof_indices(ip::IGAInterpolation) = ntuple(i->i, getnbasefunctions(ip))

#Remove dofs on edges and faces such that dofhandler can distribute dofs correctly
Ferrite.edgedof_indices(ip::IGAInterpolation{RefHexahedron}) = ntuple(_ -> (), Ferrite.nedges(ip))
Ferrite.edgedof_interior_indices(ip::IGAInterpolation) = ntuple(_ -> (), Ferrite.nedges(ip))
Ferrite.facedof_indices(ip::IGAInterpolation{RefHexahedron}) =  ntuple(_ -> (), Ferrite.nfaces(ip))
Ferrite.facedof_interior_indices(ip::IGAInterpolation{RefHexahedron}) =  ntuple(_ -> (), Ferrite.nfaces(ip))
Ferrite.facedof_indices(ip::IGAInterpolation{RefQuadrilateral}) =  ntuple(_ -> (), Ferrite.nfaces(ip))
Ferrite.facedof_interior_indices(ip::IGAInterpolation{RefQuadrilateral}) =  ntuple(_ -> (), Ferrite.nfaces(ip))

Ferrite.dirichlet_facedof_indices(::IGAInterpolation{shape, order}) where {shape, order} = Ferrite.dirichlet_facedof_indices(Bernstein{shape, order}())
Ferrite.dirichlet_edgedof_indices(::IGAInterpolation{shape, order}) where {shape, order} = Ferrite.dirichlet_edgedof_indices(Bernstein{shape, order}())
#Ferrite.dirichlet_vertexdof_indices(::IGAInterpolation{shape, order}) where {shape, order} = Ferrite.dirichlet_vertexdof_indices(Bernstein{shape, order}())



"""
    BezierCell{shape, order}(nodes)

One cell of a `BezierGrid`. `nodes` are control-point numbers. `order` is the degree.
"""
struct BezierCell{refshape,order,N} <: Ferrite.AbstractCell{refshape}
    nodes::NTuple{N,Int}
    function BezierCell{refshape,order}(nodes::NTuple{N,Int}) where {rdim, order, refshape<:RefHypercube{rdim}, N} 
        @assert(order isa Integer)
        @assert (order+1)^rdim == N
		return new{refshape,order,N}(nodes)
    end
end

function BezierCell{refshape,order,N}(nodes::NTuple{N,Int}) where {rdim, order, refshape<:RefHypercube{rdim}, N} 
    return BezierCell{refshape,order}(nodes)
end

Ferrite.geometric_interpolation(::Type{<:BezierCell{shape,order}}) where {order, shape} = IGAInterpolation{shape, order}()
Ferrite.nnodes(::BezierCell{shape, order, N}) where {order, shape, N} = N
getorders(::BezierCell{refshape, orders, N}) where {orders,refshape,N} = orders

"""
    TSplineCell{shape, order}(nodes)

One Bézier element of a T-spline `BezierGrid`. `nodes` are the global control points (T-spline
basis functions) that are non-zero on the element. Unlike `BezierCell`, the number of nodes `N` 
can be larger than the number of Bernstein functions, `(order+1)^dim`, and the nodes
have no particular ordering.
"""
struct TSplineCell{refshape,order,N} <: Ferrite.AbstractCell{refshape}
    nodes::NTuple{N,Int}
    function TSplineCell{refshape,order}(nodes::NTuple{N,Int}) where {rdim, order, refshape<:RefHypercube{rdim}, N}
        @assert(order isa Integer)
        @assert N >= (order+1)^rdim
        return new{refshape,order,N}(nodes)
    end
end
TSplineCell{refshape,order,N}(nodes::NTuple{N,Int}) where {refshape,order,N} = TSplineCell{refshape,order}(nodes)

"""
    TSplineInterpolation{shape, order, N}()

The `N` spline basis functions that are non-zero on a Bézier element of a T-spline grid.
This interpolation is only used for distributing dofs and applying boundary conditions;
the basis functions are always evaluated as the Bernstein basis (`IGAInterpolation{shape,order}`)
multiplied with the Bézier extraction operator of the cell.

Normally, this does not need to be created manually, since `add!(dh, :u, IGAInterpolation{shape,order}())`
selects the correct `TSplineInterpolation` for each cell type of a T-spline `BezierGrid`.
"""
struct TSplineInterpolation{shape, order, N} <: Ferrite.ScalarInterpolation{shape, order}
    function TSplineInterpolation{shape,order,N}() where {rdim, shape<:RefHypercube{rdim}, order, N}
        @assert N >= (order+1)^rdim
        return new{shape,order,N}()
    end
end

Ferrite.adjust_dofs_during_distribution(::TSplineInterpolation) = false
Ferrite.getnbasefunctions(::TSplineInterpolation{shape,order,N}) where {shape,order,N} = N
Ferrite.nvertices(::TSplineInterpolation{shape,order,N}) where {shape,order,N} = N
Ferrite.vertexdof_indices(::TSplineInterpolation{shape,order,N}) where {shape,order,N} = ntuple(i->(i,), N)
Ferrite.edgedof_indices(ip::TSplineInterpolation) = ntuple(_ -> (), Ferrite.nedges(ip))
Ferrite.edgedof_interior_indices(ip::TSplineInterpolation) = ntuple(_ -> (), Ferrite.nedges(ip))
Ferrite.facedof_indices(ip::TSplineInterpolation) = ntuple(_ -> (), Ferrite.nfaces(ip))
Ferrite.facedof_interior_indices(ip::TSplineInterpolation) = ntuple(_ -> (), Ferrite.nfaces(ip))
Ferrite.volumedof_interior_indices(::TSplineInterpolation) = ()
#Which basis functions are non-zero on a facet differs between cells, so facet based dirichlet conditions 
#are converted to node based conditions instead (see `add!(::ConstraintHandler, ::Dirichlet)`).
Ferrite.dirichlet_facedof_indices(ip::TSplineInterpolation) = ntuple(_ -> (), Ferrite.nfacets(ip))
Ferrite.dirichlet_edgedof_indices(ip::TSplineInterpolation) = ntuple(_ -> (), Ferrite.nedges(ip))
#The T-spline basis functions are not nodal. Placeholder coordinates are returned since Ferrite 
#requires them when setting up (unused) boundary values for node based dirichlet conditions.
Ferrite.reference_coordinates(::TSplineInterpolation{shape,order,N}) where {dim,shape<:AbstractRefShape{dim},order,N} = zeros(Vec{dim,Float64}, N)
function Ferrite.reference_shape_value(ip::TSplineInterpolation, ::Vec, ::Int)
    throw(ArgumentError("$ip can not be evaluated directly. Use BezierCellValues together with the Bézier extraction operator."))
end

Ferrite.geometric_interpolation(::Type{<:TSplineCell{shape,order,N}}) where {shape,order,N} = TSplineInterpolation{shape,order,N}()
Ferrite.nnodes(::TSplineCell{shape, order, N}) where {order, shape, N} = N
getorders(::TSplineCell{refshape, order, N}) where {order,refshape,N} = order
Ferrite.vertices(c::TSplineCell) = c.nodes

function _tspline_topology_error()
    throw(ArgumentError("The facets/edges of a TSplineCell can not be computed from its nodes. " * 
        "Use the facetsets and nodesets created by `BezierGrid(::TSplineMesh)` instead."))
end
Ferrite.facets(::TSplineCell) = _tspline_topology_error()
#Ferrite loops over the edges and faces when distributing dofs. A TSplineInterpolation has no 
#dofs on edges and faces, so empty placeholders are returned.
Ferrite.edges(c::TSplineCell) = ntuple(_ -> (), Ferrite.nedges(Ferrite.getrefshape(c)))
Ferrite.faces(c::TSplineCell) = ntuple(_ -> (), Ferrite.nfaces(Ferrite.getrefshape(c)))

"""
    _bernstein_ip(ip)

Return the Bernstein (`IGAInterpolation`) basis used on the Bézier element. Identity for other interpolations.
"""
_bernstein_ip(ip::Interpolation) = ip
_bernstein_ip(::TSplineInterpolation{shape,order}) where {shape,order} = IGAInterpolation{shape,order}()
_bernstein_ip(ip::VectorizedInterpolation{vdim}) where {vdim} = VectorizedInterpolation{vdim}(_bernstein_ip(ip.ip))

include("nurbsmesh.jl")
include("tsplines/tsplinemesh.jl")
include("mesh_generation.jl")
include("bezier_grid.jl")
include("dofhandler.jl")
include("constraints.jl")
include("splines/bezier.jl")
include("bezier_extraction.jl")
include("splines/bezier_values.jl")
include("splines/bsplines.jl")
include("VTK.jl")
include("iterators.jl")
include("apply_analytical_iga.jl")
include("L2_projection.jl")

Ferrite._mass_qr(::IGAInterpolation{shape,order}) where {shape,order}= Ferrite._mass_qr(Bernstein{shape, order}())
Ferrite._mass_qr(::Bernstein{RefQuadrilateral, 2}) = QuadratureRule{RefQuadrilateral}(2+1)
Ferrite._mass_qr(::Bernstein{RefHexahedron, 2}) = QuadratureRule{RefCube}(2+1)

#Normaly the verices function should only return the 8 corner nodes of the hexa (or 4 in 2d),
#but since the cell connectivity in IGA is different compared to normal FE elements,
#we can only distribute cells on the nodes/controlpoints
Ferrite.vertices(c::BezierCell) = c.nodes

_bernstein_ordering(::Type{<:BezierCell{shape,order}}) where {order,shape} = _bernstein_ordering(Bernstein{shape,order}())                                        

function Ferrite.faces(c::BezierCell{shape,order}) where {shape,order}
    return getindex.(Ref(c.nodes), collect.(Ferrite.dirichlet_facedof_indices(IGAInterpolation{shape,order}() )))
end

function Ferrite.edges(c::BezierCell{shape, order}) where {shape,order}
    return getindex.(Ref(c.nodes), collect.(Ferrite.dirichlet_edgedof_indices(IGAInterpolation{shape,order}() )))
end

end #end module