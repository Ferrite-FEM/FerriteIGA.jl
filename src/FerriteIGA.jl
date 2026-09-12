module FerriteIGA

using Reexport

@reexport using Tensors
@reexport using Ferrite
using Ferrite: 
    AbstractRefShape, RefHypercube, RefLine, RefQuadrilateral, RefHexahedron, getnbasefunctions,
    VectorInterpolation, VectorizedInterpolation,
    FunctionValues, GeometryMapping, MappingValues, ValuesUpdateFlags

using OrderedCollections:
    OrderedSet

using WriteVTK
using LinearAlgebra
using StaticArrays
import SparseArrays

export IGAInterpolation
export BezierExtractionOperator
export BezierCell
export BezierCoords
export VTKIGAFile
export IGACellCache, IGAFaceCache

const Optional{T} = Union{T, Nothing}
const BezierExtractionOperator{T} = Vector{SparseArrays.SparseVector{T,Int}}
const CoordsAndWeight{sdim,T} = Tuple{ <: AbstractVector{Vec{sdim,T}}, <: AbstractVector{T}}

"""
    BezierCoords

Data for one cell, from `getcoordinates`.
`x`, `w` are the NURBS points and weights.
`xb`, `wb` are the Bézier points and weights.
`beo` is the extraction operator.
"""
struct BezierCoords{dim_s,T} 
    xb   ::Vector{Vec{dim_s,T}}
    wb   ::Vector{T}
    x    ::Vector{Vec{dim_s,T}}
    w    ::Vector{T}
    beo  ::Base.RefValue{ BezierExtractionOperator{T} }
end

function resize_bezier_coord!(X::BezierCoords, N::Int)
    (; xb, wb, x, w) = X
    resize!(xb, N)
    resize!(wb, N)
    resize!(x,  N)
    resize!(w,  N)
end

zero_bezier_coord(dim, T, nnodes) = BezierCoords{dim,T}(zeros(Vec{dim,T}, nnodes), zeros(T, nnodes), zeros(Vec{dim,T}, nnodes), zeros(T, nnodes), Base.RefValue(diagonal_beo(1)))

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

include("nurbsmesh.jl")
include("mesh_generation.jl")
include("bezier_grid.jl")
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