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

The Bézier extraction operator ``\\mathbf{C}^e`` of one cell. It expresses the B-spline
basis functions that are non-zero on the cell as linear combinations of the Bernstein
polynomials of the same degree,

```math
\\mathbf{N}^e(\\boldsymbol{\\xi}) = \\mathbf{C}^e \\mathbf{B}(\\boldsymbol{\\xi}).
```

The operator is stored row-wise (one sparse vector per B-spline basis function on the
cell) in the field `C`. It is computed from the knot vectors when a [`BezierGrid`](@ref) is
created from a [`NURBSMesh`](@ref), and is obtained for a cell with
[`get_extraction_operator`](@ref) or as the `beo` field of a [`BezierCoords`](@ref).

See the [BezierGrid topic guide](@ref beziergrid_topicguide) for more details.
"""
struct BezierExtractionOperator{T}
    C::Vector{SparseVector{T,Int}} #TODO: Consider changing this to SparseMatrixCSR instead
end

include("exports.jl")

"""
    BezierCoords{sdim, T}

The geometric data of one cell of a [`BezierGrid`](@ref), returned by
`getcoordinates(grid::BezierGrid, cellid)`. In a standard Ferrite program
`getcoordinates` returns a vector of node coordinates; for an isogeometric cell more
information is needed to [`reinit!`](@ref reinit_bezier) a [`BezierCellValues`](@ref), and
`BezierCoords` bundles all of it:

 - `x::Vector{Vec{sdim,T}}`: the NURBS control points of the cell,
 - `w::Vector{T}`: the NURBS weights of the cell,
 - `xb::Vector{Vec{sdim,T}}`: the Bézier control points of the cell,
 - `wb::Vector{T}`: the Bézier weights of the cell,
 - `beo`: the [`BezierExtractionOperator`](@ref) of the cell (or `nothing` for cells
   without one).

A `BezierCoords` can be passed directly to `reinit!`:

```julia
coords = getcoordinates(grid, cellid)
reinit!(cellvalues, coords)
```

A pre-allocated `BezierCoords` can be updated in-place with
`getcoordinates!(coords, grid, cellid)`.
"""
mutable struct BezierCoords{dim_s,T} 
    const xb::Vector{Vec{dim_s,T}}
    const wb::Vector{T}
    const x::Vector{Vec{dim_s,T}}
    const w::Vector{T}
    beo::Optional{BezierExtractionOperator{T}}
end

function resize_bezier_coord!(X::BezierCoords, N::Int)
    (; xb, wb, x, w) = X
    resize!(xb, N)
    resize!(wb, N)
    resize!(x,  N)
    resize!(w,  N)
end

zero_bezier_coord(dim, T, nnodes) = BezierCoords{dim,T}(zeros(Vec{dim,T}, nnodes), zeros(T, nnodes), zeros(Vec{dim,T}, nnodes), zeros(T, nnodes), nothing)

#Base.zero(Type{BezierCoords{dim,T}}) where {dim,T} = BezierCoords

"""
    IGAInterpolation{refshape, order}()

The interpolation used for isogeometric analysis with FerriteIGA. On the reference cell it
is the tensor-product Bernstein basis of degree `order` on `refshape` (`RefLine`,
`RefQuadrilateral` or `RefHexahedron`). Combined with the Bézier extraction operators
stored in a [`BezierGrid`](@ref), it represents the B-spline or NURBS basis of the patch.

`IGAInterpolation` is used wherever a Ferrite interpolation (e.g. `Lagrange`) is used in a
standard Ferrite program: when adding fields to a `DofHandler`, and when constructing
[`BezierCellValues`](@ref) and [`BezierFacetValues`](@ref). Use `ip^dim` to obtain a
vector-valued interpolation.

The order must match the order of the cells of the grid, and is currently the same in all
parametric directions.

# Examples
```julia
ip = IGAInterpolation{RefQuadrilateral, 2}()  # scalar, quadratic
ip_u = ip^2                                   # vector-valued, quadratic
dh = DofHandler(grid)
add!(dh, :u, ip_u)
close!(dh)
```

See the [BezierCellValues topic guide](@ref beziercellvalues_topicguide) for more details.
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
    BezierCell{refshape, order}(nodes::NTuple{N, Int})

A cell of a [`BezierGrid`](@ref), i.e. one non-empty knot span of a spline patch, with
reference shape `refshape` (`RefLine`, `RefQuadrilateral` or `RefHexahedron`) and polynomial
degree `order`. The `nodes` are the indices of the `N = (order + 1)^dim` control points
whose basis functions are non-zero on the cell, ordered according to the Bernstein basis of
[`IGAInterpolation`](@ref).

Unlike a Lagrange cell, neighboring `BezierCell`s share several layers of nodes, reflecting
the higher continuity of the spline basis.

`BezierCell` can be passed to Ferrite's `generate_grid` to create a uniform B-spline grid of
a line, rectangle or box directly:

```julia
grid = generate_grid(BezierCell{RefQuadrilateral, 2}, (10, 10), Vec((0.0, 0.0)), Vec((1.0, 1.0)))
```
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