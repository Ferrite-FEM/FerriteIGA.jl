"""
    IGACellCache(dh::DofHandler, flags::UpdateFlags = UpdateFlags())

The isogeometric counterpart of Ferrite's `CellCache`, for a `DofHandler` defined on a
[`BezierGrid`](@ref). It holds pre-allocated storage for the data of one cell, and is
updated for cell `cellid` with `reinit!(cache, cellid)`. After that, the following
functions can be used:

 - `getcoordinates(cache)`: the [`BezierCoords`](@ref) of the cell,
 - `celldofs(cache)`: the global degrees of freedom of the cell,
 - `getnodes(cache)`: the node (control point) indices of the cell,
 - `cellid(cache)`: the index of the cell.

`flags` determines which of the data is updated, see Ferrite's `UpdateFlags`.

# Examples
```julia
cc = IGACellCache(dh)
for cellid in 1:getncells(grid)
    reinit!(cc, cellid)
    reinit!(cellvalues, getcoordinates(cc))
    assemble_element!(ke, cellvalues)
    assemble!(assembler, celldofs(cc), ke)
end
```
"""
struct IGACellCache{X,G<:Ferrite.AbstractGrid,DH<:Union{Ferrite.AbstractDofHandler,Nothing}}
    flags::UpdateFlags
    grid::G
    # Pretty useless to store this since you have it already for the reinit! call, but
    # needed for the CellIterator(...) workflow since the user doesn't necessarily control
    # the loop order in the cell subset.
    cellid::Base.RefValue{Int}
    nodes::Vector{Int}
    bezier_cell_data::X
    dh::DH
    dofs::Vector{Int}
end

function IGACellCache(dh::DofHandler{dim,G}, flags::UpdateFlags=UpdateFlags()) where {dim,G}
    grid = Ferrite.get_grid(dh)
    T = Ferrite.get_coordinate_eltype(grid)
    N = Ferrite.nnodes_per_cell(grid, 1)
    nodes    = zeros(Int, N)
    coords   = zero_bezier_coord(dim,T,N)
    n        = ndofs_per_cell(dh, 1)
    celldofs = zeros(Int, n)

    return IGACellCache(flags, grid, Ref(-1), nodes, coords, dh, celldofs)
end

function Ferrite.reinit!(cc::IGACellCache, i::Int)
    cc.cellid[] = i
    if cc.flags.nodes
        resize!(cc.nodes, Ferrite.nnodes_per_cell(cc.grid, i))
        Ferrite.cellnodes!(cc.nodes, cc.grid, i)
    end
    if cc.flags.coords
        resize_bezier_coord!(cc.bezier_cell_data, Ferrite.nnodes_per_cell(cc.grid, i))
        getcoordinates!(cc.bezier_cell_data, cc.grid, i)
    end
    if cc.dh !== nothing && cc.flags.dofs
        resize!(cc.dofs, ndofs_per_cell(cc.dh, i))
        celldofs!(cc.dofs, cc.dh, i)
    end
    return cc
end

# Accessor functions (TODO: Deprecate? We are so inconsistent with `getxx` vs `xx`...)
Ferrite.getnodes(cc::IGACellCache) = cc.nodes
Ferrite.getcoordinates(cc::IGACellCache) = cc.bezier_cell_data
Ferrite.celldofs(cc::IGACellCache) = cc.dofs
Ferrite.cellid(cc::IGACellCache) = cc.cellid[]


"""
    IGAFaceCache(dh::DofHandler, flags::UpdateFlags = UpdateFlags())

The isogeometric counterpart of Ferrite's `FacetCache`, for a `DofHandler` defined on a
[`BezierGrid`](@ref). It is updated with `reinit!(cache, FacetIndex(cellid, facetnr))`,
after which `getcoordinates`, `celldofs`, `getnodes` and `cellid` return the data of the
cell that the facet belongs to, see [`IGACellCache`](@ref).
"""
struct IGAFaceCache{CC}
    cc::CC  # const for julia > 1.8
    dofs::Vector{Int} # aliasing cc.dofs
    current_faceid::Base.RefValue{Int}
end

function IGAFaceCache(args...)
    cc = IGACellCache(args...)
    IGAFaceCache(cc, cc.dofs, Ref(-1))
end

function Ferrite.reinit!(fc::IGAFaceCache, face::Ferrite.BoundaryIndex)
    cellid, faceid = face
    reinit!(fc.cc, cellid)
    fc.current_faceid[] = faceid
    return nothing
end

# Delegate methods to the cell cache
for op = (:getnodes, :getcoordinates, :cellid, :celldofs)
    @eval begin
        function Ferrite.$op(fc::IGAFaceCache, args...)
            return Ferrite.$op(fc.cc, args...)
        end
    end
end
# @inline faceid(fc::FaceCache) = fc.current_faceid[]
@inline Ferrite.celldofs!(v::Vector, fc::IGAFaceCache) = celldofs!(v, fc.cc)