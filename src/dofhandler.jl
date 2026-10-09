#
# Support for T-spline grids in Ferrite's DofHandler and ConstraintHandler.
#
# The cells of a T-spline grid have different numbers of basis functions, so they have different types
# (TSplineCell{shape,order,N}). Ferrite requires a single cell type and number of dofs per SubDofHandler,
# so the cells are grouped by type, with one SubDofHandler per group.
#

_has_tspline_cells(grid::BezierGrid) = any(c -> c isa TSplineCell, grid.cells)

"""
    _localize(ip, celltype)

Return the interpolation which is distributed on cells of type `celltype`. For a `TSplineCell` with
`N` nodes, `IGAInterpolation{shape,order}` is replaced with `TSplineInterpolation{shape,order,N}`.
"""
_localize(ip::Interpolation, ::Type{<:Ferrite.AbstractCell}) = ip
_localize(ip::VectorizedInterpolation{vdim}, CT::Type{<:TSplineCell}) where {vdim} = VectorizedInterpolation{vdim}(_localize(ip.ip, CT))
_localize(::Union{IGAInterpolation{shape,order}, TSplineInterpolation{shape,order}}, ::Type{<:TSplineCell{shape,order,N}}) where {shape,order,N} = TSplineInterpolation{shape,order,N}()
function _localize(ip::Interpolation, CT::Type{<:TSplineCell})
    throw(ArgumentError("The interpolation $ip can not be used on cells of type $CT. Use an IGAInterpolation with the same shape and order as the T-spline."))
end

"""
    add!(dh::DofHandler{dim, <:BezierGrid}, name::Symbol, ip::Interpolation)

Add a field to all cells of a `BezierGrid`. For T-spline grids, which have cells with different numbers of basis functions,
one `SubDofHandler` is created per cell type, and `ip` (an `IGAInterpolation` with the same order as the T-spline) is replaced
with the corresponding `TSplineInterpolation`. Note that the number of dofs then varies between cells, so
`ndofs_per_cell(dh, cellid)` must be used instead of `ndofs_per_cell(dh)`.
"""
function Ferrite.add!(dh::DofHandler{<:Any,<:BezierGrid}, name::Symbol, ip::Interpolation)
    grid = Ferrite.get_grid(dh)
    if !_has_tspline_cells(grid)
        return invoke(Ferrite.add!, Tuple{DofHandler, Symbol, Interpolation}, dh, name, ip)
    end
    @assert !Ferrite.isclosed(dh)

    #Group cells by type
    celltypes = DataType[]
    cellsets  = OrderedSet{Int}[]
    for (cellid, cell) in enumerate(getcells(grid))
        k = findfirst(==(typeof(cell)), celltypes)
        if k === nothing
            push!(celltypes, typeof(cell))
            push!(cellsets, OrderedSet{Int}())
            k = length(celltypes)
        end
        push!(cellsets[k], cellid)
    end

    if isempty(dh.subdofhandlers)
        for cellset in cellsets
            SubDofHandler(dh, cellset)
        end
    elseif sum(sdh -> length(sdh.cellset), dh.subdofhandlers) != getncells(grid)
        error("can not add field to DofHandler with SubDofHandlers which do not cover the whole grid")
    end

    for sdh in dh.subdofhandlers
        add!(sdh, name, _localize(ip, Ferrite.getcelltype(sdh)))
    end
    return dh
end

"""
    add!(ch::ConstraintHandler{<:DofHandler{dim, <:BezierGrid}}, dbc::Dirichlet)

Add a Dirichlet boundary condition. For T-spline grids, a condition on a facetset is converted to a
condition on the control points whose basis functions are non-zero on the facets. The prescribed value
of each control point is `f` evaluated at the control point.
"""
function Ferrite.add!(ch::ConstraintHandler{<:DofHandler{<:Any,<:BezierGrid}}, dbc::Dirichlet)
    grid = Ferrite.get_grid(ch.dh)
    if _has_tspline_cells(grid) && eltype(dbc.facets) <: FacetIndex
        nodes = _facets_to_nodes(grid, dbc.facets)
        components = isempty(dbc.components) ? nothing : dbc.components
        dbc = Dirichlet(dbc.field_name, nodes, dbc.f, components)
    end
    return invoke(Ferrite.add!, Tuple{ConstraintHandler, Dirichlet}, ch, dbc)
end

"""
    _facets_to_nodes(grid::BezierGrid, facets)

Return the nodes (control points) whose basis functions are non-zero on any of the `facets`.
"""
function _facets_to_nodes(grid::BezierGrid, facets)
    nodes = OrderedSet{Int}()
    for (cellid, facetid) in facets
        cell = getcells(grid, cellid)
        if cell isa TSplineCell
            shape = Ferrite.getrefshape(cell)
            ip = IGAInterpolation{shape, getorders(cell)}()
            bernstein_facet_dofs = Ferrite.dirichlet_facetdof_indices(ip)[facetid]
            for (a, row) in pairs(grid.beo[cellid].C)
                if any(j -> abs(row[j]) > 1e-12, bernstein_facet_dofs)
                    push!(nodes, cell.nodes[a])
                end
            end
        else
            union!(nodes, Ferrite.facets(cell)[facetid])
        end
    end
    return nodes
end

#Ferrite's CellIterator does not support BezierGrid directly, so use the underlying grid when 
#filtering node sets for Dirichlet conditions.
function Ferrite.filter_dbc_set(grid::BezierGrid, fhset::AbstractSet{Int}, dbcset::AbstractSet{Int})
    return Ferrite.filter_dbc_set(grid.grid, fhset, dbcset)
end
