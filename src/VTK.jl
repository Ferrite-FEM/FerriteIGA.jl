function Ferrite.cell_to_vtkcell(::Type{<:BezierCell{RefHexahedron,order}}) where {order}
    return VTKCellTypes.VTK_BEZIER_HEXAHEDRON
end
function Ferrite.cell_to_vtkcell(::Type{<:BezierCell{RefQuadrilateral,order}}) where {order}
    return VTKCellTypes.VTK_BEZIER_QUADRILATERAL
end
function Ferrite.cell_to_vtkcell(::Type{<:BezierCell{RefLine,order}}) where {order}
    return VTKCellTypes.VTK_BEZIER_CURVE
end

# Store the Ferrite to vtk order in a cache for specific cell type
let cache = Dict{Type{<:BezierCell}, Vector{Int}}()
    global function Ferrite.nodes_to_vtkorder(cell::BezierCell{shape,order,N}) where {order,shape,N}
        celltype = typeof(cell)
        reorder = get!(cache, celltype) do 
            if shape == RefHexahedron
                igaorder = _bernstein_ordering(celltype)
                vtkorder = _vtk_ordering(celltype)

                [findfirst(ivtk-> ivtk == iiga, vtkorder) for iiga in igaorder]
            else
                collect(1:N)
            end
        end
        return cell.nodes[reorder]
    end
end

"""
    VTKIGAFile(filename, grid::BezierGrid)

Write a `BezierGrid` to a `.vtu` file.
"""
struct VTKIGAFile{VTK <: WriteVTK.DatasetFile}
    vtk::VTK
    cellnodes::Union{Vector{UnitRange{Int}}, Nothing}
    node_mapping::Union{Vector{Int}, Nothing}
end

function VTKIGAFile(filename::String, grid::BezierGrid; kwargs...)
    vtk, cellnodes, node_mapping = _create_iga_vtk_grid(filename, grid; kwargs...)
    return VTKIGAFile(vtk, cellnodes, node_mapping)
end

# Makes it possible to use the `do`-block syntax
function VTKIGAFile(f::Function, args...; kwargs...)
    vtk = VTKIGAFile(args...; kwargs...)
    try
        f(vtk)
    finally
        close(vtk)
    end
    return vtk
end

Base.close(vtk::VTKIGAFile) = (WriteVTK.vtk_save(vtk.vtk); return vtk)

function Base.show(io::IO, ::MIME"text/plain", vtk::VTKIGAFile)
    open_str = isopen(vtk.vtk) ? "open" : "closed"
    filename = vtk.vtk.path
    print(io, "VTKIGAFile for the $open_str file \"$(filename)\".")
    return nothing
end

function _create_iga_vtk_grid(filename, grid::BezierGrid{dim,C,T}; kwargs...) where {dim,C,T}
    cls = Vector{WriteVTK.MeshCell}(undef, getncells(grid))
    cellnodes = Vector{UnitRange{Int}}(undef, getncells(grid))
    ncoords = sum(Ferrite.nnodes, getcells(grid))
    coords = zeros(T, dim, ncoords)
    weights = zeros(T, dim, ncoords)
    node_mapping = zeros(Int, ncoords)
    icoord = 0

    nnodes_per_cell = 27; #Some init guess
    cell_bezier_coords = zeros(Vec{dim,T}, nnodes_per_cell)
    cell_nurbs_coords  = zeros(Vec{dim,T}, nnodes_per_cell)
    cell_bezier_weights = zeros(T, nnodes_per_cell)
    cell_nurbs_weights  = zeros(T, nnodes_per_cell)
    for cellid in 1:getncells(grid)
        cell = getcells(grid, cellid)
        CT = typeof(cell)
        vtk_celltype = Ferrite.cell_to_vtkcell(CT)
        n = length(cell.nodes)
        resize!(cell_bezier_coords, n); resize!(cell_nurbs_coords, n); resize!(cell_bezier_weights, n); resize!(cell_nurbs_weights, n);
        get_bezier_coordinates!(cell_bezier_coords, cell_bezier_weights, cell_nurbs_coords, cell_nurbs_weights, grid, cellid)
        cellnodes[cellid] = (1:n) .+ icoord
        let icoord = icoord
            vtk_cellnodes = Ferrite.nodes_to_vtkorder(CT((ntuple(i -> i + icoord, n))))
            cls[cellid] = WriteVTK.MeshCell(vtk_celltype, vtk_cellnodes)
        end
        for (x, node_idx, wb) in zip(cell_bezier_coords, cell.nodes, cell_bezier_weights)
            icoord += 1
            coords[:, icoord] = x
            node_mapping[icoord] = node_idx
            weights[icoord] = wb
        end
    end
    
    vtkfile = WriteVTK.vtk_grid(filename, coords, cls; kwargs...)
    vtkfile["RationalWeights", WriteVTK.VTKPointData()] = weights
    
    # If we in the future allow for different orders in each direction,
    # we have to add the cell orders to the vtk file
    #for p in getorders(cell)
    #    push!(cellorders, p)
    #end
    # ...and
    #vtk["HigherOrderDegrees", VTKCellData()] = [2 3 0]

    return vtkfile, cellnodes, node_mapping
end

function Ferrite.write_solution(vtk::VTKIGAFile, dh::Ferrite.AbstractDofHandler, a, suffix="")
    for fieldname in Ferrite.getfieldnames(dh)
        data = evaluate_at_discontinuous_vtkgrid_nodes_iga(dh, a, vtk.cellnodes, fieldname)
        Ferrite._vtk_write_node_data(vtk.vtk, data, string(fieldname, suffix))
    end
    return vtk
end

function Ferrite.write_projection(vtk::VTKIGAFile, proj::Ferrite.L2Projector, vals, name)
    data = evaluate_at_discontinuous_vtkgrid_nodes_iga(proj.dh, vals, vtk.cellnodes, only(Ferrite.getfieldnames(proj.dh)))
    Ferrite._vtk_write_node_data(vtk.vtk, data, name; component_names=Ferrite.component_names(eltype(vals)))
    return vtk
end

function Ferrite.write_node_data(vtk::VTKIGAFile, nodedata, name)
    data = Ferrite._map_to_discontinuous_nodes(vtk.node_mapping, nodedata)
    Ferrite._vtk_write_node_data(vtk.vtk, data, name)
    return vtk
end

function Ferrite.write_cell_data(vtk::VTKIGAFile, celldata, name)
    WriteVTK.vtk_cell_data(vtk.vtk, celldata, name)
    return vtk
end

function evaluate_at_discontinuous_vtkgrid_nodes_iga(
    dh::Ferrite.AbstractDofHandler, 
    u::AbstractVector{S}, 
    cellnodes, 
    fieldname::Symbol
) where {S}

    maybe_is_projector = !(S <: AbstractFloat)
    # Make sure the field exists
    fieldname ∈ Ferrite.getfieldnames(dh) || error("Field $fieldname not found.")
    # Figure out the return type (scalar or vector)
    field_idx = Ferrite.find_field(dh, fieldname)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    n_c = maybe_is_projector ? length(zero(S)) : Ferrite.n_components(ip)
    vtk_dim = n_c == 2 ? 3 : n_c
    
    TT = promote_type(eltype(S), Float32)
    n_vtk_nodes = maximum(maximum, cellnodes)
    
    data = fill!(Matrix{TT}(undef, vtk_dim, n_vtk_nodes), NaN)
    for sdh in dh.subdofhandlers
        field_idx = Ferrite._find_field(sdh, fieldname)
        field_idx === nothing && continue
        ip = Ferrite.getfieldinterpolation(sdh, field_idx)
        drange = Ferrite.dof_range(sdh, field_idx)
        CT = Ferrite.getcelltype(sdh)
        ip_geo = Ferrite.geometric_interpolation(CT)
        local_node_coords = Ferrite.reference_coordinates(ip_geo)
        qr = QuadratureRule{Ferrite.getrefshape(ip)}(zeros(length(local_node_coords)), local_node_coords)
        # The call to CellValues below will create either a CellValue or a BezierCellValue depending on interpolation type
        # TODO: It would be better to have an interface for creating cell values, e.g. create_cellvalue(qr, ip, ip_geo; ...)
        cv = CellValues(qr, ip, ip_geo; update_detJdV = false, update_gradients = false)
        _evaluate_at_discontinuous_vtkgrid_nodes_iga!(data, sdh, u, cv, cellnodes, drange)
    end
    return data
end

function _evaluate_at_discontinuous_vtkgrid_nodes_iga!(
    data::Matrix, sdh::Ferrite.SubDofHandler,
    u::AbstractVector{S}, cv::Ferrite.AbstractCellValues, cellnodes, drange
) where {S}
    maybe_is_projector = !(S <: AbstractFloat)
    n_base = length(drange)
    ue = zeros(S, n_base)
    dofs = zeros(Int, Ferrite.ndofs_per_cell(sdh))
    bcoords = getcoordinates(sdh.dh.grid, first(sdh.cellset))

    _reinit!(cv::Ferrite.AbstractCellValues, cell::Optional{Ferrite.AbstractCell}, coords::BezierCoords) = reinit!(cv, cell, coords.x)
    _reinit!(cv::BezierCellValues, cell::Optional{Ferrite.AbstractCell}, coords::BezierCoords) = reinit!(cv, cell, coords)

    for cellid in sdh.cellset
        cell = getcells(sdh.dh.grid, cellid)
        getcoordinates!(bcoords, sdh.dh.grid, cellid)
        celldofs!(dofs, sdh, cellid)
        _reinit!(cv, cell, bcoords)
        
        for (i, I) in pairs(drange)
            ue[i] = u[dofs[I]]
        end

        for (qp, nodeid) in pairs(cellnodes[cellid])
            #If the functions is called from a L2Projector the type should be S
            val = maybe_is_projector ? zero(S) : zero(Ferrite.shape_value_type(cv))
            @assert getnbasefunctions(cv) == length(ue)
            for i in 1:n_base
                val += shape_value(cv, qp, i) * ue[i]
            end
            
            dataview = @view data[:, nodeid]
            fill!(dataview, 0) # purge the NaN
            Ferrite.toparaview!(dataview, val)
        end
    end
    return data
end