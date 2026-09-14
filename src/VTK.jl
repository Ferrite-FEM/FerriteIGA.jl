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
    global function _iga_to_vtkorder(celltype::Type{<:BezierCell{shape,order,N}}) where {order,shape,N}
        get!(cache, celltype) do 
            if shape == RefHexahedron
                igaorder = _bernstein_ordering(celltype)
                vtkorder = _vtk_ordering(celltype)

                return [findfirst(ivtk-> ivtk == iiga, vtkorder) for iiga in igaorder]
            else
                return collect(1:N)
            end
        end
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

function _create_iga_vtk_grid(filename, grid::BezierGrid{sdim,C,T}; kwargs...) where {sdim,C,T}
    ncells = Ferrite.getncells(grid)
    
    cls = WriteVTK.MeshCell[]
    beziercoords = Vec{sdim,T}[]
    weights = T[]
    cellorders = Int[]

    # Assuming uniform grid cell types for arrays sizing
    sample_cell = first(grid.cells)
    nnodes_per_cell = Ferrite.nnodes(sample_cell)

    bcoords = zeros(Vec{sdim,T}, nnodes_per_cell)
    coords  = zeros(Vec{sdim,T}, nnodes_per_cell)
    wb = zeros(T, nnodes_per_cell)
    w  = zeros(T, nnodes_per_cell)

    cellnodes = Vector{UnitRange{Int}}(undef, ncells)
    node_mapping = Int[]

    offset = 0
    for cellid in 1:ncells
        cell = grid.cells[cellid]
        vtktype = Ferrite.cell_to_vtkcell(typeof(cell))
        reorder = _iga_to_vtkorder(typeof(cell))

        for p in getorders(cell)
            push!(cellorders, p)
        end

        get_bezier_coordinates!(bcoords, wb, coords, w, grid, cellid)

        append!(beziercoords, bcoords)
        append!(weights, wb)
        append!(node_mapping, cell.nodes)

        cnodes = (1:length(cell.nodes)) .+ offset
        cellnodes[cellid] = cnodes

        push!(cls, WriteVTK.MeshCell(vtktype, collect(cnodes[reorder])))
        offset += length(cell.nodes)
    end
    
    coords_matrix = reshape(reinterpret(T, beziercoords), (sdim, length(beziercoords)))
    vtkfile = WriteVTK.vtk_grid(filename, coords_matrix, cls; kwargs...)
    vtkfile["RationalWeights", WriteVTK.VTKPointData()] = weights
    
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
        cv = BezierCellValues(qr, ip, ip_geo; update_detJdV = false, update_gradients = false)
        _evaluate_at_discontinuous_vtkgrid_nodes_iga!(data, sdh, u, cv, cellnodes, drange)
    end
    return data
end

function _evaluate_at_discontinuous_vtkgrid_nodes_iga!(
    data::Matrix, sdh::Ferrite.SubDofHandler,
    u::AbstractVector{S}, cv::BezierCellValues, cellnodes, drange
) where {S}
    maybe_is_projector = !(S <: AbstractFloat)
    n_base = length(drange)
    ue = zeros(S, n_base)
    dofs = zeros(Int, Ferrite.ndofs_per_cell(sdh))
    bcoords = getcoordinates(sdh.dh.grid, first(sdh.cellset))
    
    for cellid in sdh.cellset
        getcoordinates!(bcoords, sdh.dh.grid, cellid)
        celldofs!(dofs, sdh, cellid)
        reinit!(cv, bcoords)
        
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