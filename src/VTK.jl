
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
struct VTKIGAFile{VTK<:WriteVTK.DatasetFile}
    vtk::VTK
	cellset::Vector{Int}
end

function VTKIGAFile(filename::String, grid::BezierGrid, cellset=1:getncells(grid); kwargs...)
    vtk = _create_iga_vtk_grid(filename, grid, cellset; kwargs...)
	cellset = sort(collect(copy(cellset)))
    return VTKIGAFile(vtk, cellset)
end

Base.close(vtk::VTKIGAFile) = WriteVTK.vtk_save(vtk.vtk)

function Base.show(io::IO, ::MIME"text/plain", vtk::VTKIGAFile)
    open_str = isopen(vtk.vtk) ? "open" : "closed"
    filename = vtk.vtk.path
    print(io, "VTKFile for the $open_str file \"$(filename)\".")
end

# Makes it possible to use the `do`-block syntax
function VTKIGAFile(f::Function, args...; kwargs...)
    vtk = VTKIGAFile(args...; kwargs...)
    try
        f(vtk)
    finally
        close(vtk)
    end
end

function _create_iga_vtk_grid(filename, grid::BezierGrid{sdim,C,T}, cellset; kwargs...) where {sdim,C,T}
	Ferrite._check_same_celltype(grid, cellset)

	cellset = collect(cellset)
	sort!(cellset)

	cell = grid.cells[first(cellset)]
	reorder = _iga_to_vtkorder(typeof(cell))
	nnodes_per_cell = Ferrite.nnodes(cell)

	#Variables for the vtk file
	cls = WriteVTK.MeshCell[]
	beziercoords = Vec{sdim,T}[]
	weights = T[]
	cellorders = Int[]

	#Variables for the iterations
	bcoords = zeros(Vec{sdim,T}, nnodes_per_cell)
	coords  = zeros(Vec{sdim,T}, nnodes_per_cell)
	wb = zeros(T, nnodes_per_cell)
	w  = zeros(T, nnodes_per_cell)

	offset = 0
	for cellid in cellset
		cell = grid.cells[cellid]

		vtktype = Ferrite.cell_to_vtkcell(typeof(cell))
		for p in getorders(cell)
			push!(cellorders, p)
		end

		get_bezier_coordinates!(bcoords, wb, coords, w, grid, cellid)

		append!(beziercoords, bcoords)
		append!(weights, wb)

		cellnodes = (1:length(cell.nodes)) .+ offset

		push!(cls, WriteVTK.MeshCell(vtktype, collect(cellnodes[reorder])))
		offset += length(cell.nodes)
	end
	
	coords = reshape(reinterpret(T, beziercoords), (sdim, length(beziercoords)))
	vtkfile = WriteVTK.vtk_grid(filename, coords, cls; kwargs...)
	vtkfile["RationalWeights",    WriteVTK.VTKPointData()] = weights
	#vtkfile["HigherOrderDegrees", WriteVTK.VTKCellData()] = reshape(cellorders, 1, length(grid.cells))
	
	return vtkfile
end


function WriteVTK.vtk_point_data(
	vtkfile::WriteVTK.DatasetFile, 
	cpvalues::Vector{<:Union{SymmetricTensor{order,dim,T,M}, 
                             Tensor{order,dimv,T,M}}}, 
	name::AbstractString, 
	grid::BezierGrid{dim,C}) where {order,dimv,dim,C,T,M}

	@assert isconcretetype(C)
	nnodes = Ferrite.nnodes(C)

	data = fill(NaN, M, nnodes*getncells(grid))  # set default value

	nodecount = 0
    for (cellid, cell) in enumerate(grid.cells)
		reorder = Ferrite_to_vtk_order(typeof(first(grid.cells)))
		nodevalues = cpvalues[collect(cell.nodes)]
		_distribute_vtk_point_data!(grid.beo[cellid], data, nodevalues, nodecount)
		nodecount += length(cell.nodes)
    end
    
	WriteVTK.vtk_point_data(vtkfile, data, name)
	
    return vtkfile
end

function Ferrite.write_solution(vtk::VTKIGAFile, dh::DofHandler, a, suffix="")
	for fieldname in Ferrite.getfieldnames(dh)
		data = _evaluate_at_geometry_nodes!(vtk, dh, a, fieldname)
		WriteVTK.vtk_point_data(vtk.vtk, data, string(fieldname, suffix))
	end
end

function Ferrite.write_projection(vtk::VTKIGAFile, proj::L2Projector, vals, name)
    data = _evaluate_at_grid_nodes_iga(proj, vals, #=vtk=# Val(true), vtk.vtk.Npts)::Matrix
    WriteVTK.vtk_point_data(vtk.vtk, data, name; component_names=Ferrite.component_names(eltype(vals)))
    return vtk
end

function _evaluate_at_geometry_nodes!(
	vtk       ::VTKIGAFile,
    dh        ::Ferrite.DofHandler{sdim}, 
    a         ::Vector{T}, 
    fieldname ::Symbol) where {T,sdim}
	# Make sure the field exists (and is not an algebraic variable)
    fieldname ∈ Ferrite.getfieldnames(dh) || error("Field $fieldname not found.")
    # Figure out the return type (scalar or vector)
    field_idx = Ferrite.find_field(dh, fieldname)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    # VTK output of solution field
    n_c = Ferrite.n_components(ip)
    vtk_dim = n_c == 2 ? 3 : n_c # VTK wants vectors padded to 3D
    # Float32 is the smallest float type supported by VTK
    TT = promote_type(T, Float32)
    nviznodes = vtk.vtk.Npts
    data = fill!(Matrix{TT}(undef, vtk_dim, nviznodes), NaN)
    # Loop over the subdofhandlers
    for sdh in dh.subdofhandlers
        # Check if this sdh contains this field, otherwise continue to the next
        field_idx = Ferrite._find_field(sdh, fieldname)
        field_idx === nothing && continue
        # Set up CellValues with the local node coords as quadrature points
        CT = getcelltype(sdh)
        ip = Ferrite.getfieldinterpolation(sdh, field_idx)
        ip_geo = Ferrite.geometric_interpolation(CT)
        local_node_coords = Ferrite.reference_coordinates(ip_geo)
        qr = QuadratureRule{getrefshape(ip)}(zeros(length(local_node_coords)), local_node_coords)
        cv = BezierCellValues(qr, ip, ip_geo^sdim; update_gradients = false, update_hessians = false, update_detJdV = false)
        drange = dof_range(sdh, field_idx)
        # Function barrier
        _nnodes_per_cell = Ferrite.nnodes_per_cell(dh.grid, first(sdh.cellset))
        nodeoffset = _nnodes_per_cell * (findfirst(==(first(sdh.cellset)), vtk.cellset) - 1)
        _evaluate_at_geometry_nodes!(data, sdh, a, cv, drange, sdh.cellset, nodeoffset)
    end
    return data
end

is_scalar_interpolaiton(::BezierCellValues{<:Ferrite.FunctionValues{DiffOrder, <:Ferrite.VectorInterpolation}}) where DiffOrder = false
is_scalar_interpolaiton(::BezierCellValues{<:Ferrite.FunctionValues{DiffOrder, <:Ferrite.ScalarInterpolation}}) where DiffOrder = true
function _evaluate_at_geometry_nodes!(data, sdh, a::Vector{T}, cv, drange, cellset, nodeoffset = 0) where {T, RT}

	dh = sdh.dh

	n_eval_points = Ferrite.getngeobasefunctions(cv)
	ncelldofs = length(drange)
	ue = zeros(eltype(a), ncelldofs)

	dofs = zeros(Int, ncelldofs)
	bcoords = getcoordinates(dh.grid, first(cellset))
    offset = nodeoffset
    for cellid in cellset
        getcoordinates!(bcoords, dh.grid, cellid)

		reinit!(cv, bcoords)

        celldofs!(dofs, sdh, cellid)
		for (i, I) in pairs(drange)
            ue[i] = a[dofs[I]]
        end

        cellnodes = (1:n_eval_points) .+ offset

        for (iqp, nodeid) in pairs(cellnodes)
            val = function_value(cv, iqp, ue)
			if data isa Matrix # VTK
                data[1:length(val), nodeid] .= val
                data[(length(val)+1):end, nodeid] .= 0 # purge the NaN
            else
                data[nodeid] = val
            end
        end

        offset += n_eval_points
    end

    return data
end

function _evaluate_at_grid_nodes_iga(
        proj::L2Projector, vals::AbstractVector{S}, ::Val{vtk}, nviznodes
    ) where {order, dim, T, M, S <: Union{Tensor{order, dim, T, M}, SymmetricTensor{order, dim, T, M}}, vtk}
    dh = proj.dh
    # The internal dofhandler in the projector is a scalar field, but the values in vals
    # can be any tensor field, however, the number of dofs should always match the length of vals
    @assert ndofs(dh) == length(vals)
    if vtk
        nout = S <: Vec{2} ? 3 : M # Pad 2D Vec to 3D
        data = fill(T(NaN), nout, nviznodes)
    else
        data = fill(T(NaN) * zero(S), nviznodes)
    end
    for sdh in dh.subdofhandlers
        ip = only(sdh.field_interpolations)
        gip = geometric_interpolation(getcelltype(sdh))
        RefShape = getrefshape(ip)
        local_node_coords = Ferrite.reference_coordinates(gip)
        qr = QuadratureRule{RefShape}(zeros(length(local_node_coords)), local_node_coords)
        cv = BezierCellValues(qr, ip, gip; update_detJdV = false, update_gradients = false)
        _evaluate_at_grid_nodes_iga!(data, cv, sdh, vals)
    end
    return data
end


function _evaluate_at_grid_nodes_iga!(data, cv::BezierCellValues, sdh::SubDofHandler, u::AbstractVector{S}) where {S}
    ue = zeros(S, getnbasefunctions(cv))
	bcoords = getcoordinates(sdh.dh.grid, first(sdh.cellset))
	dofs = zeros(Int, ndofs_per_cell(sdh))
	offset = 0
	nnodes = length(bcoords.x)
    for cellid in sdh.cellset
        getcoordinates!(bcoords, sdh.dh.grid, cellid)
		celldofs!(dofs, sdh, cellid)
		reinit!(cv, bcoords)
        for (i, I) in pairs(dofs)
            ue[i] = u[I]
        end

		cellnodes = (1:nnodes) .+ offset

        for (qp, nodeid) in pairs(cellnodes)
            # Loop manually over the shape functions since function_value
            # doesn't like scalar base functions with tensor dofs
            val = zero(S)
            for i in 1:getnbasefunctions(cv)
                val += shape_value(cv, qp, i) * ue[i]
            end

            if data isa Matrix # VTK
                dataview = @view data[:, nodeid]
                fill!(dataview, 0) # purge the NaN
                Ferrite.toparaview!(dataview, val)
            else
                data[nodeid] = val
            end
        end

		offset += nnodes
    end
    return data
end
