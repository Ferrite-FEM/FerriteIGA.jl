
function Ferrite.cell_to_vtkcell(::Type{<:BezierCell{RefHexahedron,order}}) where {order}
    return Ferrite.VTKCellTypes.VTK_BEZIER_HEXAHEDRON
end
function Ferrite.cell_to_vtkcell(::Type{<:BezierCell{RefQuadrilateral,order}}) where {order}
    return Ferrite.VTKCellTypes.VTK_BEZIER_QUADRILATERAL
end
function Ferrite.cell_to_vtkcell(::Type{<:BezierCell{RefLine,order}}) where {order}
    return Ferrite.VTKCellTypes.VTK_BEZIER_CURVE
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
    open_str = WriteVTK.isopen(vtk.vtk) ? "open" : "closed"
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
	cls = MeshCell[]
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
                             Tensor{order,dimv,T,M}, 
                             SVector{M,T}}}, 
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
    
	vtk_point_data(vtkfile, data, name)
	
    return vtkfile
end

function Ferrite.write_solution(vtk::VTKIGAFile, dh::DofHandler, a, suffix="")
	for fieldname in Ferrite.getfieldnames(dh)
		data = _evaluate_at_geometry_nodes!(vtk, dh, a, fieldname)
		vtk_point_data(vtk.vtk, data, string(fieldname, suffix))
	end
end

function Ferrite.write_projection(vtk::VTKIGAFile, proj::L2Projector, vals, name)
    data = Ferrite._evaluate_at_grid_nodes(proj, vals, #=vtk=# Val(true))::Matrix
    @assert size(data, 2) == getnnodes(Ferrite.get_grid(proj.dh))
    vtk_point_data(vtk.vtk, data, name; component_names=Ferrite.component_names(eltype(vals)))
    return vtk
end

function _evaluate_at_geometry_nodes!(
	vtk       ::VTKIGAFile,
    dh        ::Ferrite.AbstractDofHandler, 
    a         ::Vector{T}, 
    fieldname ::Symbol) where T

	Ferrite._check_same_celltype(dh.grid, vtk.cellset)
    fieldname ∈ Ferrite.getfieldnames(dh) || error("Field $fieldname not found in the dofhandler.")

    # 
    local sdh
    for _sdh in dh.subdofhandlers
        if first(vtk.cellset) in _sdh.cellset
            sdh = _sdh
            @assert Set{Int}(vtk.cellset) == _sdh.cellset
            break
        end
    end

	#
	CT = getcelltype(dh.grid, first(sdh.cellset))

    field_dim = Ferrite.n_components(sdh, fieldname)
    ncomponents = field_dim == 2 ? 3 : field_dim
    
    nviznodes = vtk.vtk.Npts
    data = fill(Float64(NaN), ncomponents, nviznodes) 

    field_idx = Ferrite.find_field(sdh, fieldname)
    field_idx === nothing && error("The field $fieldname does not exist in the subdofhandler")
	ip_geo = Ferrite.geometric_interpolation(CT)
    ip     = Ferrite.getfieldinterpolation(sdh, field_idx)
    drange = Ferrite.dof_range(sdh, fieldname)
	shape  = Ferrite.getrefshape(ip_geo)

	#
	local_node_coords = Ferrite.reference_coordinates(ip_geo)
	qr = QuadratureRule{shape}(zeros(length(local_node_coords)), local_node_coords)
	ip = Ferrite.getfieldinterpolation(sdh, field_idx)
	
	# TODO: Remove this hack when embedding works...
	RT = ip isa ScalarInterpolation ? T : Vec{Ferrite.n_components(ip),T}
	if ip isa VectorizedInterpolation
		cv = CellValues(qr, ip.ip, ip_geo)
	else
		cv = CellValues(qr, ip, ip_geo)
	end

    _evaluate_at_geometry_nodes!(data, sdh, a, cv, drange, vtk.cellset, RT)
        
    return data
end

is_scalar_interpolaiton(::BezierCellValues{<:Ferrite.FunctionValues{DiffOrder, <:Ferrite.VectorInterpolation}}) where DiffOrder = false
is_scalar_interpolaiton(::BezierCellValues{<:Ferrite.FunctionValues{DiffOrder, <:Ferrite.ScalarInterpolation}}) where DiffOrder = true
function _evaluate_at_geometry_nodes!(data, sdh, a::Vector{T}, cv, drange, cellset, ::Type{RT}) where {T, RT}

	dh = sdh.dh
	grid = dh.grid

	n_eval_points = Ferrite.getngeobasefunctions(cv)
	ncelldofs = length(drange)
	ue = zeros(eltype(a), ncelldofs)

	# TODO: Remove this hack when embedding works...
    if RT <: Vec && is_scalar_interpolaiton(cv)
        uer = reinterpret(RT, ue)
    else
        uer = ue
    end

	dofs = zeros(Int, ncelldofs)
	bcoords = getcoordinates(dh.grid, first(cellset))
    offset = 0
    for cellid in cellset
        getcoordinates!(bcoords, dh.grid, cellid)

		reinit_values!(cv, bcoords)

        celldofs!(dofs, sdh, cellid)
		for (i, I) in pairs(drange)
            ue[i] = a[dofs[I]]
        end

        cellnodes = (1:n_eval_points) .+ offset

        for (iqp, nodeid) in pairs(cellnodes)
            val = function_value(cv, iqp, uer)
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
