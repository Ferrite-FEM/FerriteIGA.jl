export BezierGrid, getweights, getweights!, get_extraction_operator, get_bezier_coordinates, get_bezier_coordinates!, get_nurbs_coordinates

"""
    BezierGrid(mesh::NURBSMesh)
    BezierGrid(grid::Grid)

Ferrite grid with NURBS weights and an extraction operator on each cell.
`BezierGrid(mesh)` builds this from a NURBS patch.
`BezierGrid(grid)` wraps an ordinary Ferrite grid.
"""
struct BezierGrid{dim,C<:Ferrite.AbstractCell,T<:Real} <: Ferrite.AbstractGrid{dim}
	grid    ::Ferrite.Grid{dim,C,T}
	weights ::Vector{Float64}
	beo     ::Vector{BezierExtractionOperator{Float64}}

	#Alias from grid
	cells::Vector{C}
    nodes::Vector{Node{dim,T}}
    cellsets::Dict{String,OrderedSet{Int}}
    nodesets::Dict{String,OrderedSet{Int}}
    facetsets::Dict{String,OrderedSet{FacetIndex}}
    vertexsets::Dict{String,OrderedSet{VertexIndex}}

	function BezierGrid(grid::Ferrite.Grid{dim,C,T}, weights::Vector{T}, beo::Vector{BezierExtractionOperator{T}}) where {dim,C,T}
		return new{dim,C,T}(
			grid, weights, beo,
			grid.cells, grid.nodes, grid.cellsets, grid.nodesets, grid.facetsets, grid.vertexsets, 
		)
	end
end

function BezierGrid(cells::Vector{C},
		nodes::Vector{Ferrite.Node{dim,T}},
		weights::AbstractVector{T},
		extraction_operator::AbstractVector{BezierExtractionOperator{T}}; 
		cellsets::Dict{String,OrderedSet{Int}}             = Dict{String,OrderedSet{Int}}(),
		nodesets::Dict{String,OrderedSet{Int}}             = Dict{String,OrderedSet{Int}}(),
		facetsets::Dict{String,OrderedSet{FacetIndex}}      = Dict{String,OrderedSet{FacetIndex}}(),
		vertexsets::Dict{String,OrderedSet{VertexIndex}}   = Dict{String,OrderedSet{VertexIndex}}()) where {dim,C,T}

	
	grid = Ferrite.Grid(cells, nodes; nodesets, cellsets, facetsets, vertexsets)

	return BezierGrid(grid, weights, extraction_operator)
end

function BezierGrid(mesh::NURBSMesh{pdim,sdim}) where {pdim,sdim}

	@assert allequal(mesh.orders)
	order = first(mesh.orders)
	N = size(mesh.IEN, 1)
    CellType = BezierCell{RefHypercube{pdim},order}
    ordering = _bernstein_ordering(CellType)
	
	cells = [CellType(Tuple(mesh.IEN[ordering,ie])) for ie in 1:getncells(mesh)]
	nodes = [Node(x)                                for x  in mesh.control_points]

    C, nbe = compute_bezier_extraction_operators(mesh.orders, mesh.knot_vectors)

	@assert nbe == length(cells)

	Cvec = bezier_extraction_to_vectors(C)

	return BezierGrid(cells, nodes, mesh.weights, Cvec)
end

function BezierGrid(grid::Ferrite.Grid{dim,C,T}) where {dim,C,T}
	weights = ones(Float64, getnnodes(grid))

	extraction_operator = BezierExtractionOperator{T}[]
	for cellid in 1:getncells(grid)
		nnodes = length(grid.cells[cellid].nodes)
		beo = diagonal_beo(nnodes)
		push!(extraction_operator, beo)
	end

	return BezierGrid(grid, weights, extraction_operator)
end

function Base.show(io::IO, ::MIME"text/plain", grid::BezierGrid)
    print(io, "$(typeof(grid)) with $(getncells(grid)) ")
    if isconcretetype(eltype(grid.cells))
        typestrs = [repr(eltype(grid.cells))]
    else
        typestrs = sort!(repr.(Set(typeof(x) for x in grid.cells)))
    end
    join(io, typestrs, '/')
    print(io, " cells and $(getnnodes(grid)) nodes (contorl points)")
end

"""
	getweights!(w::Vector{T}, grid::BezierGrid, cellid::Int) where {T} 

Returns the weights (for the nurbs interpolation) for cell with id `cellid`.
"""
Base.@propagate_inbounds function getweights!(w::Vector, grid::BezierGrid, cellid::Int)
    cell = grid.cells[cellid]
    getweights!(w, grid, cell)
end

Base.@propagate_inbounds function getweights!(w::Vector, grid::BezierGrid, cell::Ferrite.AbstractCell) 
	node_ids = Ferrite.get_node_ids(cell)
	nnodes = length(node_ids)
	@boundscheck checkbounds(Bool, w, 1:nnodes)
    @inbounds for i in 1:nnodes
        w[i] = grid.weights[node_ids[i]]
    end
end

function Ferrite.getweights(grid::BezierGrid, ic::Int)
	nodeids = collect(grid.cells[ic].nodes)
	return grid.weights[nodeids]
end

#Note: This is needed in order to get some stuff to work nicely in Ferrite.
Ferrite.get_coordinate_type(::BezierGrid{dim,C,T}) where {dim,C,T} = Vec{dim,T} 

function Ferrite.getcoordinates!(bc::BezierCoords{dim,T}, grid::BezierGrid, ic::Int) where {dim,T}
	get_bezier_coordinates!(bc.xb, bc.wb, bc.x, bc.w, grid, ic)
	bc.beo[] = grid.beo[ic]
	return bc
end

function Ferrite.getcoordinates(grid::BezierGrid{dim,C,T}, ic::Int) where {dim,C,T}

	n = Ferrite.nnodes_per_cell(grid, ic)
	w = zeros(T, n)
	wb = zeros(T, n)
	xb = zeros(Vec{dim,T}, n)
	x = zeros(Vec{dim,T}, n)
	
	bc = BezierCoords(xb, wb, x, w, Ref(grid.beo[ic]))
	getcoordinates!(bc,grid,ic)

	return bc
end

function get_bezier_coordinates!(xb::AbstractVector{Vec{dim,T}}, 
								 wb::AbstractVector{T}, 
	                             x::AbstractVector{Vec{dim,T}},  
								 w::AbstractVector{T}, 
								 grid::BezierGrid, 
								 ic::Int) where {dim,T}
    n = length(xb)
	C = grid.beo[ic]

	Ferrite.getcoordinates!(x, grid.grid, ic)
	getweights!(w, grid, ic)

	for i in 1:n
		xb[i] = zero(Vec{dim,T})
		wb[i] = zero(T)
	end

	for i in 1:n
		c_row = C[i]
		_w = w[i]
		_x = _w*x[i]
		for (j, nz_ind) in enumerate(c_row.nzind)                
			xb[nz_ind] += c_row.nzval[j] * _x
			wb[nz_ind] += c_row.nzval[j] * _w
		end
	end
	xb ./= wb

	return nothing
end

function get_bezier_coordinates(grid::BezierGrid{dim,C,T}, ic::Int) where {dim,C,T}

	n = Ferrite.nnodes_per_cell(grid, ic)
	w = zeros(T, n)
	x = zeros(Vec{dim,T}, n)
	wb = zeros(T, n)
	xb = zeros(Vec{dim,T}, n)
	
	get_bezier_coordinates!(xb, wb, x, w, grid, ic)
	return xb, wb, x, w
end

function get_nurbs_coordinates(grid::BezierGrid{dim,C,T}, cell::Int) where {dim,C,T}
    nodeidx = grid.cells[cell].nodes
    return [grid.nodes[i].x for i in nodeidx]::Vector{Vec{dim,T}}
end

"""
	get_extraction_operator(grid::BezierGrid, cellid)

Bézier extraction operator for cell `cellid`.
"""
function get_extraction_operator(grid::BezierGrid, cellid::Int)
	return grid.beo[cellid]
end
