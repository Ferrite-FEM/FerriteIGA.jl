
"""
    BezierGrid(mesh::NURBSMesh)
    BezierGrid(grid::Grid)

Ferrite grid with NURBS weights and an extraction operator on each cell.
`BezierGrid(mesh)` builds this from a NURBS patch.
`BezierGrid(grid)` wraps an ordinary Ferrite grid.
"""
struct BezierGrid{dim,C<:Ferrite.AbstractCell,T<:Real} <: Ferrite.AbstractGrid{dim}
	grid    ::Ferrite.Grid{dim,C,T}
	weights ::Vector{T}
	beo     ::Vector{Optional{BezierExtractionOperator{T}}}

	#Alias from grid
	cells::Vector{C}
    nodes::Vector{Node{dim,T}}
    cellsets::Dict{String,OrderedSet{Int}}
    nodesets::Dict{String,OrderedSet{Int}}
    facetsets::Dict{String,OrderedSet{FacetIndex}}
    vertexsets::Dict{String,OrderedSet{VertexIndex}}

	function BezierGrid(grid::Ferrite.Grid{dim,C,T}, weights::Vector{T}, beo::Vector{Optional{BezierExtractionOperator{T}}}) where {dim,C,T}
		return new{dim,C,T}(
			grid, weights, beo,
			grid.cells, grid.nodes, grid.cellsets, grid.nodesets, grid.facetsets, grid.vertexsets, 
		)
	end
end

function BezierGrid(cells::Vector{C},
		nodes::Vector{Ferrite.Node{dim,T}},
		weights::AbstractVector{T},
		extraction_operators::Union{Vector{BezierExtractionOperator{T}}, Vector{Optional{BezierExtractionOperator{T}}}}; 
		cellsets::Dict{String,OrderedSet{Int}}             = Dict{String,OrderedSet{Int}}(),
		nodesets::Dict{String,OrderedSet{Int}}             = Dict{String,OrderedSet{Int}}(),
		facetsets::Dict{String,OrderedSet{FacetIndex}}      = Dict{String,OrderedSet{FacetIndex}}(),
		vertexsets::Dict{String,OrderedSet{VertexIndex}}   = Dict{String,OrderedSet{VertexIndex}}()) where {dim,C,T}

	extraction_operators = convert(Vector{Optional{BezierExtractionOperator{T}}}, extraction_operators)
	grid = Ferrite.Grid(cells, nodes; nodesets, cellsets, facetsets, vertexsets)

	return BezierGrid(grid, weights, extraction_operators)
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
	weights = ones(T, getnnodes(grid))
	extraction_operator = Optional{BezierExtractionOperator{T}}[nothing for _ in 1:getncells(grid)]

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
	getweights!(w, grid::BezierGrid, cellid)
	getweights!(w, grid::BezierGrid, cell)

Fill `w` with the NURBS weights for the given cell (by id or cell object).
`w` must be at least as long as the number of control points on the cell.
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

"""
	get_nurbs_weights(grid::BezierGrid, cellid)

Return a new vector of NURBS weights for cell `cellid`.
"""
function get_nurbs_weights(grid::BezierGrid, ic::Int)
	nodeids = collect(grid.cells[ic].nodes)
	return grid.weights[nodeids]
end

#Note: This is needed in order to get some stuff to work nicely in Ferrite.
Ferrite.get_coordinate_type(::BezierGrid{dim,C,T}) where {dim,C,T} = Vec{dim,T} 

function Ferrite.getcoordinates!(bc::BezierCoords{dim,T}, grid::BezierGrid, ic::Int) where {dim,T}
	get_bezier_coordinates!(bc.xb, bc.wb, bc.x, bc.w, grid, ic)
	bc.beo = grid.beo[ic]
	return bc
end

function Ferrite.getcoordinates(grid::BezierGrid{dim,C,T}, ic::Int) where {dim,C,T}

	n = Ferrite.nnodes_per_cell(grid, ic)
	w = zeros(T, n)
	wb = zeros(T, n)
	xb = zeros(Vec{dim,T}, n)
	x = zeros(Vec{dim,T}, n)
	
	bc = BezierCoords(xb, wb, x, w, grid.beo[ic])
	getcoordinates!(bc,grid,ic)

	return bc
end

"""
	get_bezier_coordinates!(xb, wb, x, w, grid::BezierGrid, cellid)

In-place fill of NURBS and Bézier data for cell `cellid`.

- `x`, `w` are the NURBS control points and weights respectively.
- `xb`, `wb` are the above transformed in Bernstein basis from the Bezier extraction operator.

If the cell has no extraction operator, `xb`/`wb` is made as a copy of `x`/`w`.
"""
function get_bezier_coordinates!(xb::AbstractVector{Vec{dim,T}}, 
								 wb::AbstractVector{T}, 
	                             x::AbstractVector{Vec{dim,T}},  
								 w::AbstractVector{T}, 
								 grid::BezierGrid, 
								 ic::Int) where {dim,T}
    n = length(xb)
	
	Ferrite.getcoordinates!(x, grid.grid, ic)
	getweights!(w, grid, ic)

	C = grid.beo[ic]
	if C === nothing #This is not an IGA cell
		xb .= x
		wb .= w
	else
		transform_coords!(xb, wb, C, x, w)
	end

	return nothing
end

"""
	get_bezier_coordinates(grid::BezierGrid, cellid)

Allocate and return `(xb, wb, x, w)` for cell `cellid`:
Bézier points/weights and NURBS points/weights.
See [`get_bezier_coordinates!`](@ref).
"""
function get_bezier_coordinates(grid::BezierGrid{dim,C,T}, ic::Int) where {dim,C,T}

	n = Ferrite.nnodes_per_cell(grid, ic)
	w = zeros(T, n)
	x = zeros(Vec{dim,T}, n)
	wb = zeros(T, n)
	xb = zeros(Vec{dim,T}, n)
	
	get_bezier_coordinates!(xb, wb, x, w, grid, ic)
	return xb, wb, x, w
end

"""
	get_nurbs_coordinates(grid::BezierGrid, cellid)

Return the NURBS control-point coordinates for cell `cellid`
"""
function get_nurbs_coordinates(grid::BezierGrid{dim,C,T}, cell::Int) where {dim,C,T}
    nodeidx = grid.cells[cell].nodes
    return [grid.nodes[i].x for i in nodeidx]::Vector{Vec{dim,T}}
end

"""
	get_extraction_operator(grid::BezierGrid, cellid)

Bézier extraction operator ``C^e`` for cell `cellid`, or `nothing`
if the cell is not an IGA cell.
"""
function get_extraction_operator(grid::BezierGrid, cellid::Int)
	return grid.beo[cellid]
end
