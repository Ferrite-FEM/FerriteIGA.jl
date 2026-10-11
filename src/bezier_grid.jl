
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

"""
    BezierGrid(mesh::TSplineMesh)

Grid of the Bézier elements of a T-spline, with one `TSplineCell` per element. Cells with different 
numbers of basis functions have different types, so the grid has mixed cell types. 

The grid contains the nodesets and facetsets `"left"`, `"right"`, `"bottom"` and `"top"` for the 
boundaries of the parametric domain (`mesh.domain`). The nodesets contain the control points whose basis function is 
non-zero on the boundary.
"""
function BezierGrid(mesh::TSplineMesh{sdim,T}) where {sdim,T}
	p = mesh.order
	cells_any = [TSplineCell{RefQuadrilateral,p}(Tuple(ien)) for ien in mesh.IEN]
	CT = Union{unique(typeof.(cells_any))...}
	cells = convert(Vector{CT}, cells_any)
	nodes = [Node(x) for x in mesh.control_points]
	beos = compute_bezier_extraction_operators(mesh)

	(umin, umax), (vmin, vmax) = mesh.domain

	#A basis function is non-zero on a boundary line if the line passes through the interior of its support,
	#or if the line is at the end of the support and the knot there has multiplicity p+1 (open knot vector)
	nonzero_at_start(K, u) = K[1] < u < K[end] || K[1] == K[p+1] == u
	nonzero_at_end(K, u)   = K[1] < u < K[end] || K[2] == K[p+2] == u
	nodesets = Dict{String,OrderedSet{Int}}(name => OrderedSet{Int}() for name in ("left", "right", "bottom", "top"))
	for (A, (Ξ, Π)) in enumerate(mesh.local_knots)
		nonzero_at_start(Ξ, umin) && push!(nodesets["left"], A)
		nonzero_at_end(Ξ, umax)   && push!(nodesets["right"], A)
		nonzero_at_start(Π, vmin) && push!(nodesets["bottom"], A)
		nonzero_at_end(Π, vmax)   && push!(nodesets["top"], A)
	end

	facetsets = Dict{String,OrderedSet{FacetIndex}}(name => OrderedSet{FacetIndex}() for name in ("left", "right", "bottom", "top"))
	for (e, (ξa, ξb, ηa, ηb)) in enumerate(mesh.elements)
		ηa == vmin && push!(facetsets["bottom"], FacetIndex(e, 1))
		ξb == umax && push!(facetsets["right"],  FacetIndex(e, 2))
		ηb == vmax && push!(facetsets["top"],    FacetIndex(e, 3))
		ξa == umin && push!(facetsets["left"],   FacetIndex(e, 4))
	end

	return BezierGrid(cells, nodes, mesh.weights, beos; nodesets, facetsets)
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
	get_nurbs_weights!(w, grid::BezierGrid, cellid)
	get_nurbs_weights!(w, grid::BezierGrid, cell)

Fill `w` with the NURBS weights for the given cell (by id or cell object).
`w` must be at least as long as the number of control points on the cell.
"""
Base.@propagate_inbounds function get_nurbs_weights!(w::Vector, grid::BezierGrid, cellid::Int)
    cell = grid.cells[cellid]
    get_nurbs_weights!(w, grid, cell)
end

Base.@propagate_inbounds function get_nurbs_weights!(w::Vector, grid::BezierGrid, cell::Ferrite.AbstractCell) 
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

"""
	_nbezierpoints(grid::BezierGrid, cellid)

Number of Bézier control points of the cell. This equals the number of nodes (spline control points), 
except for T-spline cells.
"""
function _nbezierpoints(grid::BezierGrid, ic::Int)
	C = grid.beo[ic]
	return C === nothing ? Ferrite.nnodes_per_cell(grid, ic) : length(first(C.C))
end

function Ferrite.getcoordinates!(bc::BezierCoords{dim,T}, grid::BezierGrid, ic::Int) where {dim,T}
	resize_bezier_coord!(bc, _nbezierpoints(grid, ic), Ferrite.nnodes_per_cell(grid, ic))
	get_bezier_coordinates!(bc.xb, bc.wb, bc.x, bc.w, grid, ic)
	bc.beo = grid.beo[ic]
	return bc
end

function Ferrite.getcoordinates(grid::BezierGrid{dim,C,T}, ic::Int) where {dim,C,T}

	n = Ferrite.nnodes_per_cell(grid, ic)
	nb = _nbezierpoints(grid, ic)
	w = zeros(T, n)
	wb = zeros(T, nb)
	xb = zeros(Vec{dim,T}, nb)
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
  These have length `(order+1)^dim`, which differs from the length of `x` and `w` for T-spline cells.

If the cell has no extraction operator, `xb`/`wb` is made as a copy of `x`/`w`.
"""
function get_bezier_coordinates!(xb::AbstractVector{Vec{dim,T}}, 
								 wb::AbstractVector{T}, 
	                             x::AbstractVector{Vec{dim,T}},  
								 w::AbstractVector{T}, 
								 grid::BezierGrid, 
								 ic::Int) where {dim,T}
	Ferrite.getcoordinates!(x, grid.grid, ic)
	get_nurbs_weights!(w, grid, ic)

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
	nb = _nbezierpoints(grid, ic)
	w = zeros(T, n)
	x = zeros(Vec{dim,T}, n)
	wb = zeros(T, nb)
	xb = zeros(Vec{dim,T}, nb)
	
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
