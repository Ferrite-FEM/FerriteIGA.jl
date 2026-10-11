"""
	TSplineMesh(knots, order, cells; control_points, weights, domain)

A (rational) T-spline with parametric dimension 2 or 3 and polynomial degree `order` in all directions.

The T-mesh is given in *index space*:
- `knots::NTuple{dim,Vector{T}}` are the knot values of the index lines (planes in 3D) `i = 1:nᵈ` in each direction `d`.
  Consecutive knots may be equal (zero knot intervals). An open (interpolatory) boundary is obtained by repeating
  the boundary value `order÷2 + 1` times, e.g. `[0,0,1,2,3,3]` for cubic splines.
- `cells` are the boxes of the T-mesh: rectangles `(i0, i1, j0, j1)` in 2D and boxes `(i0, i1, j0, j1, k0, k1)` in 3D.
  They must tile the index domain.

One basis function is associated with each *anchor*: the vertices of the T-mesh for odd degrees, and the
centers of the cells for even degrees. The anchors are sorted with the x-index running fastest and the last index slowest, and
`control_points` and `weights` must be given in this order. If `control_points` are not given, the Greville points
of the basis functions are used, which gives the identity map from the parametric domain (see [`greville_points`](@ref)).

The parametric `domain = ((ξ0, ξ1), (η0, η1)[, (ζ0, ζ1)])` defaults to the range of the knots. A smaller domain can be used
to let the T-mesh extend outside of the domain (a *frame*), e.g. with periodically extended knots, which gives
basis functions that are not interpolatory at the boundary. The bounds of the domain must be lines (planes) of the T-mesh.
Bézier elements are only created inside the domain, and only the basis functions which are non-zero in the domain
are kept (and numbered).

Use `BezierGrid(mesh)` to create a grid of Bézier elements with extraction operators for analysis.
See [`is_analysis_suitable`](@ref) for which T-meshes give a well behaved basis.
"""
struct TSplineMesh{pdim,sdim,T}
	order::Int
	knots::NTuple{pdim,Vector{T}}
	cells::Vector{NTuple{pdim,NTuple{2,Int}}} #Boxes of the T-mesh in index space, ((i0,i1), (j0,j1), ...)
	#walls[d] contains the unit walls (edges in 2D, faces in 3D) orthogonal to direction d. A key `I` denotes the wall
	#at the index I[d] in direction d, covering the unit cell with lower corner I[e] in the other directions e.
	walls::NTuple{pdim,Set{NTuple{pdim,Int}}}
	anchors::Vector{NTuple{pdim,Int}} #Index coordinates times 2, such that cell centers are integers
	local_indices::Vector{NTuple{pdim,Vector{Int}}} #Local index vectors (order+2 indices per direction)
	local_knots::Vector{NTuple{pdim,Vector{T}}}
	control_points::Vector{Vec{sdim,T}}
	weights::Vector{T}
	elements::Vector{NTuple{pdim,NTuple{2,T}}} #Bézier elements in the parametric domain, ((ξa, ξb), (ηa, ηb), ...)
	IEN::Vector{Vector{Int}}
	domain::NTuple{pdim,NTuple{2,T}} #The parametric domain ((ξ0, ξ1), (η0, η1), ...)
end

function TSplineMesh(knots::NTuple{pdim,Vector{T}}, order::Int, cells::AbstractVector{<:NTuple{N,Int}};
		control_points::Union{Nothing,AbstractVector{<:Vec}} = nothing,
		weights::Union{Nothing,AbstractVector} = nothing,
		domain = nothing) where {pdim,T,N}

	pdim in (2, 3) || throw(ArgumentError("Only parametric dimension 2 and 3 are supported."))
	N == 2pdim || throw(ArgumentError("The cells must be given as tuples with $(2pdim) indices (i0, i1, j0, j1, ...)."))
	order >= 1 || throw(ArgumentError("The order must be at least 1."))
	boxes = [ntuple(d -> (c[2d-1], c[2d]), pdim) for c in cells]
	walls, vertices = _validate_tmesh(knots, order, boxes)
	full_domain = ntuple(d -> (knots[d][1], knots[d][end]), pdim)
	domain = domain === nothing ? full_domain : ntuple(d -> T.(domain[d]), pdim)

	all_anchors = _compute_anchors(order, boxes, vertices)
	all_indices = [_local_index_vectors(knots, order, walls, A) for A in all_anchors]
	all_local_knots = [ntuple(d -> knots[d][I[d]], pdim) for I in all_indices]
	#Basis functions with zero support (only possible outside the domain) are not considered
	nonzero = [all(K -> K[1] < K[end], Ks) for Ks in all_local_knots]
	all_elements, all_IEN = _compute_elements_and_ien(knots, walls, all_local_knots[nonzero])
	nonzero_ids = findall(nonzero)

	#Keep the elements inside the domain, and the basis functions which are non-zero on them
	keep_elements = Int[]
	for (e, bounds) in enumerate(all_elements)
		inside  = all(d -> domain[d][1] <= bounds[d][1] && bounds[d][2] <= domain[d][2], 1:pdim)
		outside = any(d -> bounds[d][2] <= domain[d][1] || bounds[d][1] >= domain[d][2], 1:pdim)
		inside || outside || throw(ArgumentError("The Bézier element $bounds crosses the boundary of the domain $domain. The bounds of the domain must be lines of the T-mesh."))
		inside && push!(keep_elements, e)
	end
	isempty(keep_elements) && throw(ArgumentError("The domain $domain does not contain any elements."))
	elements = all_elements[keep_elements]
	kept = sort!(unique(reduce(vcat, all_IEN[keep_elements])))
	renumber = Dict(A => i for (i, A) in enumerate(kept))
	IEN = [[renumber[A] for A in ien] for ien in all_IEN[keep_elements]]
	anchors = all_anchors[nonzero_ids[kept]]
	local_indices = all_indices[nonzero_ids[kept]]
	local_knots = all_local_knots[nonzero_ids[kept]]
	#Without a frame, all basis functions must be non-zero in the domain
	if length(anchors) != length(all_anchors) && domain == full_domain
		A = first(setdiff(all_anchors, anchors))
		throw(ArgumentError("The basis function anchored at index $(A./2) has a support with zero measure."))
	end

	if control_points === nothing
		control_points = _greville_points(order, local_knots)
	end
	weights = weights === nothing ? ones(T, length(anchors)) : collect(T, weights)
	length(control_points) == length(anchors) || throw(ArgumentError("The T-mesh has $(length(anchors)) anchors (basis functions), but $(length(control_points)) control points were given."))
	length(weights) == length(anchors) || throw(ArgumentError("The T-mesh has $(length(anchors)) anchors (basis functions), but $(length(weights)) weights were given."))

	sdim = length(first(control_points))
	return TSplineMesh{pdim,sdim,T}(order, knots, boxes, walls, anchors, local_indices, local_knots,
		collect(Vec{sdim,T}, control_points), weights, elements, IEN, domain)
end

"""
	TSplineMesh(mesh::NURBSMesh)

Represent a NURBS patch (2D or 3D, with equal orders in all directions) as a T-spline.
The T-mesh is the full tensor product mesh, so the basis is identical to the NURBS basis.
"""
function TSplineMesh(mesh::NURBSMesh{pdim,sdim,T}) where {pdim,sdim,T}
	allequal(mesh.orders) || throw(ArgumentError("TSplineMesh requires equal orders in all directions."))
	p = first(mesh.orders)
	s = cld(p,2)
	knots = ntuple(d -> mesh.knot_vectors[d][(s+1):(end-s)], pdim)
	n = length.(knots)
	cells = [Tuple(Iterators.flatten((I[d], I[d]+1) for d in 1:pdim)) for I in CartesianIndices(n .- 1)]

	#The anchors are sorted with the first index fastest, which is the same ordering as the NURBS control points.
	tmesh = TSplineMesh(knots, p, vec(cells))
	nbf = length.(mesh.knot_vectors) .- p .- 1
	anchor_index(A) = isodd(p) ? A .÷ 2 : (A .- 1) .÷ 2
	order_map = [LinearIndices(nbf)[anchor_index(A)...] for A in tmesh.anchors]
	@assert sort(order_map) == 1:prod(nbf)

	return TSplineMesh(knots, p, vec(cells); control_points = mesh.control_points[order_map], weights = mesh.weights[order_map])
end

function Base.show(io::IO, mesh::TSplineMesh{pdim}) where {pdim}
    print(io, "$(pdim)D TSplineMesh with order $(mesh.order), $(getncells(mesh)) Bézier elements and $(getnnodes(mesh)) control points")
end

Ferrite.getncells(mesh::TSplineMesh) = length(mesh.elements)
Ferrite.getnnodes(mesh::TSplineMesh) = length(mesh.control_points)

"""
	greville_points(mesh::TSplineMesh)

The Greville abscissae (mean of the interior knots of the local knot vectors) of each basis function,
in the same order as the control points. Using these as control points gives the identity map from the
parametric domain; for an affine geometry, map the Greville points with the affine map.
"""
greville_points(mesh::TSplineMesh) = _greville_points(mesh.order, mesh.local_knots)

function _greville_points(p::Int, local_knots::Vector{NTuple{pdim,Vector{T}}}) where {pdim,T}
	return [Vec{pdim,T}(ntuple(d -> sum(Ks[d][2:p+1])/p, pdim)) for Ks in local_knots]
end

# # #
# T-mesh construction
# # #

#Unit cells in the other directions than `d`, as ranges over the box
_other_ranges(box, d) = ntuple(e -> e == d ? (0:0) : (box[e][1]:box[e][2]-1), length(box))

function _validate_tmesh(knots::NTuple{pdim,Vector{T}}, p::Int, boxes::Vector{NTuple{pdim,NTuple{2,Int}}}) where {pdim,T}
	n = length.(knots)
	for d in 1:pdim
		issorted(knots[d]) || throw(ArgumentError("The knots in direction $d are not non-decreasing."))
		knots[d][1] < knots[d][end] || throw(ArgumentError("The knots in direction $d span an empty interval."))
	end

	#Check that the cells tile the index domain
	count = zeros(Int, (n .- 1)...)
	for box in boxes
		if !all(d -> 1 <= box[d][1] < box[d][2] <= n[d], 1:pdim)
			throw(ArgumentError("Cell $box is not a non-empty box inside the index domain $(Tuple((1, nd) for nd in n))."))
		end
		count[ntuple(d -> box[d][1]:box[d][2]-1, pdim)...] .+= 1
	end
	all(==(1), count) || throw(ArgumentError("The cells do not tile the index domain (overlaps or gaps)."))

	walls = ntuple(_ -> Set{NTuple{pdim,Int}}(), pdim)
	vertices = Set{NTuple{pdim,Int}}()
	for box in boxes
		for d in 1:pdim, i in box[d], I in CartesianIndices(_other_ranges(box, d))
			push!(walls[d], Base.setindex(Tuple(I), i, d))
		end
		for corner in Iterators.product(box...)
			push!(vertices, corner)
		end
	end

	#Two vertices on opposite sides of a face, with the same index coordinate, must be connected (only checked in 2D)
	if pdim == 2
		for ((i0, i1), (j0, j1)) in boxes
			bottom = Set(i for i in i0+1:i1-1 if (i,j0) in vertices)
			top    = Set(i for i in i0+1:i1-1 if (i,j1) in vertices)
			left   = Set(j for j in j0+1:j1-1 if (i0,j) in vertices)
			right  = Set(j for j in j0+1:j1-1 if (i1,j) in vertices)
			if !isempty(intersect(bottom, top)) || !isempty(intersect(left, right))
				throw(ArgumentError("The T-mesh is not admissible: face $((i0,i1,j0,j1)) has vertices on opposite sides with the same index coordinate, which must be connected by an edge."))
			end
		end
	end

	#The zero-length knot intervals at the boundary (the frame) must consist of full lines (planes)
	r = p÷2 + 1
	for d in 1:pdim, i in 2:r, (ii, iboundary) in ((i, 1), (n[d]-i+1, n[d]))
		knots[d][ii] == knots[d][iboundary] || continue
		full = all(I -> Base.setindex(Tuple(I), ii, d) in walls[d], CartesianIndices(ntuple(e -> e == d ? (0:0) : (1:n[e]-1), pdim)))
		full || throw(ArgumentError("The index line (plane) $ii in direction $d, at the boundary with zero knot interval, must span the full domain."))
	end

	return walls, vertices
end

function _compute_anchors(p::Int, boxes::Vector{NTuple{pdim,NTuple{2,Int}}}, vertices) where {pdim}
	anchors = if isodd(p)
		[2 .* v for v in vertices]
	else
		[ntuple(d -> box[d][1] + box[d][2], pdim) for box in boxes]
	end
	return sort!(collect(anchors), by = reverse)
end

"""
	_lines_crossing(walls, n, d, X)

Ray-cast through the T-mesh in direction `d`, along the line through the point with (doubled) index coordinates `X`.
Return the sorted indices of the lines (planes) orthogonal to `d` which contain the line point.
Points on the boundary of a wall are included, such that a line ending in a T-junction on the ray is found.
"""
function _lines_crossing(walls::NTuple{pdim,Set{NTuple{pdim,Int}}}, n::Int, d::Int, X::NTuple{pdim,Int}) where {pdim}
	#The unit walls that can contain the point: for an odd coordinate the point is inside a unit cell,
	#for an even coordinate it is on a line, so the walls on both sides contain it.
	candidates = ntuple(e -> e == d ? (0,) : (isodd(X[e]) ? ((X[e]-1)÷2,) : (X[e]÷2 - 1, X[e]÷2)), pdim)
	lines = Int[]
	for i in 1:n
		if any(I -> Base.setindex(I, i, d) in walls[d], Iterators.product(candidates...))
			push!(lines, i)
		end
	end
	return lines
end

function _local_index_vector(lines::Vector{Int}, p::Int, X::Int)
	if isodd(p)
		k = (p+1)÷2
		i = X÷2
		@assert i in lines
		left  = filter(<(i), lines)
		right = filter(>(i), lines)
		left  = left[max(1,end-k+1):end]
		right = right[1:min(k,end)]
		first_ = isempty(left) ? i : first(left)
		last_  = isempty(right) ? i : last(right)
		return [fill(first_, k - length(left)); left; i; right; fill(last_, k - length(right))]
	else
		k = p÷2 + 1
		left  = filter(x -> 2x < X, lines)
		right = filter(x -> 2x > X, lines)
		@assert !isempty(left) && !isempty(right)
		left  = left[max(1,end-k+1):end]
		right = right[1:min(k,end)]
		return [fill(first(left), k - length(left)); left; right; fill(last(right), k - length(right))]
	end
end

function _local_index_vectors(knots::NTuple{pdim,Vector{T}}, p::Int, walls, A::NTuple{pdim,Int}) where {pdim,T}
	return ntuple(d -> _local_index_vector(_lines_crossing(walls, length(knots[d]), d, A), p, A[d]), pdim)
end

#Find the Bézier elements, i.e. the regions bounded by lines (planes) of reduced continuity of all basis functions
#(the elemental T-mesh), and the basis functions that are non-zero on each element.
function _compute_elements_and_ien(knots::NTuple{pdim,Vector{T}}, walls, local_knots) where {pdim,T}
	vals = ntuple(d -> unique(knots[d]), pdim)
	nf = ntuple(d -> length(vals[d]) - 1, pdim) #Number of cells of the fine tensor product grid
	idx(d, v) = searchsortedfirst(vals[d], v)

	#cover[d][I] is true if the interface orthogonal to d at the value vals[d][I[d]], between the fine cells
	#I - e_d and I, is covered by a knot line (plane) of some basis function or a wall of the T-mesh.
	cover = ntuple(d -> falses(ntuple(e -> e == d ? nf[e]+1 : nf[e], pdim)), pdim)
	for Ks in local_knots
		ranges = ntuple(e -> idx(e, Ks[e][1]):idx(e, Ks[e][end])-1, pdim)
		for d in 1:pdim, v in Ks[d]
			cover[d][Base.setindex(ranges, idx(d, v):idx(d, v), d)...] .= true
		end
	end
	for d in 1:pdim, I in walls[d]
		ranges = ntuple(e -> e == d ? (idx(d, knots[d][I[d]]):idx(d, knots[d][I[d]])) :
		                              (idx(e, knots[e][I[e]]):idx(e, knots[e][I[e]+1])-1), pdim)
		cover[d][ranges...] .= true
	end

	#Merge fine cells not separated by a knot line (union-find)
	lin = LinearIndices(nf)
	parent = collect(1:prod(nf))
	findroot(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
	unite(x, y) = (rx = findroot(x); ry = findroot(y); rx != ry && (parent[ry] = rx))
	for I in CartesianIndices(nf), d in 1:pdim
		if I[d] < nf[d]
			J = Base.setindex(Tuple(I), I[d]+1, d)
			cover[d][J...] || unite(lin[I], lin[J...])
		end
	end

	bbox = Dict{Int,Tuple{Vector{Int},Vector{Int},Int}}() #root => (mins, maxs, ncells)
	for I in CartesianIndices(nf)
		r = findroot(lin[I])
		mins, maxs, ncells = get(bbox, r, (collect(Tuple(I)), collect(Tuple(I)), 0))
		mins .= min.(mins, Tuple(I))
		maxs .= max.(maxs, Tuple(I))
		bbox[r] = (mins, maxs, ncells + 1)
	end
	roots = sort!(collect(keys(bbox)), by = r -> reverse(bbox[r][1]))
	elements = NTuple{pdim,NTuple{2,T}}[]
	element_of_cell = zeros(Int, nf)
	for (e, r) in enumerate(roots)
		mins, maxs, ncells = bbox[r]
		ncells == prod(maxs .- mins .+ 1) || error("Internal error: a Bézier element of the T-spline is not a box.")
		push!(elements, ntuple(d -> (vals[d][mins[d]], vals[d][maxs[d]+1]), pdim))
		element_of_cell[ntuple(d -> mins[d]:maxs[d], pdim)...] .= e
	end

	IEN = [Int[] for _ in elements]
	for (A, Ks) in enumerate(local_knots)
		for e in unique(element_of_cell[ntuple(d -> idx(d, Ks[d][1]):idx(d, Ks[d][end])-1, pdim)...])
			push!(IEN[e], A)
		end
	end
	sort!.(IEN)
	return elements, IEN
end

# # #
# Analysis suitability
# # #

"""
	tjunction_extensions(mesh::TSplineMesh)

Return the extensions of all T-junctions of a 2D T-mesh, as two vectors of index space
segments: horizontal extensions `(i1, i2, j)` and vertical extensions `(i, j1, j2)`.

Following Beirão da Veiga, Buffa, Sangalli and Vázquez (2013), the extension of a T-junction spans `p` bays:
the *face extension* extends `⌊(p+1)/2⌋` bays into the face where the edge is missing, and the *edge extension*
extends `⌈(p-1)/2⌉` bays in the opposite direction. That is, `(p+1)/2` and `(p-1)/2` bays for odd degrees,
and `p/2` bays in both directions for even degrees.
"""
function tjunction_extensions(mesh::TSplineMesh{pdim}) where {pdim}
	pdim == 2 || throw(ArgumentError("T-junction extensions are only defined for 2D T-meshes. Use is_dual_compatible in 3D."))
	p = mesh.order
	vedges, hedges = mesh.walls
	n1, n2 = length.(mesh.knots)
	kf = cld(p, 2) # = ⌊(p+1)/2⌋
	ke = p ÷ 2     # = ⌈(p-1)/2⌉
	hext = NTuple{3,Int}[]
	vext = NTuple{3,Int}[]
	#The T-junctions are vertices of the T-mesh (for even degrees, the anchors are faces, not vertices)
	vertices = Set{NTuple{2,Int}}()
	for box in mesh.cells, corner in Iterators.product(box...)
		push!(vertices, corner)
	end
	for (i, j) in vertices
		(1 < i < n1 && 1 < j < n2) || continue
		left  = (i-1,j) in hedges
		right = (i,j)   in hedges
		down  = (i,j-1) in vedges
		up    = (i,j)   in vedges
		(left + right + down + up) == 3 || continue
		if !left || !right
			lines = _lines_crossing(mesh.walls, n1, 1, (2i, 2j))
			pos = findfirst(==(i), lines)
			nleft, nright = !left ? (kf, ke) : (ke, kf)
			push!(hext, (lines[max(pos-nleft, 1)], lines[min(pos+nright, end)], j))
		else
			lines = _lines_crossing(mesh.walls, n2, 2, (2i, 2j))
			pos = findfirst(==(j), lines)
			ndown, nup = !down ? (kf, ke) : (ke, kf)
			push!(vext, (i, lines[max(pos-ndown, 1)], lines[min(pos+nup, end)]))
		end
	end
	return hext, vext
end

"""
	is_analysis_suitable(mesh::TSplineMesh)

Check if the T-mesh is analysis-suitable.

In 2D, a T-mesh is analysis-suitable if no horizontal T-junction extension intersects a vertical T-junction extension
(see [`tjunction_extensions`](@ref)). Analysis-suitable T-meshes are exactly the dual-compatible T-meshes
(Beirão da Veiga, Buffa, Sangalli and Vázquez, 2013). In 3D, where T-junction extensions are not established,
dual-compatibility is used as the definition (see [`is_dual_compatible`](@ref)).

For analysis-suitable T-splines of any degree, the basis functions are linearly independent, form a partition of unity,
reproduce polynomials of degree `p` and have optimal approximation properties. They are also locally linearly
independent: exactly `(p+1)^dim` basis functions are non-zero on each Bézier element.
"""
function is_analysis_suitable(mesh::TSplineMesh{pdim}) where {pdim}
	pdim == 3 && return is_dual_compatible(mesh)
	hext, vext = tjunction_extensions(mesh)
	for (i1, i4, j) in hext, (i, j1, j4) in vext
		if i1 <= i <= i4 && j1 <= j <= j4
			return false
		end
	end
	return true
end

"""
	is_dual_compatible(mesh::TSplineMesh)

Check if the T-mesh is dual-compatible, i.e. if the tensor products of the univariate B-spline dual functionals
of the local knot vectors form a dual basis for the T-spline basis functions. This implies linear independence,
local linear independence, partition of unity (if the constants are in the space), and optimal approximation properties.

The check is done in index space: for every pair of basis functions with overlapping index ranges, there must be a
direction in which the local index vectors *overlap* (are consecutive indices of a common global index vector, i.e. every
index of one vector within the range of the other is also an index of the other) and *differ*.
In 2D, this is equivalent to the dual-compatibility of Beirão da Veiga et al. (2013), and to analysis-suitability.
"""
function is_dual_compatible(mesh::TSplineMesh{pdim}) where {pdim}
	n = length(mesh.local_indices)
	sets = [ntuple(d -> Set(I[d]), pdim) for I in mesh.local_indices]
	lo(A, d) = first(mesh.local_indices[A][d])
	hi(A, d) = last(mesh.local_indices[A][d])
	overlap(A, B, d) = all(k -> !(lo(B, d) <= k <= hi(B, d)) || k in sets[B][d], sets[A][d]) &&
	                   all(k -> !(lo(A, d) <= k <= hi(A, d)) || k in sets[A][d], sets[B][d])
	#Sweep over the anchors sorted by their first index in direction 1
	order = sortperm(1:n, by = A -> lo(A, 1))
	for (ia, A) in enumerate(order)
		for B in view(order, ia+1:n)
			lo(B, 1) > hi(A, 1) && break
			all(d -> lo(B, d) <= hi(A, d) && lo(A, d) <= hi(B, d), 1:pdim) || continue
			any(d -> mesh.local_indices[A][d] != mesh.local_indices[B][d] && overlap(A, B, d), 1:pdim) || return false
		end
	end
	return true
end

"""
	is_linearly_independent(mesh::TSplineMesh; rtol)

Numerically check if the T-spline basis functions are linearly independent, by computing the rank of the matrix of
Bézier coefficients of all basis functions. A dense matrix is used, so this is only intended for meshes of moderate size.

Analysis-suitable T-meshes are always linearly independent (for any knot values), but the opposite does not hold.
"""
function is_linearly_independent(mesh::TSplineMesh{pdim,sdim,T}; rtol = sqrt(eps(T))) where {pdim,sdim,T}
	beos = compute_bezier_extraction_operators(mesh)
	nb = (mesh.order+1)^pdim
	G = zeros(T, getnnodes(mesh), getncells(mesh)*nb)
	for (e, beo) in enumerate(beos)
		for (a, A) in enumerate(mesh.IEN[e])
			G[A, (e-1)*nb .+ (1:nb)] .= beo.C[a]
		end
	end
	return LinearAlgebra.rank(G; rtol) == getnnodes(mesh)
end
