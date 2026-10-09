"""
	TSplineMesh(knots, order, faces; control_points, weights)

A (rational) T-spline with parametric dimension 2 and polynomial degree `order` in both directions.

The T-mesh is given in *index space*:
- `knots::NTuple{2,Vector{T}}` are the knot values `ξᵢ` and `ηⱼ` of the index lines `i = 1:n1` and `j = 1:n2`.
  Consecutive knots may be equal (zero knot intervals). An open (interpolatory) boundary is obtained by repeating
  the boundary value `order÷2 + 1` times, e.g. `[0,0,1,2,3,3]` for cubic splines.
- `faces::Vector{NTuple{4,Int}}` are the rectangles `(i0, i1, j0, j1)` of the T-mesh. They must tile the index domain `[1,n1]×[1,n2]`.

One basis function is associated with each *anchor*: the vertices of the T-mesh for odd degrees, and the
centers of the faces for even degrees. The anchors are sorted row by row (bottom to top, left to right), and
`control_points` and `weights` must be given in this order. If `control_points` are not given, the Greville points
of the basis functions are used, which gives the identity map from the parametric domain (see [`greville_points`](@ref)).

Use `BezierGrid(mesh)` to create a grid of Bézier elements with extraction operators for analysis.
See [`is_analysis_suitable`](@ref) for which T-meshes give a well behaved basis.
"""
struct TSplineMesh{sdim,T}
	order::Int
	knots::NTuple{2,Vector{T}}
	faces::Vector{NTuple{4,Int}}
	hedges::Set{NTuple{2,Int}} #Unit horizontal edges: (i,j) denotes the edge (i,j)-(i+1,j)
	vedges::Set{NTuple{2,Int}} #Unit vertical edges:   (i,j) denotes the edge (i,j)-(i,j+1)
	anchors::Vector{NTuple{2,Int}} #Index coordinates times 2, such that face centers are integers
	local_knots::Vector{NTuple{2,Vector{T}}}
	control_points::Vector{Vec{sdim,T}}
	weights::Vector{T}
	elements::Vector{NTuple{4,T}} #Bézier elements in the parametric domain, (ξa, ξb, ηa, ηb)
	IEN::Vector{Vector{Int}}
end

function TSplineMesh(knots::NTuple{2,Vector{T}}, order::Int, faces::AbstractVector{NTuple{4,Int}};
		control_points::Union{Nothing,AbstractVector{<:Vec}} = nothing,
		weights::Union{Nothing,AbstractVector} = nothing) where {T}

	order >= 1 || throw(ArgumentError("The order must be at least 1."))
	faces = collect(faces)
	hedges, vedges, vertices = _validate_tmesh(knots, order, faces)

	anchors = _compute_anchors(order, faces, vertices)
	local_knots = [_local_knot_vectors(knots, order, hedges, vedges, A) for A in anchors]
	for (A, (Ξ, Π)) in zip(anchors, local_knots)
		if Ξ[1] == Ξ[end] || Π[1] == Π[end]
			throw(ArgumentError("The basis function anchored at index $(A./2) has a support with zero area."))
		end
	end
	elements, IEN = _compute_elements_and_ien(knots, hedges, vedges, local_knots)

	if control_points === nothing
		control_points = _greville_points(order, local_knots)
	end
	weights = weights === nothing ? ones(T, length(anchors)) : collect(T, weights)
	length(control_points) == length(anchors) || throw(ArgumentError("The T-mesh has $(length(anchors)) anchors (basis functions), but $(length(control_points)) control points were given."))
	length(weights) == length(anchors) || throw(ArgumentError("The T-mesh has $(length(anchors)) anchors (basis functions), but $(length(weights)) weights were given."))

	sdim = length(first(control_points))
	return TSplineMesh{sdim,T}(order, knots, faces, hedges, vedges, anchors, local_knots,
		collect(Vec{sdim,T}, control_points), weights, elements, IEN)
end

"""
	TSplineMesh(mesh::NURBSMesh)

Represent a two-dimensional NURBS patch (with equal orders in both directions) as a T-spline.
The T-mesh is the full tensor product mesh, so the basis is identical to the NURBS basis.
"""
function TSplineMesh(mesh::NURBSMesh{2,sdim,T}) where {sdim,T}
	allequal(mesh.orders) || throw(ArgumentError("TSplineMesh requires equal orders in both directions."))
	p = first(mesh.orders)
	s = cld(p,2)
	knots = ntuple(d -> mesh.knot_vectors[d][(s+1):(end-s)], 2)
	n1, n2 = length.(knots)
	faces = [(i, i+1, j, j+1) for j in 1:(n2-1) for i in 1:(n1-1)]

	#The anchors are sorted row by row, which is the same ordering as the NURBS control points.
	tmesh = TSplineMesh(knots, p, faces)
	nbf = length.(mesh.knot_vectors) .- p .- 1
	cp_index(A) = isodd(p) ? (A[1]÷2) + (A[2]÷2 - 1)*nbf[1] : (A[1]-1)÷2 + ((A[2]-1)÷2 - 1)*nbf[1]
	order_map = [cp_index(A) for A in tmesh.anchors]
	@assert sort(order_map) == 1:prod(nbf)

	return TSplineMesh(knots, p, faces; control_points = mesh.control_points[order_map], weights = mesh.weights[order_map])
end

function Base.show(io::IO, mesh::TSplineMesh)
    print(io, "TSplineMesh with order $(mesh.order), $(getncells(mesh)) Bézier elements and $(getnnodes(mesh)) control points")
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

function _greville_points(p::Int, local_knots::Vector{NTuple{2,Vector{T}}}) where T
	return [Vec{2,T}((sum(Ξ[2:p+1])/p, sum(Π[2:p+1])/p)) for (Ξ, Π) in local_knots]
end

# # #
# T-mesh construction
# # #

function _validate_tmesh(knots::NTuple{2,Vector{T}}, p::Int, faces::Vector{NTuple{4,Int}}) where T
	n1, n2 = length.(knots)
	for d in 1:2
		issorted(knots[d]) || throw(ArgumentError("The knots in direction $d are not non-decreasing."))
		knots[d][1] < knots[d][end] || throw(ArgumentError("The knots in direction $d span an empty interval."))
	end

	#Check that the faces tile the index domain
	count = zeros(Int, n1-1, n2-1)
	for (i0, i1, j0, j1) in faces
		if !(1 <= i0 < i1 <= n1 && 1 <= j0 < j1 <= n2)
			throw(ArgumentError("Face $((i0,i1,j0,j1)) is not a non-empty rectangle inside the index domain [1,$n1]×[1,$n2]."))
		end
		count[i0:i1-1, j0:j1-1] .+= 1
	end
	all(==(1), count) || throw(ArgumentError("The faces do not tile the index domain (overlaps or gaps)."))

	hedges = Set{NTuple{2,Int}}()
	vedges = Set{NTuple{2,Int}}()
	vertices = Set{NTuple{2,Int}}()
	for (i0, i1, j0, j1) in faces
		for i in i0:i1-1
			push!(hedges, (i,j0), (i,j1))
		end
		for j in j0:j1-1
			push!(vedges, (i0,j), (i1,j))
		end
		push!(vertices, (i0,j0), (i1,j0), (i1,j1), (i0,j1))
	end

	#Two vertices on opposite sides of a face, with the same index coordinate, must be connected
	for (i0, i1, j0, j1) in faces
		bottom = Set(i for i in i0+1:i1-1 if (i,j0) in vertices)
		top    = Set(i for i in i0+1:i1-1 if (i,j1) in vertices)
		left   = Set(j for j in j0+1:j1-1 if (i0,j) in vertices)
		right  = Set(j for j in j0+1:j1-1 if (i1,j) in vertices)
		if !isempty(intersect(bottom, top)) || !isempty(intersect(left, right))
			throw(ArgumentError("The T-mesh is not admissible: face $((i0,i1,j0,j1)) has vertices on opposite sides with the same index coordinate, which must be connected by an edge."))
		end
	end

	#The zero-length knot intervals at the boundary (the frame) must consist of full lines
	r = p÷2 + 1
	for i in 2:r
		for ii in (i, n1-i+1)
			knots[1][ii] == knots[1][ii == i ? 1 : n1] || continue
			all(j -> (ii,j) in vedges, 1:n2-1) || throw(ArgumentError("The vertical index line i=$ii, at the boundary with zero knot interval, must span the full domain."))
		end
	end
	for j in 2:r
		for jj in (j, n2-j+1)
			knots[2][jj] == knots[2][jj == j ? 1 : n2] || continue
			all(i -> (i,jj) in hedges, 1:n1-1) || throw(ArgumentError("The horizontal index line j=$jj, at the boundary with zero knot interval, must span the full domain."))
		end
	end

	return hedges, vedges, vertices
end

function _compute_anchors(p::Int, faces::Vector{NTuple{4,Int}}, vertices::Set{NTuple{2,Int}})
	anchors = if isodd(p)
		[(2i, 2j) for (i,j) in vertices]
	else
		[(i0+i1, j0+j1) for (i0, i1, j0, j1) in faces]
	end
	return sort!(anchors, by = A -> (A[2], A[1]))
end

"""
	_lines_crossing(hedges, vedges, n, dir, X)

Ray-cast through the T-mesh. Return the sorted indices of the index lines orthogonal to the direction `dir`
which contain the point with (doubled) index coordinate `X` in the other direction.
Vertices (end points of the lines) are included, such that a line ending in a T-junction on the ray is found.
"""
function _lines_crossing(hedges::Set{NTuple{2,Int}}, vedges::Set{NTuple{2,Int}}, n::Int, dir::Int, X::Int)
	lines = Int[]
	edges = dir == 1 ? vedges : hedges
	key(i, j) = dir == 1 ? (i,j) : (j,i)
	for i in 1:n
		hit = if iseven(X)
			(key(i, X÷2) in edges) || (key(i, X÷2 - 1) in edges)
		else
			key(i, (X-1)÷2) in edges
		end
		hit && push!(lines, i)
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

function _local_knot_vectors(knots::NTuple{2,Vector{T}}, p::Int, hedges, vedges, A::NTuple{2,Int}) where T
	n1, n2 = length.(knots)
	ilines = _lines_crossing(hedges, vedges, n1, 1, A[2])
	jlines = _lines_crossing(hedges, vedges, n2, 2, A[1])
	Ξ = knots[1][_local_index_vector(ilines, p, A[1])]
	Π = knots[2][_local_index_vector(jlines, p, A[2])]
	return (Ξ, Π)
end

#Find the Bézier elements, i.e. the regions bounded by lines of reduced continuity of all basis functions
#(the elemental T-mesh), and the basis functions that are non-zero on each element.
function _compute_elements_and_ien(knots::NTuple{2,Vector{T}}, hedges, vedges, local_knots) where T
	uvals = unique(knots[1])
	vvals = unique(knots[2])
	nu = length(uvals) - 1
	nv = length(vvals) - 1
	uidx(u) = searchsortedfirst(uvals, u)
	vidx(v) = searchsortedfirst(vvals, v)

	#Interfaces between the cells of the fine tensor product grid which are covered by a knot line
	coverV = falses(nu+1, nv) # line at uvals[a], between vvals[b] and vvals[b+1]
	coverH = falses(nu, nv+1) # line at vvals[b], between uvals[a] and uvals[a+1]
	for (Ξ, Π) in local_knots
		b0, b1 = vidx(Π[1]), vidx(Π[end])
		a0, a1 = uidx(Ξ[1]), uidx(Ξ[end])
		for ξ in Ξ
			coverV[uidx(ξ), b0:b1-1] .= true
		end
		for η in Π
			coverH[a0:a1-1, vidx(η)] .= true
		end
	end
	for (i,j) in vedges
		coverV[uidx(knots[1][i]), vidx(knots[2][j]):vidx(knots[2][j+1])-1] .= true
	end
	for (i,j) in hedges
		coverH[uidx(knots[1][i]):uidx(knots[1][i+1])-1, vidx(knots[2][j])] .= true
	end

	#Merge fine cells not separated by a knot line (union-find)
	lin = LinearIndices((nu, nv))
	parent = collect(1:nu*nv)
	findroot(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
	unite(x, y) = (rx = findroot(x); ry = findroot(y); rx != ry && (parent[ry] = rx))
	for b in 1:nv, a in 1:nu
		a < nu && !coverV[a+1, b] && unite(lin[a,b], lin[a+1,b])
		b < nv && !coverH[a, b+1] && unite(lin[a,b], lin[a,b+1])
	end

	bbox = Dict{Int,Vector{Int}}() #root => [amin, amax, bmin, bmax, ncells]
	for b in 1:nv, a in 1:nu
		r = findroot(lin[a,b])
		bb = get!(bbox, r, [a, a, b, b, 0])
		bb[1] = min(bb[1], a); bb[2] = max(bb[2], a)
		bb[3] = min(bb[3], b); bb[4] = max(bb[4], b)
		bb[5] += 1
	end
	roots = sort!(collect(keys(bbox)), by = r -> (bbox[r][3], bbox[r][1]))
	elements = NTuple{4,T}[]
	element_of_cell = zeros(Int, nu, nv)
	for (e, r) in enumerate(roots)
		amin, amax, bmin, bmax, ncells = bbox[r]
		ncells == (amax-amin+1)*(bmax-bmin+1) || error("Internal error: a Bézier element of the T-spline is not rectangular.")
		push!(elements, (uvals[amin], uvals[amax+1], vvals[bmin], vvals[bmax+1]))
		element_of_cell[amin:amax, bmin:bmax] .= e
	end

	IEN = [Int[] for _ in elements]
	for (A, (Ξ, Π)) in enumerate(local_knots)
		for e in unique(element_of_cell[uidx(Ξ[1]):uidx(Ξ[end])-1, vidx(Π[1]):vidx(Π[end])-1])
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

Return the extensions of all T-junctions of the T-mesh, as two vectors of index space
segments: horizontal extensions `(i1, i2, j)` and vertical extensions `(i, j1, j2)`.

Following Beirão da Veiga, Buffa, Sangalli and Vázquez (2013), the extension of a T-junction spans `p` bays:
the *face extension* extends `⌊(p+1)/2⌋` bays into the face where the edge is missing, and the *edge extension* 
extends `⌈(p-1)/2⌉` bays in the opposite direction. That is, `(p+1)/2` and `(p-1)/2` bays for odd degrees, 
and `p/2` bays in both directions for even degrees.
"""
function tjunction_extensions(mesh::TSplineMesh)
	p = mesh.order
	(; hedges, vedges) = mesh
	n1, n2 = length.(mesh.knots)
	kf = cld(p, 2) # = ⌊(p+1)/2⌋
	ke = p ÷ 2     # = ⌈(p-1)/2⌉
	hext = NTuple{3,Int}[]
	vext = NTuple{3,Int}[]
	#The T-junctions are vertices of the T-mesh (for even degrees, the anchors are faces, not vertices)
	vertices = Set{NTuple{2,Int}}()
	for (i0, i1, j0, j1) in mesh.faces
		push!(vertices, (i0,j0), (i1,j0), (i1,j1), (i0,j1))
	end
	for (i, j) in vertices
		(1 < i < n1 && 1 < j < n2) || continue
		left  = (i-1,j) in hedges
		right = (i,j)   in hedges
		down  = (i,j-1) in vedges
		up    = (i,j)   in vedges
		(left + right + down + up) == 3 || continue
		if !left || !right
			lines = _lines_crossing(hedges, vedges, n1, 1, 2j)
			pos = findfirst(==(i), lines)
			nleft, nright = !left ? (kf, ke) : (ke, kf)
			push!(hext, (lines[max(pos-nleft, 1)], lines[min(pos+nright, end)], j))
		else
			lines = _lines_crossing(hedges, vedges, n2, 2, 2i)
			pos = findfirst(==(j), lines)
			ndown, nup = !down ? (kf, ke) : (ke, kf)
			push!(vext, (i, lines[max(pos-ndown, 1)], lines[min(pos+nup, end)]))
		end
	end
	return hext, vext
end

"""
	is_analysis_suitable(mesh::TSplineMesh)

Check if the T-mesh is analysis-suitable, i.e. if no horizontal T-junction extension intersects
a vertical T-junction extension (see [`tjunction_extensions`](@ref)).
For analysis-suitable T-splines of any degree, the basis functions are linearly independent, form a partition of unity,
reproduce polynomials of degree `p` and have optimal approximation properties. Analysis-suitable T-meshes are exactly the
dual-compatible T-meshes (Beirão da Veiga, Buffa, Sangalli and Vázquez, 2013), so the basis functions are also locally 
linearly independent: exactly `(p+1)^2` basis functions are non-zero on each Bézier element.
"""
function is_analysis_suitable(mesh::TSplineMesh)
	hext, vext = tjunction_extensions(mesh)
	for (i1, i4, j) in hext, (i, j1, j4) in vext
		if i1 <= i <= i4 && j1 <= j <= j4
			return false
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
function is_linearly_independent(mesh::TSplineMesh{sdim,T}; rtol = sqrt(eps(T))) where {sdim,T}
	beos = compute_bezier_extraction_operators(mesh)
	nb = (mesh.order+1)^2
	G = zeros(T, getnnodes(mesh), getncells(mesh)*nb)
	for (e, beo) in enumerate(beos)
		for (a, A) in enumerate(mesh.IEN[e])
			G[A, (e-1)*nb .+ (1:nb)] .= beo.C[a]
		end
	end
	return LinearAlgebra.rank(G; rtol) == getnnodes(mesh)
end
