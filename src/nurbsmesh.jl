"""
    NURBSMesh(knot_vectors::NTuple{pdim, Vector{T}}, orders::NTuple{pdim, Int}, control_points::Vector{Vec{sdim, T}}, [weights::Vector{T}])

A single tensor-product B-spline or NURBS patch with parametric dimension `pdim` and
spatial dimension `sdim`.

**Arguments:**
* `knot_vectors`: one (non-decreasing) knot vector per parametric direction.
* `orders`: the polynomial degree in each parametric direction.
* `control_points`: the control points. They are ordered with the first parametric
  direction running fastest, i.e. control point ``(i, j)`` of a 2D patch with ``n_1``
  basis functions in the first direction has index ``i + (j - 1) n_1``.
* `weights`: the NURBS weight of each control point. Defaults to all ones, which gives a
  B-spline patch.

The number of control points must be ``\\prod_d (\\mathrm{length}(\\Xi_d) - p_d - 1)``.
The cells (elements) of the patch are the non-empty knot spans; knot spans of zero length
(from repeated knots) are skipped.

A `NURBSMesh` describes the geometry and the spline space. For analysis, it is converted to
a [`BezierGrid`](@ref). Simple patches can be created with [`generate_nurbs_patch`](@ref),
and refined with [`knotinsertion`](@ref), [`orderelevation`](@ref) and
[`smoothnesselevation`](@ref).

# Examples
A quarter of a circle with radius 1 as a quadratic NURBS curve:
```julia
knot_vectors = ([0.0, 0.0, 0.0, 1.0, 1.0, 1.0],)
control_points = [Vec((1.0, 0.0)), Vec((1.0, 1.0)), Vec((0.0, 1.0))]
weights = [1.0, 1 / sqrt(2), 1.0]
mesh = NURBSMesh(knot_vectors, (2,), control_points, weights)
```

See the [NURBSMesh topic guide](@ref nurbsmesh_topicguide) for more details.
"""
struct NURBSMesh{pdim,sdim,T}
	knot_vectors::NTuple{pdim,Vector{T}}
	orders::NTuple{pdim,Int}
	control_points::Vector{Vec{sdim,T}}
	weights::Vector{T}
	IEN::Matrix{Int}
	INN::Matrix{Int}

	function NURBSMesh(knot_vectors::NTuple{pdim,Vector{T}}, orders::NTuple{pdim,Int},
									control_points::Vector{Vec{sdim,T}}, 
									weights::AbstractVector{T}=ones(T, length(control_points))) where {pdim,sdim,T}

		pdim==3 && sdim==2 ? error("A 3d geometry can not exist in 2d") : nothing

		nbasefuncs = length.(knot_vectors) .- orders .- 1
		nel, nnp, nen, INN, IEN = generate_nurbs_meshdata(orders, nbasefuncs)
		
		@assert(prod(nbasefuncs)==maximum(IEN))

		#Remove elements which are zero length
		to_remove = Int[]
		for e in 1:nel
			nurbs_coords = [INN[IEN[end,e],d] for d in 1:pdim]
			for d in 1:pdim
				if knot_vectors[d][nurbs_coords[d]] == knot_vectors[d][nurbs_coords[d]+1]
					push!(to_remove, e)
					break
				end
			end
		end
		to_keep = setdiff(collect(1:nel), to_remove)

		IEN = IEN[:, to_keep] #IEN = IEN[end:-1:1, to_keep]
		new{pdim,sdim,T}(knot_vectors,orders,control_points,weights,IEN,INN)
	end

end

function Base.show(io::IO, ::MIME"text/plain", grid::NURBSMesh)
    print(io, "NURBSMesh and with order $(grid.orders) and $(getncells(grid)) cells and $(getnnodes(grid)) nodes (control points) ")
end

Ferrite.getncells(mesh::NURBSMesh) = size(mesh.IEN, 2)
Ferrite.getnnodes(mesh::NURBSMesh) = maximum(mesh.IEN) 

function Ferrite.getcoordinates(mesh::NURBSMesh, ie::Int)
	return mesh.control_points[mesh.IEN[:,ie]]
end

function _copy_nurbsmesh_data(mesh::NURBSMesh)
	return (map(copy, mesh.knot_vectors), mesh.orders, copy(mesh.control_points), copy(mesh.weights))
end

"""
    knotinsertion(mesh::NURBSMesh, ξ; dir::Int)

h-refinement: return a new [`NURBSMesh`](@ref) with the knot value `ξ` inserted once in
parametric direction `dir`. The geometry and the parametrization are unchanged, but the
spline space is enriched. Inserting a new value splits a knot span into two cells, with
continuity ``C^{p-1}`` across the new knot. Inserting an existing value instead increases its
multiplicity ``m`` and lowers the continuity there to ``C^{p-m}``.

The input mesh is not modified, see [`knotinsertion!`](@ref) for the in-place version.

# Examples
```julia
mesh = generate_nurbs_patch(:rectangle, (1, 1), 2; size = (1.0, 1.0))
mesh = knotinsertion(mesh, 0.0; dir = 1)  # two cells in the first direction
```
"""
function knotinsertion(mesh::NURBSMesh, ξ::T; dir::Int) where T
	knot_vectors, orders, control_points, weights = _copy_nurbsmesh_data(mesh)
	knotinsertion!(knot_vectors, orders, control_points, weights, ξ; dir)
	return NURBSMesh(knot_vectors, orders, control_points, weights)
end

"""
    orderelevation(mesh::NURBSMesh; dir::Int)

p-refinement: return a new [`NURBSMesh`](@ref) with the polynomial degree raised by one in
parametric direction `dir`. The geometry and parametrization are unchanged. The
multiplicity of every knot is increased by one, so the continuity across existing knots is
preserved.

Note that [`BezierGrid`](@ref) currently requires the same degree in all directions, so
elevate the degree in every direction before creating the grid.

The input mesh is not modified, see [`orderelevation!`](@ref) for the in-place version.
"""
function orderelevation(mesh::NURBSMesh; dir::Int)
	knot_vectors, orders, control_points, weights = _copy_nurbsmesh_data(mesh)
	orders = orderelevation!(knot_vectors, orders, control_points, weights; dir)
	return NURBSMesh(knot_vectors, orders, control_points, weights)
end

"""
    smoothnesselevation(mesh::NURBSMesh, new_knots::AbstractVector; dir::Int)

k-refinement: return a new [`NURBSMesh`](@ref) where the degree is first raised by one in
parametric direction `dir`, and then each value in `new_knots` is inserted once. Since the
knots are inserted *after* the degree elevation, the continuity across the new knots is
``C^{p}``, where ``p`` is the original degree. This gives a smoother spline space than doing
the two operations in the opposite order.

The input mesh is not modified, see [`smoothnesselevation!`](@ref) for the in-place version.
"""
function smoothnesselevation(mesh::NURBSMesh, new_knots::AbstractVector{T}; dir::Int) where T
	knot_vectors, orders, control_points, weights = _copy_nurbsmesh_data(mesh)
	orders = smoothnesselevation!(knot_vectors, orders, control_points, weights, new_knots; dir)
	return NURBSMesh(knot_vectors, orders, control_points, weights)
end

"""
    eval_parametric_coordinate(mesh::NURBSMesh{pdim}, ξ::Vec{pdim})

Evaluate the B-spline geometry map of the patch at the parametric coordinate `ξ`, i.e.
return the physical point ``\\boldsymbol{x}(\\boldsymbol{\\xi}) = \\sum_A N_A(\\boldsymbol{\\xi}) \\boldsymbol{X}_A``.

!!! warning
    The NURBS weights are currently not taken into account, so the result is only exact for
    B-spline patches (all weights equal to one). Furthermore, all basis functions are
    evaluated, so the function is slow for large patches.
"""
function eval_parametric_coordinate(mesh::NURBSMesh{pdim,sdim}, ξ::Vec{pdim}) where {pdim,sdim}
	bspline = BSplineBasis(mesh.knot_vectors, mesh.orders)
	@assert getnbasefunctions(bspline) == length(mesh.control_points)
	x = zero(Vec{sdim,Float64})
	for i in 1:getnbasefunctions(bspline)
		N = Ferrite.reference_shape_value(bspline, ξ, i)
		x += N*mesh.control_points[i]
	end
	return x
end

"""
    parent_to_parametric_map(mesh::NURBSMesh{pdim}, cellid::Int, xi::Vec{pdim})

Map the coordinate `xi` in the parent (reference) domain ``[-1, 1]^{pdim}`` of cell
`cellid` to the corresponding coordinate in the parametric (knot) domain of the patch. The
result is returned as a `Vector` with `pdim` entries.

Together with [`eval_parametric_coordinate`](@ref) this can be used to find the physical
position of a point on a cell.
"""
function parent_to_parametric_map(nurbsmesh::NURBSMesh{pdim}, cellid::Int, xi::Vec{pdim}) where {pdim}
	Ξ = nurbsmesh.knot_vectors
	_ni = nurbsmesh.INN[nurbsmesh.IEN[end,cellid],1:pdim]

	#Map to parametric domain from parent domain
	ξηζ = [0.5*((Ξ[d][_ni[d]+1] - Ξ[d][_ni[d]])*xi[d] + (Ξ[d][_ni[d]+1] + Ξ[d][_ni[d]])) for d in 1:pdim]
end

"""
	generate_nurbs_meshdata(orders, nbf)

Computes connectivity arrays for a nurbs patch, based on the order and number of basefunctions of the patch.
"""
function generate_nurbs_meshdata(orders::NTuple{dim,Int}, nbf::NTuple{dim,Int}) where dim

	nel = prod(nbf .- orders) #(n-p)*(m-q)*(l-r)
	nnp = prod(nbf) #n*m*l
	nen = prod(orders.+1) #(p+1)*(q+1)*(r+1)

	INN = zeros(Int, nnp, dim)
	IEN = zeros(Int, nen, nel)

	A = 0; e = 0
    dims = 1:dim

    for i in Tuple.(CartesianIndices(nbf))
        A += 1
        INN[A, dims] .= i
        if all(i .>= (orders.+1))
            e+=1
			for loc in Tuple.(CartesianIndices(orders.+1))
				loc = loc .- 1
				B = A
				b = 1
                for d in dim:-1:1
                    _d = dims[1:d-1]
					B -= loc[d] * prod(nbf[_d])
                    b += loc[d] * prod(orders[_d] .+ 1)
				end
                IEN[b,e] = B
            end
        end
	end
	IEN .= reverse(IEN, dims=1)

	return nel, nnp, nen, INN, IEN
end

"""
    knotinsertion!(knot_vectors, orders, control_points, weights, ξ; dir::Int)

h-refinement: insert the knot value `ξ` once in parametric direction `dir`, using Boehm's
algorithm [Piegl1997](@cite). The knot vector `knot_vectors[dir]`, the `control_points` and
the `weights` are updated in-place such that the geometry is unchanged. Returns `nothing`.

The arrays are those that define a [`NURBSMesh`](@ref) (and are e.g. obtained by copying the
fields of an existing mesh). After refinement, create a new `NURBSMesh` from the arrays so
that the element connectivity is updated.

See [`knotinsertion`](@ref) for a version that operates on a `NURBSMesh`.

# Examples
```julia
kv = (copy(mesh.knot_vectors[1]), copy(mesh.knot_vectors[2]))
cp = copy(mesh.control_points)
w = copy(mesh.weights)
for ξ in (-0.5, 0.0, 0.5)
    knotinsertion!(kv, mesh.orders, cp, w, ξ; dir = 1)
end
mesh = NURBSMesh(kv, mesh.orders, cp, w)
```
"""
function knotinsertion!(knot_vectors::NTuple{pdim,Vector{T}}, orders::NTuple{pdim,Int}, control_points::Vector{Vec{sdim,T}}, weights::Vector{T}, ξᴺ::T; dir::Int) where {pdim,sdim,T}

	Ξ = knot_vectors[dir]
	p = orders[dir]
	n = length(Ξ) - p - 1 
	n_other = length(control_points) ÷ n 

	#finds the index of the knot span where the new knot will be
	k = findfirst(>(ξᴺ), Ξ) - 1
	new_cps = Vector{Vec{sdim,T}}(undef, (n+1)*n_other)
	new_ws  = Vector{T}(undef, (n+1)*n_other)

	for r in 1:n_other
		stride   = dir == 1 ? 1 : n_other
		old_base = dir == 1 ? (r - 1) * n + 1 : r
		new_base = dir == 1 ? (r - 1) * (n + 1) + 1 : r

		#changes nothing and only adds to new knot vector when the new knot hasnt been reached yet
		for i in 1:k-p
			idx = new_base + (i - 1) * stride
			src = old_base + (i - 1) * stride
			new_cps[idx] = control_points[src]
			new_ws[idx]  = weights[src]
		end

		# inserts the new knot via NURBS weighted blend
		for j in 1:p
			i   = k - p + j
			α   = (ξᴺ - Ξ[i]) / (Ξ[i+p] - Ξ[i])
			src = old_base + (i - 1) * stride
			src_prev = src - stride
			wi  = α * weights[src] + (1 - α) * weights[src_prev]
			idx = new_base + (i - 1) * stride
			new_ws[idx]  = wi
			new_cps[idx] = (α * weights[src] * control_points[src] + (1 - α) * weights[src_prev] * control_points[src_prev]) / wi
		end

		# shifts the control points and weights after the insertion span
		for i in k+1:n+1
			idx = new_base + (i - 1) * stride
			src = old_base + (i - 2) * stride
			new_cps[idx] = control_points[src]
			new_ws[idx]  = weights[src]
		end
	end

	insert!(knot_vectors[dir], k+1, ξᴺ)
	copyto!(resize!(control_points, length(new_cps)), new_cps)
	copyto!(resize!(weights, length(new_ws)), new_ws)
	return nothing
end

"""
    orderelevation!(knot_vectors, orders, control_points, weights; dir::Int) -> orders

p-refinement: raise the polynomial degree by one in parametric direction `dir`, using
algorithm A5.9 of [Piegl1997](@cite). The knot vector `knot_vectors[dir]`, the
`control_points` and the `weights` are updated in-place such that the geometry is
unchanged. Since `orders` is a tuple, the new orders are returned.

After refinement, create a new [`NURBSMesh`](@ref) from the arrays. See
[`orderelevation`](@ref) for a version that operates on a `NURBSMesh`.
"""
function orderelevation!(knot_vectors::NTuple{pdim,Vector{T}}, orders::NTuple{pdim,Int}, control_points::Vector{Vec{sdim,T}}, weights::Vector{T}; dir::Int) where {pdim,sdim,T}

	Ξ = knot_vectors[dir]
	p = orders[dir]
	n_ctrl = length(Ξ) - p - 1
	n_other = length(control_points) ÷ n_ctrl
	t = 1
	ph = p + t
	ph2 = ph ÷ 2
	dim = sdim + 1
	m = length(Ξ) - 1

	#binomial coefficients for elevating each bezier segment by one degree
	bezalfs = zeros(T, ph + 1, p + 1)
	bezalfs[1, 1] = one(T)
	bezalfs[ph + 1, p + 1] = one(T)
	for i in 1:ph2
		inv = one(T) / T(binomial(ph, i))
		mpi = min(p, i)
		for j in max(0, i - t):mpi
			bezalfs[i + 1, j + 1] = inv * T(binomial(p, j)) * T(binomial(t, i - j))
		end
	end

	#uses symmetry to fill the rest of the coefficient table
	for i in (ph2 + 1):(ph - 1)
		mpi = min(p, i)
		for j in max(0, i - t):mpi
			bezalfs[i + 1, j + 1] = bezalfs[ph - i + 1, p - j + 1]
		end
	end

	#buffers reused while walking the knot vector
	bpts = Vector{Vector{T}}(undef, p + 1)
	ebpts = Vector{Vector{T}}(undef, ph + 1)
	Nextbpts = Vector{Vector{T}}(undef, max(p - 1, 0))
	alfs = Vector{T}(undef, max(p - 1, 0))
	Qw = Vector{Vector{T}}(undef, n_ctrl * 3)
	Uh = Vector{T}(undef, length(Ξ) * 3)
	new_cps = Vector{Vec{sdim,T}}(undef, 0)
	new_ws = Vector{T}(undef, 0)
	newΞ = Ξ

	#iterate over the adjacent directions to the elevated direction
	for r in 1:n_other
		stride = dir == 1 ? 1 : n_other
		old_base = dir == 1 ? (r - 1) * n_ctrl + 1 : r

		#load one row as homogeneous control points
		Pw = Vector{Vector{T}}(undef, n_ctrl)
		for i in 1:n_ctrl
			idx = old_base + (i - 1) * stride
			wi = weights[idx]
			hw = Vector{T}(undef, dim)
			for d in 1:sdim
				hw[d] = wi * control_points[idx][d]
			end
			hw[dim] = wi
			Pw[i] = hw
		end

		#set up for the first bezier span
		mh = ph
		kind = ph + 1
		rr = -1
		a = p
		b = p + 1
		cind = 2
		ua = Ξ[1]
		Qw[1] = Pw[1]
		Uh[1:(ph + 1)] .= ua
		bpts .= Pw[1:(p + 1)]

		#walk the knot vector span by span
		while b < m
			#find the next knot and its multiplicity
			i = b
			while b < m && Ξ[b + 2] == Ξ[b + 1]
				b += 1
			end
			mul = b - i + 1
			mh = mh + mul + t
			ub = Ξ[b + 1]
			oldr = rr
			rr = p - mul
			#decide which elevated points to keep from this span
			lbz = oldr > 0 ? (oldr + 2) ÷ 2 : 1
			rbz = rr > 0 ? ph - (rr + 1) ÷ 2 : ph

			#extract the span as a bezier via knot insertion
			if rr > 0
				numer = ub - ua
				for k in p:-1:(mul + 1)
					alfs[k - mul] = numer / (Ξ[a + k + 1] - ua)
				end
				for j in 1:rr
					save = rr - j
					s = mul + j
					for k in p:-1:s
						α = alfs[k - s + 1]
						bpts[k + 1] = α .* bpts[k + 1] .+ (1 - α) .* bpts[k]
					end
					Nextbpts[save + 1] = bpts[p + 1]
				end
			end

			#elevate the bezier segment by one degree
			for ii in lbz:ph
				ebpts[ii + 1] = zeros(T, dim)
				mpi = min(p, ii)
				for j in max(0, ii - t):mpi
					ebpts[ii + 1] .+= bezalfs[ii + 1, j + 1] .* bpts[j + 1]
				end
			end

			#remove the knot when the previous span had repeated knots
			if oldr > 1
				first = kind - 2
				last = kind
				den = ub - ua
				bet = (ub - Uh[kind]) / den
				for tr in 1:(oldr - 1)
					ii = first
					jj = last
					kj = jj - kind + 1
					while jj - ii > tr
						if ii < cind
							alf = (ub - Uh[ii + 1]) / (ua - Uh[ii + 1])
							Qw[ii + 1] = alf .* Qw[ii + 1] .+ (1 - alf) .* Qw[ii]
						end
						if jj >= lbz
							if jj - tr <= kind - ph + oldr
								gam = (ub - Uh[jj - tr + 1]) / den
								ebpts[kj + 1] = gam .* ebpts[kj + 1] .+ (1 - gam) .* ebpts[kj + 2]
							else
								ebpts[kj + 1] = bet .* ebpts[kj + 1] .+ (1 - bet) .* ebpts[kj + 2]
							end
						end
						ii += 1
						jj -= 1
						kj -= 1
					end
					first -= 1
					last += 1
				end
			end

			#load the knot ua into the new knot vector
			if a != p
				cnt = ph - oldr
				Uh[(kind + 1):(kind + cnt)] .= ua
				kind += cnt
			end
			#load the new control points from this span
			for j in lbz:rbz
				Qw[cind] = ebpts[j + 1]
				cind += 1
			end
			if b < m
				#load control points for the next span
				rr > 0 && (bpts[1:rr] .= Nextbpts[1:rr])
				bpts[(rr + 1):(p + 1)] .= Pw[(b - p + rr + 1):(b + 1)]
				a = b
				b += 1
				ua = ub
			else
				#load the end knot into the new knot vector
				Uh[(kind + 1):(kind + ph + 1)] .= ub
			end
		end

		nh = mh - ph - 1
		count_cpts = nh + 1
		if r == 1
			newΞ = Uh[1:(count_cpts + ph + 1)]
			resize!(new_cps, count_cpts * n_other)
			resize!(new_ws, count_cpts * n_other)
		end
		
		#convert back from homogeneous and store this row
		new_base = dir == 1 ? (r - 1) * count_cpts + 1 : r
		for i in 1:count_cpts
			hw = Qw[i]
			wi = hw[dim]
			cp = Vec(ntuple(d -> hw[d] / wi, sdim))
			idx = new_base + (i - 1) * stride
			new_cps[idx] = cp
			new_ws[idx] = wi
		end
	end

	#updates the elevated knots, cps, and weights for the final elevated ones
	copyto!(resize!(knot_vectors[dir], length(newΞ)), newΞ)
	copyto!(resize!(control_points, length(new_cps)), new_cps)
	copyto!(resize!(weights, length(new_ws)), new_ws)
	return ntuple(i -> i == dir ? orders[i] + 1 : orders[i], pdim)
end

"""
    smoothnesselevation!(knot_vectors, orders, control_points, weights, new_knots; dir::Int) -> orders

k-refinement: raise the degree by one in parametric direction `dir` with
[`orderelevation!`](@ref), then insert each value of `new_knots` with
[`knotinsertion!`](@ref). The new knots have continuity ``C^{p}``, where ``p`` is the
original degree. The arrays are updated in-place, and the new orders are returned.

After refinement, create a new [`NURBSMesh`](@ref) from the arrays. See
[`smoothnesselevation`](@ref) for a version that operates on a `NURBSMesh`.
"""
function smoothnesselevation!(knot_vectors::NTuple{pdim,Vector{T}}, orders::NTuple{pdim,Int}, control_points::Vector{Vec{sdim,T}}, weights::Vector{T}, new_knots::AbstractVector{T}; dir::Int) where {pdim,sdim,T}
	orders = orderelevation!(knot_vectors, orders, control_points, weights; dir=dir)
	for ξᴺ in new_knots
		knotinsertion!(knot_vectors, orders, control_points, weights, ξᴺ; dir=dir)
	end
	return orders
end
