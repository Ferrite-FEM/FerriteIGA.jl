export NURBSMesh, parent_to_parametric_map, eval_parametric_coordinate
export knotinsertion!, orderelevation!, smoothnesselevation!

"""
	NURBSMesh(knot_vectors, orders, control_points, weights=ones(...))

A NURBS or B-spline patch. `knot_vectors` and `orders` give the basis in each direction.
`control_points` give the shape. `weights` are all ones for a B-spline.
"""
struct NURBSMesh{pdim,sdim,T} #<: Ferrite.AbstractGrid
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
const getncontrolponits = Ferrite.getnnodes

function Ferrite.getcoordinates(mesh::NURBSMesh, ie::Int)
	return mesh.control_points[mesh.IEN[:,ie]]
end

"""
	eval_parametric_coordinate(mesh::NURBSMesh, ξ)

Physical point corresponding to the parameter `ξ` on the patch.

TODO: This function is currently very in-effecient for large domain.
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

	#Possible faster algorithm:
	#=
	knot_spans = ntuple(pdim) do i
		_find_span(n[i], p[i], ξ[i], Ξ[i])
	end

	shape_values = ntuple(pdim) do i
		_eval_nonzero_bspline_values!(N[i], first(knot_span[i]), orders[i], ξ[i], knot_vectors[i] )
	end

	x = zero(Vec{sdim,Float64})
	for k in knot_spans[3]
		Nk = N[k]
		for j in knot_spans[2]
			Nj = N[k]
			N = Nj*Nk
			for i in knot_spans[1]
				N *= N[i]
				index = CartesianIndex(i,j,k)
				x += N * control_points[index]
			end
		end
	end =#
	
end

"""
	parent_to_parametric_map(nurbsmesh::NURBSMesh{pdim,sdim}, cellid::Int, xi::Vec{pdim})

Given a coordinate for a cell in the parent domain `xi`, this functions returns the coordinate in 
the parametric domain.
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
	knotinsertion!(knot_vectors, orders, control_points, weights, ξ; dir)

h-refinement. Insert the knot `ξ` once in parametric direction `dir` (Boehm).
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
	orderelevation!(knot_vectors, orders, control_points, weights; dir)

p-refinement. Raise the polynomial degree by one in direction `dir`.
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
	smoothnesselevation!(knot_vectors, orders, control_points, weights, new_knots; dir)

k-refinement. Raise the degree by one in direction `dir`, then insert each knot in `new_knots`.
This is a pointwise smoothness elevation. 
"""
function smoothnesselevation!(knot_vectors::NTuple{pdim,Vector{T}}, orders::NTuple{pdim,Int}, control_points::Vector{Vec{sdim,T}}, weights::Vector{T}, new_knots::AbstractVector{T}; dir::Int) where {pdim,sdim,T}
	orders = orderelevation!(knot_vectors, orders, control_points, weights; dir=dir)
	for ξᴺ in new_knots
		knotinsertion!(knot_vectors, orders, control_points, weights, ξᴺ; dir=dir)
	end
	return orders
end
