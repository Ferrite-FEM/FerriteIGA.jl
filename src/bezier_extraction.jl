"""
    transform_basis!(C::BezierExtractionOperator, Nout::AbstractMatrix, Nin::AbstractMatrix, [w=nothing]; stride::Int=1)

Transform basis function values or their derivatives from the Bernstein basis (`Nin`) to the B-spline or NURBS basis (`Nout`) using the extraction operator `C`.

When the optional NURBS weights `w` are provided, the transformation incorporates them element-wise as ``N_i^{\\text{out}} = w_i \\sum_j C_{i,j} N_j^{\\text{in}}``. 

The `stride` argument allows the transformation of vector-valued basis functions. 
"""
function transform_basis!(
    C::BezierExtractionOperator,
    Nout::AbstractMatrix{T},
    Nin::AbstractMatrix{T},
    w = nothing;
    stride::Int = 1,
) where {T}
    #The number of spline functions (rows of `Nout`) can differ from the number of 
    #Bernstein functions (rows of `Nin`), e.g. for T-splines.
    @assert size(Nout, 2) == size(Nin, 2)
    @assert size(Nout, 1) == length(C.C) * stride
    @assert w === nothing || length(w) == length(C.C)

    #Vectorized interpolations (ip^vdim): transform the scalar basis once instead of every vector entry
    if stride > 1 && _is_vectorized_entry(T, stride)
        return _transform_vectorized_basis!(C, Nout, Nin, w, stride)
    end

    @inbounds for col in axes(Nout, 2)
        for i in eachindex(C.C)
            Ci = C.C[i]
            wi = w === nothing ? one(T) : w[i]
            nzind = Ci.nzind
            nzval = Ci.nzval
            
            for d in 1:stride
                oi = (i - 1) * stride + d
                s = zero(T)
                for k in eachindex(nzind)
                    j = nzind[k]
                    ii = (j - 1) * stride + d
                    s += nzval[k] * Nin[ii, col]
                end
                Nout[oi, col] = wi * s
            end
        end
    end
    return nothing
end

# Vectorized interpolation `ip^vdim`: the basis function `(j-1)*vdim + d` is `e_d Nⱼ` (values), `e_d ⊗ ∇Nⱼ`
# (gradients) or `e_d ⊗ ∇∇Nⱼ` (hessians). The scalar part of basis function j is component/row 1 of entry
# `(j-1)*vdim + 1`. These are only used when the entry type matches the stride (vdim == spatial dimension).
_is_vectorized_entry(::Type{<:Vec{dim}}, stride) where {dim} = dim == stride
_is_vectorized_entry(::Type{<:Tensor{2,dim}}, stride) where {dim} = dim == stride
_is_vectorized_entry(::Type{<:Tensor{3,dim}}, stride) where {dim} = dim == stride
_is_vectorized_entry(::Type, stride) = false

@inline _scalar_part(v::Vec) = v[1]
@inline _scalar_part(t::Tensor{2,dim}) where {dim} = Vec{dim}(k -> t[1, k])
@inline _scalar_part(t::Tensor{3,dim}) where {dim} = Tensor{2,dim}((k, l) -> t[1, k, l])

@inline _vectorized_entry(::Type{<:Vec{dim}}, s::Number, d::Int) where {dim} = Vec{dim}(k -> k == d ? s : zero(s))
@inline _vectorized_entry(::Type{<:Tensor{2,dim}}, g::Vec{dim,T}, d::Int) where {dim,T} = Tensor{2,dim}((k, l) -> k == d ? g[l] : zero(T))
@inline _vectorized_entry(::Type{<:Tensor{3,dim}}, h::Tensor{2,dim,T}, d::Int) where {dim,T} = Tensor{3,dim}((k, l, m) -> k == d ? h[l, m] : zero(T))

function _transform_vectorized_basis!(C::BezierExtractionOperator, Nout::AbstractMatrix{T}, Nin::AbstractMatrix{T}, w, stride::Int) where {T}
    @inbounds for col in axes(Nout, 2)
        for i in eachindex(C.C)
            Ci = C.C[i]
            nzind = Ci.nzind
            nzval = Ci.nzval
            s = zero(_scalar_part(Nin[1, col]))
            for k in eachindex(nzind)
                s += nzval[k] * _scalar_part(Nin[(nzind[k] - 1) * stride + 1, col])
            end
            if w !== nothing
                s *= w[i]
            end
            for d in 1:stride
                Nout[(i - 1) * stride + d, col] = _vectorized_entry(T, s, d)
            end
        end
    end
    return nothing
end

"""
    transform_coords!(xb::AbstractVector, wb::AbstractVector, C::BezierExtractionOperator, x::AbstractVector, w::AbstractVector)

Compute the local Bézier control points `xb` and corresponding weights `wb` for an element, given the global B-spline/NURBS control coordinates `x` and weights `w`.

The transformation projects the global control points into homogeneous space, applies the extraction operator `C`, and projects the results back into Cartesian coordinates. The operation mutates `xb` and `wb` in-place according to:

``\\mathbf{x}_b^{(j)} = \\frac{1}{w_b^{(j)}} \\sum_i C_{i,j} w_i \\mathbf{x}_i``

``w_b^{(j)} = \\sum_i C_{i,j} w_i``
"""
function transform_coords!(
    xb::AbstractVector{Vec{dim,T}},
    wb::AbstractVector{T},
    C::BezierExtractionOperator,
    x::AbstractVector{Vec{dim,T}},
    w::AbstractVector{T},
) where {dim,T}

	C = C.C
    n = length(C)

	#`x` and `w` are the spline control points of the cell, and `xb` and `wb` are the Bézier control points.
	@assert length(x) == n
    @assert length(xb) == length(wb)
    @assert length(w) == length(x)
    fill!(xb, zero(Vec{dim,T}))
    fill!(wb, zero(T))

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

    return xb, wb
end

function bezier_extraction_to_vectors(Ce::AbstractVector{<:AbstractMatrix})
    T = Float64
    nbeos = length(Ce)

    beos = BezierExtractionOperator{T}[]
    for ie in 1:nbeos
        beo = bezier_extraction_to_vector(Ce[ie])
        push!(beos, beo)
    end
    return beos
end

function bezier_extraction_to_vector(Ce::AbstractMatrix{T}) where T

    Cvecs = Vector{SparseArrays.SparseVector{T,Int}}()
    for r in 1:size(Ce,1)
        ce = Ce[r,:]
        push!(Cvecs, SparseArrays.sparsevec(ce))
    end

    return BezierExtractionOperator(Cvecs)
end

function beo2matrix(m::BezierExtractionOperator{T}) where T

    m2 = zeros(T, length(m), length(first(m)))
    for r in 1:length(m)
        for c in 1:length(m[r])
            m2[r,c] = m[r][c]
        end
    end
    return m2
end

"""
	compute_bezier_extraction_operators2(orders::NTuple{pdim,Int}, knots::NTuple{pdim,Vector{T}})

Computes the bezier extraction operator in each parametric direction, and uses the kron operator to combine them.
"""
@generated function compute_bezier_extraction_operators(orders::NTuple{pdim,Int}, knots::NTuple{pdim,Vector{T}}) where {pdim,T} 
	quote
		#Get bezier extraction vector in each dimension
		#
		Ce = Vector{SparseArrays.SparseMatrixCSC{Float64,Int64}}[]
		nels = Int[]
		for d in 1:pdim
			_Ce, _nel = _compute_bezier_extraction_operators(orders[d], knots[d])
			push!(Ce, _Ce)
			push!(nels, _nel)
		end
		
		#Tensor prodcut of the bezier extraction operators
		#
		C = Vector{eltype(first(Ce))}()
		Base.Cartesian.@nloops $pdim i (d)->1:nels[d] begin
			#kron not defined for 1d, so special case for pdim==1
			if $pdim == 1
				_C = Ce[1][i_1]
			else
				_C = Base.Cartesian.@ncall $pdim kron (d)->Ce[$pdim-d+1][i_{$pdim-d+1}]
			end
			push!(C, _C)
		end

		#Reorder
		#
		@assert allequal(orders)
		ordering = _bernstein_ordering(Bernstein{RefHypercube{pdim}, orders[1]}())
		nel = prod(nels)
		C_reorder = Vector{eltype(first(Ce))}()
		for i in 1:nel
			push!(C_reorder, C[i][ordering,ordering])
		end

		return C_reorder, nel
	end
end

function _compute_bezier_extraction_operators(p::Int, knot::Vector{T}) where T
	a = p+1
	b = a+1
	nb = 1
	m = length(knot)
	C = [Matrix(Diagonal(ones(T,p+1)))]

	while b < m
		push!(C, Matrix(Diagonal(ones(T,p+1))))
		i = b
		while b<m && knot[b+1]==knot[b]
			 b+=1;
		end;
		mult = b-i+1

		if mult < p + 1
			
			α = zeros(T,p)
			for j in p:-1:(mult+1)
				α[j-mult] = (knot[b]-knot[a])/(knot[a+j]-knot[a])
			end
			r = p-mult

			for j in 1:r
				save = r-j+1
				s = mult+j

				for k in (p+1):-1:(s+1)
					C[nb][:,k] = α[k-s]*C[nb][:,k] + (1.0-α[k-s])*C[nb][:,k-1]
				end
				if b<m
					C[nb+1][save:(j+save),save] = C[nb][(p-j+1):(p+1), p+1]
				end
			end
			nb += 1
			if b<m
				a=b
				b+=1
			end
		end
	end

	C = SparseArrays.sparse.(C[1:nb])
	return C, nb
end
"""
	_insert_knot(U, c, p, u)

Insert the knot `u` once in the knot vector `U` of a univariate B-spline function with coefficients `c` (Boehm's algorithm).
"""
function _insert_knot(U::Vector{T}, c::Vector{T}, p::Int, u::T) where T
	k = searchsortedlast(U, u) #U[k] <= u < U[k+1]
	n = length(c)
	cnew = zeros(T, n+1)
	for i in 1:n+1
		if i <= k-p
			cnew[i] = c[i]
		elseif i <= k
			α = (u - U[i])/(U[i+p] - U[i])
			cnew[i] = α*c[i] + (1-α)*c[i-1]
		else
			cnew[i] = c[i-1]
		end
	end
	Unew = copy(U)
	insert!(Unew, k+1, u)
	return Unew, cnew
end

"""
	_bspline_bezier_coefficients(Ξ, p, a, b)

Coefficients of the univariate B-spline basis function with the local knot vector `Ξ` (`p+2` knots),
restricted to the interval `[a,b]`, in the Bernstein basis of degree `p` on `[a,b]`. The interval must 
not contain knots of `Ξ` in its interior.

This is univariate Bézier extraction of a single (T-spline) basis function (Scott et al. 2011): the local knot vector is 
extended to an open knot vector, and the knots `a` and `b` are inserted until they have multiplicity `p`.
"""
function _bspline_bezier_coefficients(Ξ::AbstractVector{T}, p::Int, a::T, b::T) where T
	@assert length(Ξ) == p+2
	@assert Ξ[1] <= a < b <= Ξ[end]
	@assert !any(ξ -> a < ξ < b, Ξ)

	#Extended knot vector, where the first and last knot have multiplicity p+1
	nt = p + 1 - count(==(Ξ[1]), Ξ)
	ne = p + 1 - count(==(Ξ[end]), Ξ)
	U = T[fill(Ξ[1], nt); Ξ; fill(Ξ[end], ne)]
	c = zeros(T, length(U) - p - 1)
	c[nt+1] = one(T)

	for u in (a, b)
		for _ in (count(==(u), U) + 1):p
			U, c = _insert_knot(U, c, p, u)
		end
	end

	k = findlast(==(a), U)
	@assert U[k+1] == b
	return c[(k-p):k]
end

"""
	compute_bezier_extraction_operators(mesh::TSplineMesh)

Compute the Bézier extraction operator of each element of the T-spline. The rows of the operator for element `e`
correspond to the basis functions `mesh.IEN[e]`, and the columns to the Bernstein basis functions 
(ordered as in `IGAInterpolation`).
"""
function compute_bezier_extraction_operators(mesh::TSplineMesh{sdim,T}) where {sdim,T}
	p = mesh.order
	ordering = _bernstein_ordering(Bernstein{RefQuadrilateral,p}())
	beos = BezierExtractionOperator{T}[]
	for (e, (ξa, ξb, ηa, ηb)) in enumerate(mesh.elements)
		rows = SparseVector{T,Int}[]
		for A in mesh.IEN[e]
			Ξ, Π = mesh.local_knots[A]
			c1 = _bspline_bezier_coefficients(Ξ, p, ξa, ξb)
			c2 = _bspline_bezier_coefficients(Π, p, ηa, ηb)
			push!(rows, sparsevec(kron(c2, c1)[ordering]))
		end
		push!(beos, BezierExtractionOperator(rows))
	end
	return beos
end
