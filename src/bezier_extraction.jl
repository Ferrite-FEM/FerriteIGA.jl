
function transform_basis!(
    C::BezierExtractionOperator,
    Nout::AbstractMatrix{T},
    Nin::AbstractMatrix{T},
    w = nothing;
    stride::Int = 1,
) where {T}
    @assert size(Nout) == size(Nin)
    @assert size(Nout, 1) == length(C.C) * stride
    @assert w === nothing || length(w) == length(C.C)

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

function transform_coords!(
    xb::AbstractVector{Vec{dim,T}},
    wb::AbstractVector{T},
    C::BezierExtractionOperator,
    x::AbstractVector{Vec{dim,T}},
    w::AbstractVector{T},
) where {dim,T}

	C = C.C
    n = length(C)

	@assert length(xb) == n
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