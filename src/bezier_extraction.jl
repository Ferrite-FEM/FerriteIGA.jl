
function bezier_extraction_to_vectors(Ce::AbstractVector{<:AbstractMatrix})
    T = Float64
    nbe = length(Ce)

    Cvecs = [Vector{SparseArrays.SparseVector{T,Int}}() for _ in 1:nbe]
    for ie in 1:nbe
        cv = bezier_extraction_to_vector(Ce[ie])
        Cvecs[ie] = cv
    end
    return Cvecs
end

function bezier_extraction_to_vector(Ce::AbstractMatrix{T}) where T

    Cvecs = Vector{SparseArrays.SparseVector{T,Int}}()
    for r in 1:size(Ce,1)
        ce = Ce[r,:]
        push!(Cvecs, SparseArrays.sparsevec(ce))
    end

    return Cvecs
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
	compute_bezier_points(bezier_points::AbstractVector{T2}, Ce::BezierExtractionOperator{T}, control_points::AbstractVector{T2}; dim::Int=1)

Given a BezierExtractionOperator and control points for a cell, compute the bezier controlpoints.
"""
function compute_bezier_points!(bezier_points::Vector{T2}, Ce::BezierExtractionOperator, control_points::AbstractVector{T2}; dim::Int=1) where {T2}

	n_points = length(first(Ce))
	@boundscheck (length(control_points) == length(Ce)*dim)
	@boundscheck (length(bezier_points) == n_points*dim)

	Base.@inbounds for i in 1:length(bezier_points)
		bezier_points[i] = zero(T2)
	end

	for ic in 1:length(control_points)÷dim
		ce = Ce[ic]

		for i in 1:n_points
			for d in 1:dim
				bezier_points[(i-1)*dim + d] += ce[i] * control_points[(ic-1)*dim + d]
			end
		end
	end

	return nothing
end
#const bezier_extraction! = compute_bezier_points!

"""
	compute_bezier_points(Ce::BezierExtractionOperator{T}, control_points::AbstractVector{T2}; dim::Int=1)

Given a BezierExtractionOperator and control points for a cell, compute the bezier controlpoints.
"""
function compute_bezier_points(Ce::BezierExtractionOperator{T}, control_points::AbstractVector{T2}; dim::Int=1) where {T2,T}

	bezierpoints = zeros(T2, length(first(Ce))*dim)
	compute_bezier_points!(bezierpoints, Ce, control_points, dim=dim)

	return bezierpoints

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

	#The last C-matrix is not used
	#pop!(C)
	C = SparseArrays.sparse.(C[1:nb])
	return C, nb
end

function diagonal_beo(N::Int)
    beo = BezierExtractionOperator{Float64}(undef, N)
    for i in 1:N
        V = SparseArrays.sparsevec([i], [1.0], N)
        beo[i] = V
    end
    return beo
end

function combine_beo(a::BezierExtractionOperator{Float64}, b::BezierExtractionOperator{Float64})
    la = length(a)
    lb = length(b)
    ntotal = la+lb
    beo = BezierExtractionOperator{Float64}()
    for i in 1:la
        newrow = SparseArrays.sparsevec(a[i].nzind, a[i].nzval, ntotal)
        push!(beo, newrow)
    end
    for i in 1:lb
        newrow = SparseArrays.sparsevec(b[i].nzind .+ la, b[i].nzval, ntotal)
        push!(beo, newrow)
    end
    return beo
end