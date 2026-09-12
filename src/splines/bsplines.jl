"""
BSplineBasis{dim,T,order} <: Ferrite.ScalarInterpolation{RefHypercube{dim}, order}
    
    Not really used in FE-codes. We use it for testing in FerriteIGA.jl
"""
struct BSplineBasis{dim,T,order} <: Ferrite.ScalarInterpolation{RefHypercube{dim}, order}
	knot_vector::NTuple{dim,Vector{T}}

    function BSplineBasis(knots::NTuple{dim,Vector{T}}, order::NTuple{dim,Int}) where {dim,T} 
        @assert(length(order)==dim)
        @assert all( issorted.(knots) )
		return new{dim,T,Tuple(order)}(knots)
    end

    function BSplineBasis(knots::Vector{T}, order::Int) where {T} 
		return BSplineBasis((knots,), (order,))
    end
    
end

getnbasefunctions_dim(basis::BSplineBasis{dim,T,order}) where {dim,T,order} = ntuple(i -> length(basis.knot_vector[i]) - order[i] - 1, dim)
Ferrite.getnbasefunctions(basis::BSplineBasis{dim,T,order}) where {dim,T,order} = prod(getnbasefunctions_dim(basis))

function Ferrite.reference_shape_value(b::BSplineBasis{dim,T,order}, xi::Vec{dim,T2}, i::Int) where {dim,T,T2,order}
    @assert( i <= Ferrite.getnbasefunctions(b))

    _n = getnbasefunctions_dim(b)
    
    indecies = CartesianIndices(_n)[i]
    val = one(T2)
    for i in 1:dim
        val *= _bspline_basis_value_alg1(order[i], b.knot_vector[i], indecies[i], xi[i])
    end
    return val
end

"""
Evaluate a single B-spline basis function
"""
function _bspline_basis_value_alg1(p::Int, knot::AbstractVector{<:Real}, i::Int, xi::T) where {T}

    if p > 0
        # First term
        denom1 = knot[i+p] - knot[i]
        term1 = if iszero(denom1)
            zero(T)
        else
            (xi - knot[i]) * _bspline_basis_value_alg1(p-1, knot, i, xi) / denom1
        end

        # Second term
        denom2 = knot[i+p+1] - knot[i+1]
        term2 = if iszero(denom2)
            zero(T)
        else
            (knot[i+p+1] - xi) * _bspline_basis_value_alg1(p-1, knot, i+1, xi) / denom2
        end

        return term1 + term2
    else
        # Base case (indicator function)
        if knot[i+1] < knot[end]
            return (knot[i] <= xi < knot[i+1]) ? one(T) : zero(T)
        else
            return (knot[i] <= xi) ? one(T) : zero(T)
        end
    end
end

"""
Algorithm for calculating one basis functions value
From NURBS-book, alg2.4
"""
function _bspline_basis_value_alg2(p::Int,U::Vector{T},i::Int,u::T2) where {T,T2}
	i -=1
    m = length(U)-1
    N = zeros(T2, 10)
    if ((i==0 && u==U[0+1] || (i==m-p-1 && u==U[end])))
        return 1.0
    end
    if (u < U[i+1] || u>=U[i+p+1+1])
        return 0.0
    end
    for j in 1:(p+1)
        if u >= U[i+j] && u < U[i+j+1]
            N[j] = 1.0
        else
            N[j] = 0.0
        end
    end
    for k in (1:p) .+1
        if N[0+1] == 0.0
            saved = 0.0
        else
            saved = ((u-U[i+1])*N[0+1])/(U[i+k]-U[i+1])
        end
        for j in (0:(p-k+1)) .+1
            Uleft = U[i+j+1]
            Uright = U[i+j+k+1-1]
            if N[j+1] == 0.0
                N[j] = saved;
                saved = 0.0
            else
                temp = N[j+1]/(Uright-Uleft)
                N[j] = saved + (Uright-u)*temp
                saved = (u-Uleft)*temp
            end
        end
    end
    return N[1]
end

"""
Algorithm that evaluates the knotspan in which u (xi) lies.
From NURBS-book alg2.1
"""
function _find_span(n::Int ,p::Int,u::T,U::Vector{T}) where T
    @assert(!(u<U[1] || u>U[end]))

    if u==U[end]
        span = length(U)-1
        while U[span]==U[span+1]
            span+=-1
        end
    else
        low = 0
        high = length(U)
        mid = round(Int64,(low + high)/2)

        while u < U[mid] || u >= U[mid + 1]
	        if u < U[mid]; high = mid;
	        else; low = mid;
	        end
	        mid = round(Int64,(low + high)/2)
        end
        span = round(Int64,(low + high)/2)
    end
    return span
end

