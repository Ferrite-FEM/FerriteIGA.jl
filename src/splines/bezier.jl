export Bernstein

"""
    Bernstein{shape, order}()

Bernstein polynomials of degree `order` on the reference `shape`.
`IGAInterpolation` uses this basis on each cell after Bézier extraction.
"""
struct Bernstein{shape, order} <: Ferrite.ScalarInterpolation{shape, order}
    function Bernstein{shape,order}() where {rdim, shape<:RefHypercube{rdim}, order} 
        @assert order isa Int
        return new{shape,order}()
    end
end

Ferrite.adjust_dofs_during_distribution(::Bernstein) = true
Ferrite.adjust_dofs_during_distribution(::Bernstein{<:Any, 2}) = false
Ferrite.adjust_dofs_during_distribution(::Bernstein{<:Any, 1}) = false

Ferrite.vertexdof_indices(ip::Bernstein{refshape,order}) where {refshape,order} = _compute_vertexdof_indices(ip)
Ferrite.edgedof_indices(ip::Bernstein{refshape,order}) where {refshape,order} = _compute_edgedof_indices(ip)
Ferrite.facedof_indices(ip::Bernstein{refshape,order}) where {refshape,order} = _compute_facedof_indices(ip)

# # #
#   Bernstein line, order 2
# # #
Ferrite.getnbasefunctions(::Bernstein{RefLine,2}) = 3

Ferrite.vertexdof_indices(::Bernstein{RefLine,2}) = ((1,),(2,))
Ferrite.edgedof_indices(::Bernstein{RefLine,2}) = ((1,2,3),)
Ferrite.facedof_indices(::Bernstein{RefLine,2}) = ()

function Ferrite.reference_shape_value(ip::Bernstein{RefLine,2}, _ξ::Vec{1}, i::Int)
    ξ = 0.5*(_ξ[1] + 1.0)
    i == 1 && return (1-ξ)^2
    i == 2 && return ξ^2
    i == 3 && return 2ξ*(1 - ξ) 
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end


# # #
#   Bernstein Quadrilateral, order 2
# # #
Ferrite.getnbasefunctions(::Bernstein{RefQuadrilateral,2}) = 9
Ferrite.vertexdof_indices(::Bernstein{RefQuadrilateral,2}) = ((1,),(2,),(3,),(4,))
Ferrite.edgedof_indices(::Bernstein{RefQuadrilateral,2}) = ((1,2, 5), (2,3, 6), (3,4, 7), (4,1, 8))
Ferrite.edgedof_interior_indices(::Bernstein{RefQuadrilateral,2}) = ((5,), (6,), (7,), (8,))
Ferrite.facedof_indices(::Bernstein{RefQuadrilateral,2}) = ((1,2,3,4,5,6,7,8,9),)
Ferrite.facedof_interior_indices(::Bernstein{RefQuadrilateral,2}) = (9,)

function Ferrite.reference_shape_value(ip::Bernstein{RefQuadrilateral,2}, _ξ::Vec{2}, i::Int)
    ξ, η = _ξ
    i == 1 && return 0.0625((1 - η)^2)*((1 - ξ)^2)
    i == 2 && return 0.0625((1 + ξ)^2)*((1 - η)^2)
    i == 3 && return 0.0625((1 + η)^2)*((1 + ξ)^2)
    i == 4 && return 0.0625((1 + η)^2)*((1 - ξ)^2)
    i == 5 && return 0.125(1 + ξ)*((1 - η)^2)*(1 - ξ)
    i == 6 && return 0.125(1 + η)*((1 + ξ)^2)*(1 - η)
    i == 7 && return 0.125(1 + ξ)*((1 + η)^2)*(1 - ξ)
    i == 8 && return 0.125(1 + η)*((1 - ξ)^2)*(1 - η)
    i == 9 && return 0.25(1 + η)*(1 + ξ)*(1 - η)*(1 - ξ)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

# # #
#   Bernstein Hexahedron, order 2
# # #
Ferrite.getnbasefunctions(::Bernstein{RefHexahedron,2}) = 27
Ferrite.vertexdof_indices(::Bernstein{RefHexahedron,2}) = (
    (1,),(2,),(3,),(4,),(5,),(6,),(7,),(8,)
)
Ferrite.facedof_indices(::Bernstein{RefHexahedron,2}) = (
    (1,4,3,2, 12,11,10,9, 21),
    (1,2,6,5, 9,18,13,17, 22),
    (2,3,7,6, 10,19,14,18, 23),
    (3,4,8,7, 11,20,15,19, 24),
    (1,5,8,4, 17,16,20,12, 25),
    (5,6,7,8, 13,14,15,16, 26),
)
Ferrite.facedof_interior_indices(::Bernstein{RefHexahedron,2}) = (
    (21,), (22,), (23,), (24,), (25,), (26,),
)
Ferrite.edgedof_indices(::Bernstein{RefHexahedron,2}) = (
    (1,2, 9),
    (2,3, 10),
    (3,4, 11),
    (4,1, 12),
    (5,6, 13),
    (6,7, 14),
    (7,8, 15),
    (8,5, 16),
    (1,5, 17),
    (2,6, 18),
    (3,7, 19),
    (4,8, 20),
)
Ferrite.edgedof_interior_indices(::Bernstein{RefHexahedron,2}) = (
    (9,), (10,), (11,), (12,), (13,), (14,), (15,), (16,), (17), (18,), (19,), (20,)
)
Ferrite.volumedof_interior_indices(::Bernstein{RefHexahedron,2}) = (27,)

function Ferrite.reference_shape_value(ip::Bernstein{RefHexahedron,2}, _ξ::Vec{3}, i::Int)
    ξ, η, ζ = _ξ
    i == 1 && return 0.015625((1 - ζ)^2)*((1 - η)^2)*((1 - ξ)^2)
    i == 2 && return 0.015625((1 + ξ)^2)*((1 - ζ)^2)*((1 - η)^2)
    i == 3 && return 0.015625((1 + η)^2)*((1 + ξ)^2)*((1 - ζ)^2)
    i == 4 && return 0.015625((1 + η)^2)*((1 - ζ)^2)*((1 - ξ)^2)
    i == 5 && return 0.015625((1 + ζ)^2)*((1 - η)^2)*((1 - ξ)^2)
    i == 6 && return 0.015625((1 + ζ)^2)*((1 + ξ)^2)*((1 - η)^2)
    i == 7 && return 0.015625((1 + ζ)^2)*((1 + η)^2)*((1 + ξ)^2)
    i == 8 && return 0.015625((1 + ζ)^2)*((1 + η)^2)*((1 - ξ)^2)
    i == 9 && return 0.03125(1 + ξ)*((1 - ζ)^2)*((1 - η)^2)*(1 - ξ)
    i == 10 && return 0.03125(1 + η)*((1 + ξ)^2)*((1 - ζ)^2)*(1 - η)
    i == 11 && return 0.03125(1 + ξ)*((1 + η)^2)*((1 - ζ)^2)*(1 - ξ)
    i == 12 && return 0.03125(1 + η)*((1 - ζ)^2)*((1 - ξ)^2)*(1 - η)
    i == 13 && return 0.03125(1 + ξ)*((1 + ζ)^2)*((1 - η)^2)*(1 - ξ)
    i == 14 && return 0.03125(1 + η)*((1 + ζ)^2)*((1 + ξ)^2)*(1 - η)
    i == 15 && return 0.03125(1 + ξ)*((1 + ζ)^2)*((1 + η)^2)*(1 - ξ)
    i == 16 && return 0.03125(1 + η)*((1 + ζ)^2)*((1 - ξ)^2)*(1 - η)
    i == 17 && return 0.03125(1 + ζ)*((1 - η)^2)*((1 - ξ)^2)*(1 - ζ)
    i == 18 && return 0.03125(1 + ζ)*((1 + ξ)^2)*((1 - η)^2)*(1 - ζ)
    i == 19 && return 0.03125(1 + ζ)*((1 + η)^2)*((1 + ξ)^2)*(1 - ζ)
    i == 20 && return 0.03125(1 + ζ)*((1 + η)^2)*((1 - ξ)^2)*(1 - ζ)
    i == 21 && return 0.0625(1 + η)*(1 + ξ)*((1 - ζ)^2)*(1 - η)*(1 - ξ)
    i == 22 && return 0.0625(1 + ζ)*(1 + ξ)*((1 - η)^2)*(1 - ζ)*(1 - ξ)
    i == 23 && return 0.0625(1 + ζ)*(1 + η)*((1 + ξ)^2)*(1 - ζ)*(1 - η)
    i == 24 && return 0.0625(1 + ζ)*(1 + ξ)*((1 + η)^2)*(1 - ζ)*(1 - ξ)
    i == 25 && return 0.0625(1 + ζ)*(1 + η)*((1 - ξ)^2)*(1 - ζ)*(1 - η)
    i == 26 && return 0.0625(1 + η)*(1 + ξ)*((1 + ζ)^2)*(1 - η)*(1 - ξ)
    i == 27 && return 0.125(1 + ζ)*(1 + η)*(1 + ξ)*(1 - ζ)*(1 - η)*(1 - ξ)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

function Ferrite.reference_shape_value(ip::Bernstein{refshape,order}, _ξ::Vec{dim}, i::Int) where {dim,refshape<:Ferrite.AbstractRefShape{dim},order}
    _compute_bezier_reference_shape_value(ip,_ξ, i)
end

function _compute_bezier_reference_shape_value(ip::Bernstein{shape,order}, ξ::Vec{dim,T}, i::Int) where {dim,shape<:Ferrite.AbstractRefShape{dim},order,T} 
    _n = ntuple(i->order+1, dim)
    ordering = _bernstein_ordering(ip)
    basefunction_indeces = CartesianIndices(_n)[ordering[i]]

    val = one(T)
    for i in 1:dim
        val *= _bernstein_basis_recursive(order, basefunction_indeces[i], ξ[i])
    end
    return val
end

function _compute_vertexdof_indices(::Bernstein{RefLine,order}) where order
    ((1,),(2,))
end

function _compute_edgedof_indices(ip::Bernstein{RefLine,order}) where {order}
    ind = ntuple(i-> i, order+1)
    return (ind,)
end

function _compute_facedof_indices(ip::Bernstein{RefLine,order}) where {order}
    ()
end

function _compute_vertexdof_indices(::Bernstein{RefQuadrilateral,order}) where order
    ((1,),(2,),(3,),(4,),)
end

function _compute_facedof_indices(ip::Bernstein{RefQuadrilateral,order}) where {order}
    ind = ntuple(i-> i, (order+1)^2)
    return (ind,)
end

function _compute_edgedof_indices(ip::Bernstein{RefQuadrilateral,order}) where {order}
    faces = Tuple[]
    orders = (order, order)

    ind = reshape([findfirst(i->i==j, _bernstein_ordering(ip)) for j in 1:prod(orders.+1)], (orders.+1)...)
    
    #Order for 1d interpolation.
    order_1d1 = _bernstein_ordering(Bernstein{RefLine,order}())
    order_1d2 = _bernstein_ordering(Bernstein{RefLine,order}())

    push!(faces, Tuple(ind[:,1][order_1d1])) #bot
    push!(faces, Tuple(ind[end,:][order_1d2]))# right
    push!(faces, Tuple(ind[:,end][order_1d1]))# top
    push!(faces, Tuple(ind[1,:][order_1d2]))# left
    
    return Tuple(faces) 
end

function _compute_volumedof_indices(ip::Bernstein{RefQuadrilateral,order}) where {order}
    orders = (order, order)
    ind = reshape([findfirst(i->i==j, _bernstein_ordering(ip)) for j in 1:prod(orders.+1)], (orders.+1)...)
    volumedofs = ind[2:end-1, 2:end-1]
    return Tuple(volumedofs) 
end

function _compute_vertexdof_indices(::Bernstein{RefHexahedron,order}) where order
    ((1,),(2,),(3,),(4,),(5,),(6,),(7,),(8,))
end

function _compute_facedof_indices(ip::Bernstein{RefHexahedron,order}) where {order}
    faces = Tuple[]
    orders = (order, order, order)

    ind = reshape([findfirst(i->i==j, _bernstein_ordering(ip)) for j in 1:prod(orders.+1)], (orders.+1)...)
    
    #Order for 2d interpolation.
    order_2d = _bernstein_ordering(Bernstein{RefQuadrilateral,order}())
    
    push!(faces, Tuple(ind[:,:,1][order_2d])) # bottom
    push!(faces, Tuple(ind[:,1,:][order_2d]))   # front
    push!(faces, Tuple(ind[end,:,:][order_2d])) # right
    push!(faces, Tuple(ind[:,end,:][order_2d])) # back
    push!(faces, Tuple(ind[1,:,:][order_2d]))   # left
    push!(faces, Tuple(ind[:,:,end][order_2d]))   # top

    return Tuple(faces)
end

function _compute_edgedof_indices(ip::Bernstein{RefHexahedron,order}) where {order}
    edges = Tuple[]
    orders = (order, order, order)

    ind = reshape([findfirst(i->i==j, _bernstein_ordering(ip)) for j in 1:prod(orders.+1)], (orders.+1)...)

    #Order for 1d interpolation.
    order_1d = _bernstein_ordering(Bernstein{1,order}())

    # bottom
    push!(edges, Tuple(ind[:,1,1][order_1d]))
    push!(edges, Tuple(ind[end,:,1][order_1d]))
    push!(edges, Tuple(ind[:,end,1][order_1d]))
    push!(edges, Tuple(ind[1,:,1][order_1d]))

    # top
    push!(edges, Tuple(ind[:,1,end][order_1d]))
    push!(edges, Tuple(ind[end,:,end][order_1d]))
    push!(edges, Tuple(ind[:,end,end][order_1d]))
    push!(edges, Tuple(ind[1,:,end][order_1d]))

    # verticals
    push!(edges, Tuple(ind[1,1,:][order_1d]))
    push!(edges, Tuple(ind[end,1,:][order_1d]))
    push!(edges, Tuple(ind[1,end,:][order_1d]))
    push!(edges, Tuple(ind[end,end,:][order_1d]))

    return Tuple(edges)
end

function _compute_volumedof_indices(ip::Bernstein{RefHexahedron,order}) where {order}
    orders = (order, order, order)
    ind = reshape([findfirst(i->i==j, _bernstein_ordering(ip)) for j in 1:prod(orders.+1)], (orders.+1)...)
    volumedofs = ind[2:end-1, 2:end-1, 2:end-1]
    return Tuple(volumedofs) 
end

function _bernstein_basis_recursive(p::Int, i::Int, xi::T) where T
	if i == 1 && p == 0
		return one(T)
	elseif i < 1 || i > p + 1
		return zero(T)
	else
        return 0.5 * (1 - xi) * _bernstein_basis_recursive(p - 1, i, xi) + 0.5 * (1 + xi) * _bernstein_basis_recursive(p - 1, i - 1, xi)
    end
end

function _bernstein_basis_derivative_recursive(p::Int, i::Int, xi::T) where T
    return p * (_bernstein_basis_recursive(p - 1, i - 1, xi) - _bernstein_basis_recursive(p - 1, i, xi))
end


# # #
# Bernstein generation any order
# # # 
Ferrite.getnbasefunctions(::Bernstein{shape,order}) where {rdim, shape <: AbstractRefShape{rdim}, order}= (order+1)^rdim
    
function Ferrite.reference_coordinates(ip::Bernstein{shape,order}) where {rdim, shape <: AbstractRefShape{rdim}, order}

    T = Float64
    nbasefunks_dim = ntuple(i->order+1, rdim)
    nbasefuncs = prod(nbasefunks_dim)
    
    coords = Vec{rdim,T}[]

    ordering = _bernstein_ordering(ip)
    ranges = [range(T(-1.0), stop=T(1.0), length=nbasefunks_dim[i]) for i in 1:rdim]

    inds = CartesianIndices(nbasefunks_dim)
    for i in 1:nbasefuncs
        ind = inds[ordering[i]]
        x = Vec{rdim,T}(d -> ranges[d][ind[d]])
        push!(coords, x)
    end

    return coords
end

"""
    _bernstein_ordering(::Bernstein)

Return the ordering of the bernstein basis base-functions, from a "CartesianIndices-ordering".
"""
function _bernstein_ordering(::Bernstein{RefLine, order}) where {order}
    ind = reshape(1:prod(order+1), order+1)
    ordering = Int[]

    # Corners
    push!(ordering, ind[1])
    push!(ordering, ind[end])

    # Interior
    append!(ordering, ind[2:end-1])

    return ordering
end

function _bernstein_ordering(::Bernstein{RefQuadrilateral,order}) where {order}
    orders = (order, order)
    ind = reshape(1:prod(orders .+ 1), (orders .+ 1)...)
    ordering = Int[]

    # Corners
    push!(ordering, ind[1,1])
    push!(ordering, ind[end,1])
    push!(ordering, ind[end,end])
    push!(ordering, ind[1,end])

    # edges
    append!(ordering, ind[2:end-1,1])
    append!(ordering, ind[end, 2:end-1])
    append!(ordering, ind[2:end-1,end])
    append!(ordering, ind[1, 2:end-1])

    # inner dofs
    append!(ordering, ind[2:end-1, 2:end-1])

    return ordering
end

function _bernstein_ordering(::Bernstein{RefHexahedron, order}) where {order}
    orders = (order, order, order)
    ind = reshape(1:prod(orders .+ 1), (orders .+ 1)...)
    ordering = Int[]
    
    # Corners, bottom
    push!(ordering, ind[1,  1,  1])
    push!(ordering, ind[end,1,  1])
    push!(ordering, ind[end,end,1])
    push!(ordering, ind[1  ,end,1])

    # Corners, top
    push!(ordering, ind[1,  1  ,end])
    push!(ordering, ind[end,1,  end])
    push!(ordering, ind[end,end,end])
    push!(ordering, ind[1,  end,end])

    # edges, bottom 
    append!(ordering, ind[2:end-1,      1,1])
    append!(ordering, ind[end,    2:end-1,1])
    append!(ordering, ind[2:end-1,end,    1])
    append!(ordering, ind[1,      2:end-1,1])

    # edges, top
    append!(ordering, ind[2:end-1,1,end])
    append!(ordering, ind[end, 2:end-1,end])
    append!(ordering, ind[2:end-1,end,end])
    append!(ordering, ind[1, 2:end-1,end])

    # edges, mid
    append!(ordering, ind[1,     1, 2:end-1])
    append!(ordering, ind[end,   1, 2:end-1])
    append!(ordering, ind[end, end, 2:end-1])
    append!(ordering, ind[1,   end, 2:end-1])

    # Faces (vtk orders left face first, but Ferrite orders bottom first)
    append!(ordering, ind[2:end-1, 2:end-1, 1][:])   # bottom
    append!(ordering, ind[2:end-1, 1, 2:end-1][:])   # front
    append!(ordering, ind[end, 2:end-1, 2:end-1][:]) # right
    append!(ordering, ind[2:end-1, end, 2:end-1][:]) # back
    append!(ordering, ind[1, 2:end-1, 2:end-1][:])   # left
    append!(ordering, ind[2:end-1, 2:end-1, end][:]) # top

    # Inner dofs
    append!(ordering, ind[2:end-1, 2:end-1, 2:end-1][:])

    return ordering
end

#Almost the same orderign as _bernstein_ordering, but some changes for faces and edges
#The ordering is the same as in VTK: https://blog.kitware.com/wp-content/uploads/2020/03/Implementation-of-rational-Be%CC%81zier-cells-into-VTK-Report.pdf.
function _vtk_ordering(c::Type{<:BezierCell{RefLine}})
    _bernstein_ordering(c)
end

function _vtk_ordering(c::Type{<:BezierCell{RefQuadrilateral}}) 
    _bernstein_ordering(c)
end

function _vtk_ordering(::Type{<:BezierCell{RefHexahedron, order}}) where {order}
    orders = (order, order, order)

    ind = reshape(1:prod(orders .+ 1), (orders .+ 1)...)

    # Corners, bottom
    ordering = Int[]
    push!(ordering, ind[1,1,1])
    push!(ordering, ind[end,1,1])
    push!(ordering, ind[end,end,1])
    push!(ordering, ind[1,end,1])
    # Corners, top
    push!(ordering, ind[1,1,end])
    push!(ordering, ind[end,1,end])
    push!(ordering, ind[end,end,end])
    push!(ordering, ind[1,end,end])

    # edges, bottom
    append!(ordering, ind[2:end-1,1,1])
    append!(ordering, ind[end, 2:end-1,1])
    append!(ordering, ind[2:end-1,end,1])
    append!(ordering, ind[1, 2:end-1,1])
    # edges, top
    append!(ordering, ind[2:end-1,1,end])
    append!(ordering, ind[end, 2:end-1,end])
    append!(ordering, ind[2:end-1,end,end])
    append!(ordering, ind[1, 2:end-1,end])
    # edges, mid
    append!(ordering, ind[1,1, 2:end-1])
    append!(ordering, ind[end,1, 2:end-1])
    append!(ordering, ind[1, end, 2:end-1])
    append!(ordering, ind[end, end, 2:end-1])

    # Faces (vtk orders left face first, but Ferrite orders bottom first)
    append!(ordering, ind[1, 2:end-1, 2:end-1][:])   # left
    append!(ordering, ind[end, 2:end-1, 2:end-1][:]) # right
    append!(ordering, ind[2:end-1, 1, 2:end-1][:])   # front
    append!(ordering, ind[2:end-1, end, 2:end-1][:]) # back
    append!(ordering, ind[2:end-1, 2:end-1, 1][:])   # bottom
    append!(ordering, ind[2:end-1, 2:end-1, end][:]) # top

    # Inner dofs
    append!(ordering, ind[2:end-1, 2:end-1, 2:end-1][:])

    return ordering
end
