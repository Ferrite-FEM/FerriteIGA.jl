function Ferrite.default_geometric_interpolation(::IGAInterpolation{shape, order}) where {order, dim, shape <: AbstractRefShape{dim}}
    return VectorizedInterpolation{dim}(IGAInterpolation{shape, order}())
end
function Ferrite.default_geometric_interpolation(::VectorizedInterpolation{vdim,shape,order,IGAInterpolation{shape, order}}) where {vdim,order, dim, shape <: AbstractRefShape{dim}}
    return VectorizedInterpolation{dim}(IGAInterpolation{shape, order}())
end
function Ferrite.default_geometric_interpolation(::Bernstein{shape, order}) where {order, dim, shape <: AbstractRefShape{dim}}
    return VectorizedInterpolation{dim}(Bernstein{shape, order}())
end

"""
    BezierCellValues(qr, ip)

Shape values on an IGA cell. Call `reinit!(cv, coords)` first. Then `shape_value` and
`shape_gradient` return NURBS values.
"""
mutable struct BezierCellValues{FV, GM, QR, T, detJ} <: Ferrite.AbstractCellValues
    const bezier_values::FV # FunctionValues
    const nurbs_values::FV  # FunctionValues
    const geo_mapping::GM   # GeometryMapping
    const qr::QR            # QuadratureRule
    const detJdV::detJ      # AbstractVector{<:Number} or Nothing
    current_beo::BezierExtractionOperator{T}
    const current_w::Vector{T}
end

"""
    BezierFacetValues(fqr, ip)

Shape values on an IGA face. Call `reinit!(fv, coords, faceid)` first.
"""
mutable struct BezierFacetValues{FV, GM, FQR, dim, T, V_FV<:AbstractVector{FV}, V_GM<:AbstractVector{GM}} <: Ferrite.AbstractFacetValues
    const bezier_values::V_FV # FunctionValues
    const nurbs_values::V_FV  # FunctionValues
    const geo_mapping::V_GM   # GeometryMapping
    const fqr::FQR            # QuadratureRule
    const detJdV::Vector{T}
    const normals::Vector{Vec{dim,T}}
    current_beo::BezierExtractionOperator{T}
    const current_w::Vector{T}
    current_facet::Int
end

Ferrite.shape_value_type(cv::BezierCellValues) = Ferrite.shape_value_type(cv.bezier_values)
Ferrite.shape_value_type(cv::BezierFacetValues) = Ferrite.shape_value_type(cv.bezier_values[Ferrite.getcurrentfacet(cv)])

Ferrite.shape_gradient_type(cv::BezierFacetValues) = Ferrite.shape_gradient_type(cv.bezier_values[Ferrite.getcurrentfacet(cv)])
Ferrite.shape_gradient_type(cv::BezierCellValues) = Ferrite.shape_gradient_type(cv.bezier_values)

BezierCellAndFacetValues{T,CV} = Union{BezierCellValues{T,CV}, BezierFacetValues{T,CV}}

Ferrite.nfacets(fv::BezierFacetValues) = length(fv.geo_mapping)
Ferrite.getnormal(fv::BezierFacetValues, iqp::Int) = fv.normals[iqp]
Ferrite.function_interpolation(cv::BezierCellValues) = Ferrite.function_interpolation(cv.bezier_values)
Ferrite.geometric_interpolation(cv::BezierCellValues) = Ferrite.geometric_interpolation(cv.geo_mapping)
Ferrite.function_interpolation(cv::BezierFacetValues) = Ferrite.function_interpolation(cv.bezier_values[1])
Ferrite.geometric_interpolation(cv::BezierFacetValues) = Ferrite.geometric_interpolation(cv.geo_mapping[1])

function BezierCellValues(::Type{T}, qr::QuadratureRule, ip_fun::Interpolation, ip_geo::VectorizedInterpolation, ::ValuesUpdateFlags{FunDiffOrder, GeoDiffOrder, DetJdV}) where {T, FunDiffOrder, GeoDiffOrder, DetJdV}

    geo_mapping = GeometryMapping{GeoDiffOrder}(T, ip_geo.ip, qr)
    fun_values = FunctionValues{FunDiffOrder}(T, ip_fun, qr, ip_geo)
    detJdV = DetJdV ? fill(T(NaN), length(Ferrite.getweights(qr))) : nothing

    undef_beo = BezierExtractionOperator(Vector{SparseArrays.SparseVector{T,Int}}(undef,0))
    undef_w   = NaN .* zeros(Float64, Ferrite.getngeobasefunctions(geo_mapping))

    return BezierCellValues(
        fun_values, 
        deepcopy(fun_values), 
        geo_mapping, qr, detJdV, undef_beo, undef_w)
end

function BezierCellValues(qr::QuadratureRule, ip::Interpolation, args...; kwargs...) 
    return BezierCellValues(Float64, qr, ip, args...; kwargs...)
end
function BezierCellValues(::Type{T}, qr, ip::Interpolation, ip_geo::ScalarInterpolation; kwargs...) where T
    return BezierCellValues(T, qr, ip, VectorizedInterpolation(ip_geo); kwargs...)
end
function BezierCellValues(::Type{T}, qr, ip, ip_geo::VectorizedInterpolation = Ferrite.default_geometric_interpolation(ip); kwargs...) where T
    return BezierCellValues(T, qr, ip, ip_geo, Ferrite.ValuesUpdateFlags(ip; kwargs...))
end

function BezierFacetValues(::Type{T}, fqr::FacetQuadratureRule, ip_fun::Interpolation, ip_geo::VectorizedInterpolation{sdim}, ::ValuesUpdateFlags{FunDiffOrder, GeoDiffOrder, DetJdV}) where {T, sdim, FunDiffOrder, GeoDiffOrder, DetJdV}
    @assert DetJdV
    geo_mapping = [GeometryMapping{GeoDiffOrder}(T, ip_geo.ip, qr) for qr in fqr.facet_rules]
    fun_values = [FunctionValues{FunDiffOrder}(T, ip_fun, qr, ip_geo) for qr in fqr.facet_rules]
    max_nquadpoints = maximum(qr->length(Ferrite.getweights(qr)), fqr.facet_rules)
    detJdV  = fill(T(NaN), max_nquadpoints)
    normals = fill(zero(Vec{sdim, T}) * T(NaN), max_nquadpoints)
    undef_beo = BezierExtractionOperator(Vector{SparseArrays.SparseVector{T,Int}}(undef,0))
    undef_w   = NaN .* zeros(Float64, Ferrite.getngeobasefunctions(first(geo_mapping)))
    return BezierFacetValues(
        fun_values, 
        deepcopy(fun_values), 
        geo_mapping, fqr, detJdV, normals, undef_beo, undef_w, -1)
end

function BezierFacetValues(qr::FacetQuadratureRule, ip::Interpolation, args...; kwargs...) 
    return BezierFacetValues(Float64, qr, ip, args...; kwargs...)
end
function BezierFacetValues(::Type{T}, qr, ip::Interpolation, ip_geo::ScalarInterpolation; kwargs...) where T
    return BezierFacetValues(T, qr, ip, VectorizedInterpolation(ip_geo); kwargs...)
end
function BezierFacetValues(::Type{T}, qr, ip, ip_geo::VectorizedInterpolation = Ferrite.default_geometric_interpolation(ip); kwargs...) where T
    return BezierFacetValues(T, qr, ip, ip_geo, Ferrite.ValuesUpdateFlags(ip; kwargs...))
end


#Intercept construction of CellValues called with IGAInterpolation
function Ferrite.CellValues(
    ::Type{T}, 
    qr::QuadratureRule, 
    ip_fun::Union{IGAInterpolation, VectorizedInterpolation{<:Any, <:Any, <:Any, <: IGAInterpolation}}, 
    ip_geo::VectorizedInterpolation; 
    update_gradients::Bool = true, update_hessians::Bool = false, update_detJdV::Bool = true) where T 

    cv = BezierCellValues(T, qr, ip_fun, ip_geo; update_gradients, update_hessians, update_detJdV)

    return cv
end

#Intercept construction of FacetValues called with IGAInterpolation
function Ferrite.FacetValues(
    ::Type{T}, 
    qr::FacetQuadratureRule, 
    ip_fun::Union{IGAInterpolation, VectorizedInterpolation{<:Any, <:Any, <:Any, <: IGAInterpolation}}, 
    ip_geo::VectorizedInterpolation; 
    update_gradients::Bool = true,
    update_hessians::Bool  = false) where T 

    cv = BezierFacetValues(T, qr, ip_fun, ip_geo; update_gradients, update_hessians)

    return cv
end

Ferrite.getnbasefunctions(fv::BezierFacetValues)            = getnbasefunctions(fv.nurbs_values[Ferrite.getcurrentfacet(fv)])
Ferrite.getnbasefunctions(cv::BezierCellValues)            = getnbasefunctions(cv.nurbs_values)
Ferrite.getngeobasefunctions(fv::BezierFacetValues)         = Ferrite.getngeobasefunctions(fv.geo_mapping[Ferrite.getcurrentfacet(fv)])
Ferrite.getngeobasefunctions(cv::BezierCellValues)         = Ferrite.getngeobasefunctions(cv.geo_mapping)
Ferrite.getnquadpoints(bcv::BezierCellValues)                      = Ferrite.getnquadpoints(bcv.qr)
Ferrite.getnquadpoints(bcv::BezierFacetValues)                      = Ferrite.getnquadpoints(bcv.fqr, Ferrite.getcurrentfacet(bcv))
Ferrite.getdetJdV(bv::BezierCellAndFacetValues, q_point::Int)       = bv.detJdV[q_point]
Ferrite.shape_value(bcv::BezierCellValues, qp::Int, i::Int) = Ferrite.shape_value(bcv.nurbs_values, qp, i)
Ferrite.shape_gradient(bcv::BezierCellValues, q_point::Int, i::Int) = Ferrite.shape_gradient(bcv.nurbs_values, q_point, i)
Ferrite.geometric_value(cv::BezierCellValues, q_point::Int, i::Int) = Ferrite.geometric_value(cv.geo_mapping, q_point, i)

Ferrite.shape_value(fv::BezierFacetValues, qp::Int, i::Int)          = shape_value(fv.nurbs_values[Ferrite.getcurrentfacet(fv)], qp, i)
Ferrite.shape_gradient(fv::BezierFacetValues, q_point::Int, i::Int)  = shape_gradient(fv.nurbs_values[Ferrite.getcurrentfacet(fv)], q_point, i)
Ferrite.geometric_value(fv::BezierFacetValues, q_point::Int, i::Int) = Ferrite.geometric_value(fv.geo_mapping[Ferrite.getcurrentfacet(fv)], q_point, i)

Ferrite.get_fun_values(fv::BezierFacetValues) = @inbounds fv.nurbs_values[Ferrite.getcurrentfacet(fv)]
Ferrite.get_fun_values(cv::BezierCellValues) = @inbounds cv.nurbs_values

Ferrite.shape_hessian(fv::BezierFacetValues, q_point::Int, i::Int) = shape_hessian(fv.nurbs_values[Ferrite.getcurrentfacet(fv)], q_point, i)
Ferrite.shape_hessian(cv::BezierCellValues, q_point::Int, i::Int) = shape_hessian(cv.nurbs_values, q_point, i)

Ferrite.getcurrentfacet(fv::BezierFacetValues) = fv.current_facet
function Ferrite.set_current_facet!(fv::BezierFacetValues, face_nr::Int)
    checkbounds(Bool, 1:Ferrite.nfacets(fv), face_nr) || throw(ArgumentError("Face index out of range."))
    fv.current_facet = face_nr
end

function set_bezier_operator!(bcv::BezierCellAndFacetValues, beo::BezierExtractionOperator{T}) where T 
    bcv.current_beo=beo
end

function set_bezier_operator!(bcv::BezierCellAndFacetValues, beo::BezierExtractionOperator{T}, w::Vector{T}) where T 
    bcv.current_w   .= w
    bcv.current_beo=beo
end

#This function can be called when we know that the weights are all equal to one.
function Ferrite.spatial_coordinate(cv::BezierCellAndFacetValues, iqp::Int, xb::Vector{<:Vec{dim,T}}) where {dim,T}
    wb = ones(T, length(xb))
    x = spatial_coordinate(cv, iqp, (xb, wb))
    return x
end

function Ferrite.spatial_coordinate(cv::BezierCellAndFacetValues, iqp::Int, bcoords::BezierCoords)
    x = spatial_coordinate(cv, iqp, (bcoords.xb, bcoords.wb))
    return x
end

function Ferrite.spatial_coordinate(cv::BezierCellAndFacetValues, iqp::Int, (xb, wb)::CoordsAndWeight{sdim,T}) where {sdim,T}
    nbasefunks = Ferrite.getngeobasefunctions(cv)
    @boundscheck Ferrite.checkquadpoint(cv, iqp)
    W = 0.0
    x = zero(Vec{sdim,T})
    for i in 1:nbasefunks
        N = Ferrite.geometric_value(cv, iqp, i)
        x += N * wb[i] * xb[i]
        W += N * wb[i]
    end
    x /= W
    return x
end

Ferrite.reinit!(cv::BezierCellValues, bc::BezierCoords) = reinit!(cv, nothing, bc)

function Ferrite.reinit!(cv::BezierCellValues, ::Union{Ferrite.AbstractCell, Nothing}, bc::BezierCoords)
    if bc.beo === nothing
        #TODO: If bc.beo === nothing, perhaps we should reinit the normal way like CellValues does?
        error("The bezierextraction matrix is === nothing.")
    end
    set_bezier_operator!(cv, bc.beo, bc.w)
    return reinit!(cv, (bc.xb, bc.wb))
end

function Ferrite.reinit!(cv::BezierCellValues, (x,w)::CoordsAndWeight)
    n_geom_basefuncs = Ferrite.getngeobasefunctions(cv.geo_mapping)
    @assert isa(Ferrite.mapping_type(cv.bezier_values), Ferrite.IdentityMapping)
    @assert checkbounds(Bool, x, 1:n_geom_basefuncs)
    @assert checkbounds(Bool, w, 1:n_geom_basefuncs)

    _bezier_transform(cv.nurbs_values, cv.bezier_values, cv.current_beo, cv.current_w)

    for (q_point, gauss_w) in enumerate(Ferrite.getweights(cv.qr))
        mapping = Ferrite.calculate_mapping(cv.geo_mapping, q_point, x, w)
        
        if cv.detJdV !== nothing
            detJ = Ferrite.calculate_detJ(Ferrite.getjacobian(mapping))
            detJ > 0.0 || Ferrite.throw_detJ_not_pos(detJ)
            cv.detJdV[q_point] = detJ * gauss_w
        end
        _compute_intermidiate!(cv.nurbs_values, cv.geo_mapping, q_point, w)
        if cv.detJdV !== nothing
            Ferrite.apply_mapping!(cv.nurbs_values, q_point, mapping)
        end
    end
    return nothing
end

Ferrite.reinit!(fv::BezierFacetValues, bc::BezierCoords, facet_nr::Int) = reinit!(fv, nothing, bc, facet_nr)

function Ferrite.reinit!(cv::BezierFacetValues, _, bc::BezierCoords, face_nr::Int)
    set_bezier_operator!(cv, bc.beo, bc.w)
    return reinit!(cv, (bc.xb, bc.wb), face_nr)
end

function Ferrite.reinit!(fv::BezierFacetValues, (x,w)::CoordsAndWeight, face_nr::Int)
    Ferrite.set_current_facet!(fv, face_nr) 
    geo_mapping   = fv.geo_mapping[face_nr]
    bezier_values = fv.bezier_values[face_nr]
    nurbs_values  = fv.nurbs_values[face_nr]

    n_geom_basefuncs = Ferrite.getngeobasefunctions(geo_mapping)
    @assert isa(Ferrite.mapping_type(bezier_values), Ferrite.IdentityMapping)
    @assert checkbounds(Bool, x, 1:n_geom_basefuncs)
    @assert checkbounds(Bool, w, 1:n_geom_basefuncs)

    _bezier_transform(nurbs_values, bezier_values, fv.current_beo, fv.current_w)

    for (q_point, gauss_w) in enumerate(Ferrite.getweights(fv.fqr, face_nr))
        mapping = Ferrite.calculate_mapping(geo_mapping, q_point, x, w)
       
        J = Ferrite.getjacobian(mapping)
        weight_norm = Ferrite.weighted_normal(J, Ferrite.getrefshape(geo_mapping.ip), face_nr)
        detJ = norm(weight_norm)
        detJ > 0.0 || Ferrite.throw_detJ_not_pos(detJ)
        @inbounds fv.detJdV[q_point] = detJ * gauss_w
        @inbounds fv.normals[q_point] = weight_norm / detJ

        _compute_intermidiate!(nurbs_values, geo_mapping, q_point, w)
        Ferrite.apply_mapping!(nurbs_values, q_point, mapping)
    end
    return nothing
end


function _bezier_transform(nurbs::FunctionValues{DIFFORDER}, bezier::FunctionValues{DIFFORDER}, Cbe::BezierExtractionOperator{T}, w::Optional{Vector{T}}) where {T,DIFFORDER}
    @assert DIFFORDER < 3
    vdim = length(first(nurbs.Nξ))
   
    transform_basis!(Cbe, nurbs.Nξ, bezier.Nξ, w, stride=vdim)
    if DIFFORDER > 0
        transform_basis!(Cbe, nurbs.dNdξ, bezier.dNdξ, w, stride=vdim)
    end
    if DIFFORDER > 1
        transform_basis!(Cbe, nurbs.d2Ndξ2, bezier.d2Ndξ2, w, stride=vdim)
    end
end

function _compute_intermidiate!(nurbs_values::FunctionValues{0}, geom_values::GeometryMapping, q_point::Int, w::Vector{T}) where {T}
    W = zero(T)
    for j in 1:Ferrite.getngeobasefunctions(geom_values)
        W += w[j]*geom_values.M[j, q_point]
    end
    for j in 1:getnbasefunctions(nurbs_values)
        nurbs_values.Nξ[j,q_point] = nurbs_values.Nξ[j, q_point]/W
    end
end

function _compute_intermidiate!(nurbs_values::FunctionValues{DIFFORDER}, geom_values::GeometryMapping, q_point::Int, w::Vector{T}) where {T,DIFFORDER}
    dim = Ferrite.sdim_from_gradtype(Ferrite.shape_gradient_type(nurbs_values))
    @assert DIFFORDER < 3 "Diff order > 2 not supported"

    W = zero(T)
    dWdξ = zero(Vec{dim,T})
    d2Wdξ2 = zero(Tensor{2,dim,T})
    for j in 1:Ferrite.getngeobasefunctions(geom_values)
        W      += w[j]*geom_values.M[j, q_point]
        if DIFFORDER > 0
            dWdξ   += w[j]*geom_values.dMdξ[j, q_point]
        end
        if DIFFORDER > 1
            d2Wdξ2 += w[j]*geom_values.d2Mdξ2[j, q_point]
        end
    end

    #uses tensor products for vector valued spline functions or scalar otherwise
    is_vector_valued = (first(nurbs_values.Nξ) isa Vec)
    for j in 1:getnbasefunctions(nurbs_values)

        if DIFFORDER > 1
            if is_vector_valued
                _B      = nurbs_values.Nξ[j, q_point]
                _dBdξ   = nurbs_values.dNdξ[j, q_point]
                _d²Bdξ² = nurbs_values.d2Ndξ2[j, q_point]
                tmp = _dBdξ⊗dWdξ
                tmp = permutedims(tmp, (1,3,2))
                tmp = Tensor{3,dim}(tmp)

                Fij = _dBdξ*W - _B⊗dWdξ
                S = W^2
                Fij_k = (_d²Bdξ²*W + _dBdξ⊗dWdξ) - (tmp + _B⊗d2Wdξ2)
                S_k = 2*W*dWdξ
                    
                nurbs_values.d2Ndξ2[j, q_point] = (Fij_k*S - Fij⊗S_k)/S^2
            else
                _B      = nurbs_values.Nξ[j, q_point]
                _dBdξ   = nurbs_values.dNdξ[j, q_point]
                _d²Bdξ² = nurbs_values.d2Ndξ2[j, q_point]

                S = W^2
                Fi = _dBdξ*W - _B⊗dWdξ
                Fi_j = (_d²Bdξ²*W + _dBdξ⊗dWdξ) - (dWdξ⊗_dBdξ + _B⊗d2Wdξ2)
                S_j = 2*W*dWdξ
                nurbs_values.d2Ndξ2[j, q_point] = (Fi_j*S - Fi⊗S_j)/S^2
            end
        end
        
        if DIFFORDER > 0
            if is_vector_valued
                nurbs_values.dNdξ[j, q_point] = (nurbs_values.dNdξ[j, q_point] * W - (nurbs_values.Nξ[j, q_point] ⊗ dWdξ)) / W^2
            else
                nurbs_values.dNdξ[j, q_point] = (nurbs_values.dNdξ[j, q_point] * W - nurbs_values.Nξ[j, q_point] * dWdξ) / W^2
            end
        end
        nurbs_values.Nξ[j,q_point] = nurbs_values.Nξ[j, q_point]/W
    end
end

function Ferrite.calculate_mapping(geo_mapping::Ferrite.GeometryMapping{0}, q_point, x::Vector{Vec{sdim,T}}, w::Vector{T}) where {sdim,T}
    return Ferrite.MappingValues(nothing, nothing)
end

@inline _getrdim(geomapping::Ferrite.GeometryMapping) = length(first(geomapping.dMdξ))
function Ferrite.calculate_mapping(geo_mapping::Ferrite.GeometryMapping{1}, q_point, x::Vector{Vec{sdim,T}}, w::Vector{T}) where {sdim,T}
    rdim = _getrdim(geo_mapping)

    W = zero(T)
    dWdξ = zero(Vec{rdim,T})
    for j in 1:Ferrite.getngeobasefunctions(geo_mapping)
        W      += w[j]*geo_mapping.M[j, q_point]
        dWdξ   += w[j]*geo_mapping.dMdξ[j, q_point]
    end
    
    fecv_J = zero(Ferrite.otimes_returntype(eltype(x), eltype(geo_mapping.dMdξ)))
    for j in 1:Ferrite.getngeobasefunctions(geo_mapping)
        dRdξ = (geo_mapping.dMdξ[j, q_point]*W - geo_mapping.M[j, q_point]*dWdξ)/W^2
        fecv_J += x[j] ⊗ (w[j]*dRdξ)
    end
    return Ferrite.MappingValues(fecv_J, nothing)
end

function Ferrite.calculate_mapping(geo_mapping::Ferrite.GeometryMapping{2}, q_point, x::Vector{Vec{sdim,T}}, w::Vector{T}) where {sdim,T}
    dim = rdim = _getrdim(geo_mapping)
    @assert rdim == sdim
    
    W = zero(T)
    dWdξ = zero(Vec{dim,T})
    d²Wdξ² = zero(Tensor{2,dim,T})
    for j in 1:Ferrite.getngeobasefunctions(geo_mapping)
        W      += w[j]*geo_mapping.M[     j, q_point]
        dWdξ   += w[j]*geo_mapping.dMdξ[  j, q_point]
        d²Wdξ² += w[j]*geo_mapping.d2Mdξ2[j, q_point]
    end
    
    J = zero(Tensor{2,dim,T})
    H = zero(Tensor{3,dim,T})
    for j in 1:Ferrite.getngeobasefunctions(geo_mapping)
        dRdξ = (geo_mapping.dMdξ[j, q_point]*W - geo_mapping.M[j, q_point]*dWdξ)/W^2
        J += x[j] ⊗ (w[j]*dRdξ)
        Fi_j = (geo_mapping.d2Mdξ2[j, q_point]*W +geo_mapping.dMdξ[j, q_point]⊗dWdξ) - (dWdξ⊗geo_mapping.dMdξ[j, q_point] + geo_mapping.M[j, q_point]*d²Wdξ²)
        S_j = 2*W*dWdξ
        S = W^2
        Fi = geo_mapping.dMdξ[j, q_point]*W - geo_mapping.M[j, q_point]*dWdξ
        d²Rdξ² = (Fi_j*S - Fi⊗S_j)/S^2
        H += x[j] ⊗ (w[j]*d²Rdξ²)
    end
    return Ferrite.MappingValues(J, H)
end


function Base.show(io::IO, m::MIME"text/plain", fv::BezierFacetValues)
    println(io, "BezierFacetValues with")
    nqp = getnquadpoints.(fv.fqr.facet_rules)
    fip = Ferrite.function_interpolation(fv)
    gip = Ferrite.geometric_interpolation(fv)
    if all(n==first(nqp) for n in nqp)
        println(io, "- Quadrature rule with ", first(nqp), " points per face")
    else
        println(io, "- Quadrature rule with ", tuple(nqp...), " points on each face")
    end
    print(io, "- Function interpolation: "); show(io, m, fip)
    println(io)
    print(io, "- Geometric interpolation: "); show(io, m, gip)
end

function Base.show(io::IO, d::MIME"text/plain", cv::BezierCellValues)
    ip_geo = geometric_interpolation(cv)
    ip_fun = Ferrite.function_interpolation(cv)
    rdim = Ferrite.getrefdim(ip_geo)
    vdim = isa(shape_value(cv, 1, 1), Vec) ? length(shape_value(cv, 1, 1)) : 0
    GradT = Ferrite.shape_gradient_type(cv)
    sdim = GradT === nothing ? nothing : Ferrite.sdim_from_gradtype(GradT)
    vstr = vdim==0 ? "scalar" : "vdim=$vdim"
    print(io, "BezierCellValues(", vstr, ", rdim=$rdim, and sdim=$sdim): ")
    print(io, getnquadpoints(cv), " quadrature points")
    print(io, "\n Function interpolation: "); show(io, d, ip_fun)
    print(io, "\nGeometric interpolation: ");
    sdim === nothing ? show(io, d, ip_geo) : show(io, d, ip_geo^sdim)
end
