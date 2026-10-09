#
# Constraints for spline (IGA) degrees of freedom on a BezierGrid.
#
# Spline dofs are coefficients of basis functions, not nodal values, and the basis functions
# that are non-zero on a boundary span several facets. The constraints below therefore work with
# the control points (nodes) of the grid, and are computed directly when they are added to the
# ConstraintHandler.
#

"""
    _node_component_dofs(dh, fieldname, nodes)

Return a dictionary mapping each node (control point) in `nodes` to the dofs of the field `fieldname`,
one dof per component.
"""
function _node_component_dofs(dh::DofHandler, fieldname::Symbol, nodes)
    node_dofs = Dict{Int,Vector{Int}}()
    grid = Ferrite.get_grid(dh)
    for sdh in dh.subdofhandlers
        fieldname in sdh.field_names || continue
        field_idx = Ferrite.find_field(sdh, fieldname)
        drange = Ferrite.dof_range(sdh, field_idx)
        vdim = Ferrite.n_components(sdh, field_idx)
        for cellid in sdh.cellset
            cell = getcells(grid, cellid)
            any(A -> A in nodes && !haskey(node_dofs, A), cell.nodes) || continue
            dofs = celldofs(dh, cellid)
            for (a, A) in pairs(cell.nodes)
                (A in nodes && !haskey(node_dofs, A)) || continue
                node_dofs[A] = [dofs[drange[(a-1)*vdim + c]] for c in 1:vdim]
            end
        end
    end
    return node_dofs
end

_evaluate_bc_function(f, x) = applicable(f, x) ? f(x) : f(x, 0.0)

"""
    ProjectedIGADirichlet(fieldname::Symbol, facets, f::Function; components = nothing, qr_order = -1)

Dirichlet boundary condition for a spline field on a `BezierGrid`, where the prescribed values are computed by
an L2 projection of `f` onto the spline functions that are non-zero on the `facets`. Compared to `Dirichlet`,
which evaluates `f` at the control points, the projection is exact whenever `f` can be represented by the spline
basis on the boundary, and is generally more accurate.

`f(x)` should return a scalar for scalar fields, or one value per constrained component. The projection is
computed once, when the condition is added to the `ConstraintHandler`, so `f` can not depend on time
(a function `f(x, t)` is evaluated at `t = 0`). `components` selects which components of a vector field to constrain
(all by default), and `qr_order` sets the order of the facet quadrature (by default `2p + 2`).

The projection is solved for all facets at once (a global boundary mass matrix), since the spline functions
span several facets.
"""
struct ProjectedIGADirichlet
    field_name::Symbol
    facets::OrderedSet{FacetIndex}
    f::Function
    components::Union{Nothing,Vector{Int}}
    qr_order::Int
end

function ProjectedIGADirichlet(field_name::Symbol, facets, f::Function; components = nothing, qr_order::Int = -1)
    components = components === nothing ? nothing : collect(Int, components)
    return ProjectedIGADirichlet(field_name, OrderedSet{FacetIndex}(facets), f, components, qr_order)
end

"""
    add!(ch::ConstraintHandler{<:DofHandler{dim, <:BezierGrid}}, bc::ProjectedIGADirichlet)

Compute the L2 projection of the boundary values and prescribe the corresponding dofs.
The values are not changed by `update!`.
"""
function Ferrite.add!(ch::ConstraintHandler{<:DofHandler{<:Any,<:BezierGrid}}, bc::ProjectedIGADirichlet)
    dh = ch.dh
    grid = Ferrite.get_grid(dh)
    nn = getnnodes(grid)

    #Boundary mass matrix and right hand side, assembled per control point (node)
    I = Int[]
    J = Int[]
    V = Float64[]
    F = Matrix{Float64}(undef, 0, 0)
    components = Int[]
    for sdh in dh.subdofhandlers
        bc.field_name in sdh.field_names || continue
        facets = [facet for facet in bc.facets if facet[1] in sdh.cellset]
        isempty(facets) && continue

        field_idx = Ferrite.find_field(sdh, bc.field_name)
        ip = Ferrite.getfieldinterpolation(sdh, field_idx)
        vdim = Ferrite.n_components(sdh, field_idx)
        components = bc.components === nothing ? collect(1:vdim) : bc.components
        all(c -> 1 <= c <= vdim, components) || throw(ArgumentError("components $(components) not within range of field :$(bc.field_name) ($vdim dimension(s))"))
        isempty(F) && (F = zeros(nn, length(components)))

        ip_scalar = _bernstein_ip(Ferrite.get_base_interpolation(ip))
        ip_scalar isa IGAInterpolation || throw(ArgumentError("ProjectedIGADirichlet is only supported for IGA (spline) fields, got $ip."))
        p = Ferrite.getorder(ip_scalar)
        shape = Ferrite.getrefshape(ip_scalar)
        qr_order = bc.qr_order > 0 ? bc.qr_order : 2p + 2
        fv = BezierFacetValues(FacetQuadratureRule{shape}(qr_order), ip_scalar)

        for (cellid, facetid) in facets
            coords = getcoordinates(grid, cellid)
            reinit!(fv, coords, facetid)
            nodes = getcells(grid, cellid).nodes
            n = getnbasefunctions(fv)
            for qp in 1:getnquadpoints(fv)
                dΓ = getdetJdV(fv, qp)
                x = spatial_coordinate(fv, qp, coords)
                g = _evaluate_bc_function(bc.f, x)
                length(g) == length(components) || throw(ArgumentError("f(x) must return $(length(components)) value(s), one for each constrained component."))
                for a in 1:n
                    Ra = shape_value(fv, qp, a)
                    for k in eachindex(components)
                        F[nodes[a], k] += Ra * g[k] * dΓ
                    end
                    for b in 1:n
                        push!(I, nodes[a]); push!(J, nodes[b]); push!(V, Ra * shape_value(fv, qp, b) * dΓ)
                    end
                end
            end
        end
    end
    isempty(F) && error("No overlap between ProjectedIGADirichlet and the fields in the ConstraintHandler's DofHandler")

    #Solve for the control points whose basis functions are non-zero on the boundary
    Mb = sparse(I, J, V, nn, nn)
    d = [Mb[i,i] for i in 1:nn]
    boundary_nodes = findall(>(1e-12 * maximum(d)), d)
    U = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(Mb[boundary_nodes, boundary_nodes])) \ F[boundary_nodes, :]

    node_dofs = _node_component_dofs(dh, bc.field_name, Set(boundary_nodes))
    for (k, A) in pairs(boundary_nodes)
        for (ic, c) in pairs(components)
            Ferrite.add_prescribed_dof!(ch, node_dofs[A][c], U[k, ic])
        end
    end
    return ch
end

"""
    PeriodicIGADirichlet(fieldname::Symbol, mirror_facets, image_facets; components = nothing, tol = 1e-10)

Periodic boundary condition for a spline field on a `BezierGrid`. The control points on the mirror side
(those whose basis functions are non-zero on `mirror_facets`) are tied to the control points on the image side,
`u_mirror = u_image`, with affine constraints.

The two sides must be identical up to a translation: each control point on the mirror side must match exactly one
control point (with the same weight) on the image side, after translating with the difference of the centroids of
the control points on the two sides. Otherwise an error is thrown. The matching uses a simple O(N²) search, with the
relative tolerance `tol`.

Note that only the control points on the boundary are tied, so the solution is periodic with C⁰ continuity.
Periodic conditions in several directions (e.g. a doubly periodic square) can be added one by one; constraints are
chained such that no nested affine constraints are created. A `Dirichlet` condition added afterwards on a tied
dof replaces that tie, so add Dirichlet conditions first.
"""
struct PeriodicIGADirichlet
    field_name::Symbol
    mirror::OrderedSet{FacetIndex}
    image::OrderedSet{FacetIndex}
    components::Union{Nothing,Vector{Int}}
    tol::Float64
end

function PeriodicIGADirichlet(field_name::Symbol, mirror, image; components = nothing, tol::Real = 1e-10)
    components = components === nothing ? nothing : collect(Int, components)
    return PeriodicIGADirichlet(field_name, OrderedSet{FacetIndex}(mirror), OrderedSet{FacetIndex}(image), components, Float64(tol))
end

"""
    add!(ch::ConstraintHandler{<:DofHandler{dim, <:BezierGrid}}, bc::PeriodicIGADirichlet)

Match the control points on the two sides and add the affine constraints `u_mirror = u_image`.
"""
function Ferrite.add!(ch::ConstraintHandler{<:DofHandler{<:Any,<:BezierGrid}}, bc::PeriodicIGADirichlet)
    dh = ch.dh
    grid = Ferrite.get_grid(dh)
    mirror_nodes = collect(_facets_to_nodes(grid, bc.mirror))
    image_nodes = collect(_facets_to_nodes(grid, bc.image))
    n = length(mirror_nodes)
    if n != length(image_nodes)
        throw(ArgumentError("The sides of the PeriodicIGADirichlet are not identical: the mirror side has $n control points, and the image side has $(length(image_nodes))."))
    end

    xm = [grid.nodes[A].x for A in mirror_nodes]
    xi = [grid.nodes[A].x for A in image_nodes]
    offset = sum(xi)/n - sum(xm)/n
    diam = maximum(x -> norm(x - first(xm)), [xm; xi])
    tol = bc.tol * max(diam, one(diam))

    #Match each control point on the mirror side with one on the image side
    image_of = zeros(Int, n)
    used = falses(n)
    for (k, x) in pairs(xm)
        candidates = findall(l -> norm(x + offset - xi[l]) <= tol, 1:n)
        if length(candidates) != 1 || used[only(candidates)]
            throw(ArgumentError("The sides of the PeriodicIGADirichlet are not identical: no unique match for the control point at $x on the mirror side."))
        end
        l = only(candidates)
        if abs(grid.weights[mirror_nodes[k]] - grid.weights[image_nodes[l]]) > bc.tol * grid.weights[mirror_nodes[k]]
            throw(ArgumentError("The sides of the PeriodicIGADirichlet are not identical: the weights of the control points at $x and $(xi[l]) differ."))
        end
        used[l] = true
        image_of[k] = l
    end

    node_dofs = _node_component_dofs(dh, bc.field_name, Set([mirror_nodes; image_nodes]))
    isempty(node_dofs) && error("No overlap between PeriodicIGADirichlet and the fields in the ConstraintHandler's DofHandler")
    vdim = length(first(values(node_dofs)))
    components = bc.components === nothing ? collect(1:vdim) : bc.components
    all(c -> 1 <= c <= vdim, components) || throw(ArgumentError("components $(components) not within range of field :$(bc.field_name) ($vdim dimension(s))"))
    for (k, A) in pairs(mirror_nodes)
        B = image_nodes[image_of[k]]
        for c in components
            _add_periodic_pair!(ch, node_dofs[A][c], node_dofs[B][c])
        end
    end
    return ch
end

#Follow constraints of the type `d = 1.0 * d'` (without inhomogeneity) to the representative dof
function _periodic_representative(ch::ConstraintHandler, d::Int)
    while true
        i = findfirst(==(d), ch.prescribed_dofs)
        i === nothing && return d
        coeffs = ch.dofcoefficients[i]
        if coeffs !== nothing && length(coeffs) == 1 && only(coeffs).second == 1 && iszero(ch.inhomogeneities[i])
            d = only(coeffs).first
        else
            return d
        end
    end
end

_is_free_dof(ch::ConstraintHandler, d::Int) = findfirst(==(d), ch.prescribed_dofs) === nothing

#Tie the dofs `s` and `m`, such that no nested affine constraints are created (which Ferrite does not support).
function _add_periodic_pair!(ch::ConstraintHandler, s::Int, m::Int)
    rs = _periodic_representative(ch, s)
    rm = _periodic_representative(ch, m)
    rs == rm && return ch
    if _is_free_dof(ch, rs)
        slave, master = rs, rm
    elseif _is_free_dof(ch, rm)
        slave, master = rm, rs
    else
        return ch #Both dofs are already prescribed (e.g. by Dirichlet conditions)
    end
    #Redirect constraints which depend on the dof that is now constrained
    for coeffs in ch.dofcoefficients
        coeffs === nothing && continue
        for (k, (d, w)) in pairs(coeffs)
            d == slave && (coeffs[k] = master => w)
        end
    end
    Ferrite.add_prescribed_dof!(ch, slave, zero(eltype(ch.inhomogeneities)), [master => one(eltype(ch.inhomogeneities))])
    return ch
end

#Ferrite's ProjectedDirichlet does not support spline (IGA) fields
function Ferrite.add!(::ConstraintHandler{<:DofHandler{<:Any,<:BezierGrid}}, ::Ferrite.ProjectedDirichlet)
    throw(ArgumentError("ProjectedDirichlet is not supported on a BezierGrid. Use ProjectedIGADirichlet instead."))
end
