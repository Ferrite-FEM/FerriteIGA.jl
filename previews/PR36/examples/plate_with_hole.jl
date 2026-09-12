using Ferrite, FerriteIGA, LinearAlgebra

function integrate_element!(ke::AbstractMatrix, C::SymmetricTensor{4,2}, cv)
    n_basefuncs = getnbasefunctions(cv)

    δɛ = [zero(SymmetricTensor{2,2,Float64}) for i in 1:n_basefuncs]
    for q_point in 1:getnquadpoints(cv)

        for i in 1:n_basefuncs
            δɛ[i] = symmetric(shape_gradient(cv, q_point, i))
        end

        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs
            for j in 1:n_basefuncs
                ke[i, j] += (δɛ[i] ⊡ C ⊡ δɛ[j]) * dΩ
            end
        end
    end
end;

function integrate_traction_force!(fe::AbstractVector, t::Vec{2}, fv)
    n_basefuncs = getnbasefunctions(fv)

    for q_point in 1:getnquadpoints(fv)
        dA = getdetJdV(fv, q_point)
        for i in 1:n_basefuncs
            δu = shape_value(fv, q_point, i)
            fe[i] += t ⋅ δu * dA
        end
    end
end;

function assemble_problem(dh::DofHandler, grid, cv, fv, stiffmat, traction)

    f = zeros(ndofs(dh))
    K = allocate_matrix(dh)
    assembler = start_assemble(K, f)

    n = getnbasefunctions(cv)
    celldofs = zeros(Int, n)
    fe = zeros(n)     # element force vector
    ke = zeros(n, n)  # element stiffness matrix

    n = Ferrite.nnodes_per_cell(grid)
	w = zeros(Float64, n)
	x = zeros(Vec{2}, n)
	wb = zeros(Float64, n)
	xb = zeros(Vec{2}, n)

    # Assemble internal forces
    for cellid in 1:getncells(grid)
        fill!(fe, 0.0)
        fill!(ke, 0.0)
        celldofs!(celldofs, dh, cellid)

        # In a normal finite elment code, this is the point where we usually get the coordinates of the element `X = getcoordinates(grid, cellid)`. In this case, however,
        # we also require the cell weights, and we need to transform them to the bezier mesh.
        extr = get_extraction_operator(grid, cellid) # Extraction operator
        get_bezier_coordinates!(xb,wb,x,w,grid,cellid) #Nurbs coords
        set_bezier_operator!(cv, extr, w)
        reinit!(cv, (xb,wb)) ## Reinit cellvalues by passsing both bezier coords and weights
        integrate_element!(ke, stiffmat, cv)

        assemble!(assembler, celldofs, ke, fe)
    end

    # Assamble external forces
    for (cellid, faceid) in getfacetset(grid, "left")
        fill!(fe, 0.0)

        celldofs!(celldofs, dh, cellid)

        beziercoords = getcoordinates(grid, cellid)
        reinit!(fv, beziercoords, faceid)

        integrate_traction_force!(fe, traction, fv)
        f[celldofs] += fe
    end


    return K, f
end;

function get_material(; E, ν)
    λ = E*ν / ((1 + ν) * (1 - 2ν))
    μ = E / (2(1 + ν))
    δ(i,j) = i == j ? 1.0 : 0.0
    g(i,j,k,l) = λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))

    return SymmetricTensor{4, 2}(g)
end;

function calculate_stress(dh, cv::BezierCellValues, C::SymmetricTensor{4,2}, u::Vector{Float64})

    celldofs = zeros(Int, ndofs_per_cell(dh))

    #Store the stresses in each qp for all cells
    cellstresses = Vector{SymmetricTensor{2,2,Float64,3}}[]

    for cellid in 1:getncells(dh.grid)

        bc = getcoordinates(dh.grid, cellid)

        reinit!(cv, bc)
        celldofs!(celldofs, dh, cellid)

        ue = u[celldofs]
        qp_stresses = SymmetricTensor{2,2,Float64,3}[]
        for qp in 1:getnquadpoints(cv)
            ε = symmetric(function_gradient(cv, qp, ue))
            σ = C ⊡ ε
            push!(qp_stresses, σ)
        end
        push!(cellstresses, qp_stresses)
    end

    return cellstresses
end;

order = 2 # order of the NURBS
nels = (20, 10) # Number of elements
nurbsmesh = generate_nurbs_patch(:plate_with_hole, nels, order);

grid = BezierGrid(nurbsmesh);

addnodeset!(grid, "right", (x) -> x[1] ≈ -0.0)
addfacetset!(grid, "left", (x) -> x[1] ≈ -4.0)
addfacetset!(grid, "bot", (x) -> x[2] ≈ 0.0)
addfacetset!(grid, "right", (x) -> x[1] ≈ 0.0);

ip_geo = IGAInterpolation{RefQuadrilateral,order}()
ip_u = ip_geo^2
qr_cell = QuadratureRule{RefQuadrilateral}(4)
qr_face = FacetQuadratureRule{RefQuadrilateral}(3);

cv = BezierCellValues(qr_cell, ip_u)
fv = BezierFacetValues(qr_face, ip_u);

dh = DofHandler(grid)
add!(dh, :u, ip_u)
close!(dh);

ch = ConstraintHandler(dh)
dbc1 = Dirichlet(:u, getfacetset(grid, "bot"), (x, t) -> 0.0, 2)
dbc2 = Dirichlet(:u, getfacetset(grid, "right"), (x, t) -> 0.0, 1)
add!(ch, dbc1)
add!(ch, dbc2)
close!(ch)
update!(ch, 0.0);

stiffmat = get_material(E=100, ν=0.3)
traction = Vec((-10.0, 0.0))
K, f = assemble_problem(dh, grid, cv, fv, stiffmat, traction);

apply!(K, f, ch)
u = K \ f;

cellstresses = calculate_stress(dh, cv, stiffmat, u);

VTKIGAFile("plate_with_hole.vtu", grid) do vtk
    write_solution(vtk, dh, u)
    #IGA.write_projections(vtk, projector, σ_nodes, "σ")
end;

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
