using Ferrite, FerriteIGA, LinearAlgebra, Plots

function stiffness_and_mass_matrix!(ke, me, cv)
    for qp in 1:getnquadpoints(cv)
        dV = getdetJdV(cv,qp)
        for i in 1:getnbasefunctions(cv)
            for j in 1:getnbasefunctions(cv)
                ke[i,j] += shape_gradient(cv, qp, i) ⋅ shape_gradient(cv, qp, j) * dV
                me[i,j] += shape_value(cv, qp, i) * shape_value(cv, qp, j) * dV
            end
        end
    end
end;

function compute_eigenvalues(grid, cellvalues, ip)

    #Create dofhandler
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    #initlize matrices
    nbf = getnbasefunctions(ip)
    ke = zeros(nbf,nbf)
    me = zeros(nbf,nbf)
    K = allocate_matrix(dh)
    M = allocate_matrix(dh)

    #assemble system
    assembler_K = start_assemble(K)
    assembler_M = start_assemble(M)
    for cellid in 1:getncells(grid)
        fill!(ke, 0.0)
        fill!(me, 0.0)
        coords = getcoordinates(grid, cellid)
        dofs = celldofs(dh, cellid)

        reinit!(cellvalues, coords)
        stiffness_and_mass_matrix!(ke, me, cellvalues)

        assemble!(assembler_K, dofs, ke)
        assemble!(assembler_M, dofs, me)
    end

    #apply BC
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(grid, "left"), x->(0.0,)))
    add!(ch, Dirichlet(:u, getfacetset(grid, "right"), x->(0.0,)))
    close!(ch)
    apply!(K, ch)

    #Solve generlized eigenvalue problem
    λ, ϕ = eigen(Matrix(K), Matrix(M)) #Note, we need to convert to full matrices
    ω = sqrt.(λ)
    return ω
end;

const N = 999    # Number of desired dofs
const order = 2  # Order of polynomial
const L = 1.0    # Length of beam
const nel_fem = (N-1) ÷ order; # Number of elements
grid = generate_grid(Line, (nel_fem,), Vec(0.0), Vec(L))

ip = Lagrange{RefLine,order}()
qr = QuadratureRule{RefLine}(6)
cellvalues = CellValues(qr, ip)

ω_fem = compute_eigenvalues(grid, cellvalues, ip);

const nel_iga = N-order # Number of elements
grid = generate_grid(BezierCell{RefLine,order}, (nel_iga,), Vec(0.0), Vec(L))

ip = IGAInterpolation{RefLine,order}()
qr = QuadratureRule{RefLine}(6)
cellvalues = BezierCellValues(qr, ip)

ω_iga = compute_eigenvalues(grid, cellvalues, ip);

analytical_ω_f(n) = π*n
ω_analytical = analytical_ω_f.(1:N)

normalised_ω_fem  = ω_fem ./ ω_analytical
normalised_ω_iga  = ω_iga ./ ω_analytical

n_range = range(0,1,N)
fig = plot(; title="Eigenvalues", ylabel = "Normalised eigenvalue", xlabel="Normalised eigen number")
plot!(fig, n_range, normalised_ω_fem, label="Quadratic FEM")
plot!(fig, n_range, normalised_ω_iga, label="Quadratic IGA")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
