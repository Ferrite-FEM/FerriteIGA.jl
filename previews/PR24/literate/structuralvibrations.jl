# # [Structural vibrations](@id structural_vibrations)

# In this example, we will compare traditional Finite Element Method (FEM) with Isogeometric Analysis (IGA) for analysing structural vibration problems. We aim to replicate the results from __Cottrell, J. et al, Isogeometric analysis of structural vibrations, In Computer Methods in Applied Mechanics and Engineering [https://doi.org/10.1016/j.cma.2005.09.027](https://doi.org/10.1016/j.cma.2005.09.027)__ 

# Consider the structural vibrations of an elastic rod of unit length, whose natural frequencies and modes are goverened by:

# ```math
# u_{,xx} + \omega_n^2 u = 0, \quad x \in (0,1)
# ```
# subjected to boundary conditions:

# ```math
# u(0) = u(1) = 0
# ```
# where $u$ is the displacement field, and $\omega_n$ is n:th natural frequency. 

# After discretisation, we formulate the generalised eigenvalue problem, which allows us to solve for the natural frequencies and modes $\boldsymbol{\phi}_n$:
# ```math
# (\boldsymbol{K} - \omega^2_n \boldsymbol{M}) \boldsymbol{\phi}_n = 0.
# ```
# Here, $\boldsymbol{K}$ and $\boldsymbol{M}$ are the standard stiffness and mass matrices.

# ## Main code

# Per usual, we first load the relevant packages
using Ferrite, FerriteIGA, LinearAlgebra, Plots

# Next we define the element routine used to compute the element stiffness and mass matrices
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

# We also create a function for computing the natural frequencies $\omega_n$.
# The input for the function is the grid (either a FEM or IGA mesh), and the corresponing cellvalues and interpolations (either Lagrange or IGAInterpolaiton).
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

# We can now solve for the natural frequencies for both the FEM and IGA case. We will use quadratic shape functions and a total of 999 DOFs.

# ## Finite element solution

const N = 999    # Number of desired dofs
const order = 2  # Order of polynomial
const L = 1.0    # Length of beam
const nel_fem = (N-1) ÷ order; # Number of elements
@assert(isodd(N)) #src
grid = generate_grid(Line, (nel_fem,), Vec(0.0), Vec(L))

ip = Lagrange{RefLine,order}()
qr = QuadratureRule{RefLine}(6)
cellvalues = CellValues(qr, ip)

ω_fem = compute_eigenvalues(grid, cellvalues, ip);

# ## Isogeometric solution

const nel_iga = N-order # Number of elements
grid = generate_grid(BezierCell{RefLine,order}, (nel_iga,), Vec(0.0), Vec(L))

ip = IGAInterpolation{RefLine,order}()
qr = QuadratureRule{RefLine}(6)
cellvalues = BezierCellValues(qr, ip)

ω_iga = compute_eigenvalues(grid, cellvalues, ip);

# ## Results
# Now, we can plot and compare the normalised solutions. From the results, we observe that the accuracy of the FEM solution diminishes significantly for $n/N>0.5$. This highlights the advantageous properties of IGA for structural dynamics problems.

analytical_ω_f(n) = π*n
ω_analytical = analytical_ω_f.(1:N)

normalised_ω_fem  = ω_fem ./ ω_analytical
normalised_ω_iga  = ω_iga ./ ω_analytical

n_range = range(0,1,N)
fig = plot(; title="Eigenvalues", ylabel = "Normalised eigenvalue", xlabel="Normalised eigen number")
plot!(fig, n_range, normalised_ω_fem, label="Quadratic FEM")
plot!(fig, n_range, normalised_ω_iga, label="Quadratic IGA")

using Test                      #src
@test ω_iga[10] ≈ 10π atol=1e-4 #src