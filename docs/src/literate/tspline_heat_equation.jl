# # Transient heat equation on a graded T-spline mesh
#
# In this example we solve the transient heat equation on a rectangle, using a T-spline mesh
# that is fine at the top and coarse at the bottom. The top boundary is heated at a few hot spots, which
# gives steep temperature gradients close to the top. Further down, the temperature is smooth, so
# large elements are enough there.
#
# With NURBS, refining the top of the domain would also refine the bottom, since knot lines span the whole patch.
# With T-splines, the vertical lines of the fine top layers can end in *T-junctions*, so the refinement stays local.
# See the [T-splines](@ref) page for background.
#
# The strong form of the problem is
# ```math
#  \frac{\partial u}{\partial t} - \nabla \cdot (k \nabla u) = 0  \quad \textbf{x} \in \Omega,
# ```
# where ``u`` is the temperature and ``k`` the heat conductivity. The domain is ``\Omega = [0, 2] \times [0, 1]``. On the top boundary,
# the temperature is prescribed as
# ```math
#  u(x, 1, t) = \min(t/t_\mathrm{ramp}, 1) \, \frac{1 - \cos(8 \pi x / 2)}{2},
# ```
# i.e. four hot spots which are ramped up in time. The bottom boundary is kept at ``u = 0`` and the left and right
# boundaries are insulated (zero heat flux). The initial temperature is ``u = 0``.
#
# The weak form is: find ``u`` such that
# ```math
#  \int_\Omega \delta u \frac{\partial u}{\partial t} \, \mathrm{d}\Omega + \int_\Omega k \nabla \delta u \cdot \nabla u \, \mathrm{d}\Omega = 0 \quad \forall \delta u.
# ```
# Inserting the spline approximation gives the semi-discrete system ``\boldsymbol{M} \dot{\boldsymbol{u}} + \boldsymbol{K} \boldsymbol{u} = \boldsymbol{0}``, which we
# integrate in time with the implicit Euler method,
# ```math
#  (\boldsymbol{M} + \Delta t \boldsymbol{K}) \, \boldsymbol{u}_{n+1} = \boldsymbol{M} \boldsymbol{u}_n.
# ```

#md # !!! note
#md #     It is expected that the reader is familiar with the [Infinite plate with hole](@ref) example and with
#md #     Ferrite's DofHandler, ConstraintHandler and CellValues.

# Start by loading the necessary packages. `WriteVTK` is used to create a ParaView collection (`.pvd` file) of the time steps.
using Ferrite, FerriteIGA, SparseArrays, LinearAlgebra
using WriteVTK: paraview_collection

# ## Mesh
# We create a graded T-spline mesh with `generate_tspline_patch(:graded_rectangle, ...)`. The rectangle is divided
# into horizontal layers, from bottom to top. Each layer has its own number of elements in the x- and y-direction, and its own height.
# Here, the number of columns is doubled in each layer, and the layers get thinner towards the top.
order = 2
nels_x = [8, 16, 32]           # columns in each layer, bottom to top
nels_y = [2, 4, 8]             # rows in each layer
layer_heights = [0.5, 0.3, 0.2]
mesh = generate_tspline_patch(:graded_rectangle, nels_x, nels_y, order;
    size = (2.0, 1.0), cornerpos = (0.0, 0.0), layer_heights)

# The vertical lines of a layer end in T-junctions at the interface to the coarser layer below. Since all horizontal
# lines span the whole width, all T-junctions are vertical, and the mesh is analysis-suitable. Therefore, the basis
# functions are linearly independent, form a partition of unity, and have the same approximation properties as NURBS.
#is_analysis_suitable(mesh)

# As for NURBS patches, the T-spline is converted to a `BezierGrid` for the analysis. Each cell is a Bézier element
# with a Bézier extraction operator. Note that the elements are bounded by all lines of reduced continuity of the
# basis functions, so there are more Bézier elements than faces in the T-mesh.
grid = BezierGrid(mesh)

# The grid contains the facetsets (and nodesets) `"left"`, `"right"`, `"bottom"` and `"top"` of the rectangle.
# Note that `addfacetset!` can not be used for T-spline grids, since the facets of a T-spline cell can not be
# computed from its control points.

# ## Cell values and degrees of freedom
# The cell values and the dof handler are created exactly as for NURBS. The interpolation is the Bernstein basis
# `IGAInterpolation`, and `reinit!` applies the Bézier extraction to obtain the T-spline basis functions.
ip = IGAInterpolation{RefQuadrilateral, order}()
qr = QuadratureRule{RefQuadrilateral}(order + 1)
cellvalues = BezierCellValues(qr, ip)

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

# ## Boundary conditions
# The Dirichlet conditions are added on the facetsets. For T-splines (and NURBS), the degrees of freedom are not nodal values,
# so the condition is applied to the control points whose basis functions are non-zero on the boundary, with the value of the
# function evaluated at the control points. Since the boundary is open (interpolatory), this is exact for the
# zero temperature, and a good approximation of the smooth temperature at the top.
k = 1.0e-1      # heat conductivity
t_ramp = 0.1
top_temperature(x, t) = min(t / t_ramp, 1.0) * (1 - cos(8π * x[1] / 2)) / 2

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> top_temperature(x, t)))
add!(ch, Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> 0.0))
close!(ch);

# ## Assembly
# We assemble the mass matrix ``\boldsymbol{M}`` and the conductivity matrix ``\boldsymbol{K}``. The assembly loop is the same
# as for a standard finite element code. The only difference is that we get the cell coordinates as `BezierCoords` (containing both the
# control points, the weights and the extraction operator) with `getcoordinates`, which are passed to `reinit!`.
#
# In general, the number of T-spline basis functions can vary between the elements (for meshes which are not analysis-suitable),
# so we ask the cell values for the number of basis functions after `reinit!`, instead of using `ndofs_per_cell(dh)`.
function assemble_element!(Me, Ke, cellvalues, k)
    n = getnbasefunctions(cellvalues)
    for qp in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, qp)
        for i in 1:n
            δu = shape_value(cellvalues, qp, i)
            ∇δu = shape_gradient(cellvalues, qp, i)
            for j in 1:n
                u = shape_value(cellvalues, qp, j)
                ∇u = shape_gradient(cellvalues, qp, j)
                Me[i, j] += δu * u * dΩ
                Ke[i, j] += k * (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return Me, Ke
end

function assemble_global(cellvalues, dh, k)
    grid = dh.grid
    M = allocate_matrix(dh)
    K = allocate_matrix(dh)
    M_assembler = start_assemble(M)
    K_assembler = start_assemble(K)
    for cellid in 1:getncells(grid)
        bezier_coords = getcoordinates(grid, cellid)
        reinit!(cellvalues, bezier_coords)

        n = getnbasefunctions(cellvalues)
        Me = zeros(n, n)
        Ke = zeros(n, n)
        assemble_element!(Me, Ke, cellvalues, k)

        dofs = celldofs(dh, cellid)
        assemble!(M_assembler, dofs, Me)
        assemble!(K_assembler, dofs, Ke)
    end
    return M, K
end

M, K = assemble_global(cellvalues, dh, k);

# ## Time stepping
# With the implicit Euler method, the system matrix ``\boldsymbol{A} = \boldsymbol{M} + \Delta t \boldsymbol{K}`` is constant.
# The Dirichlet conditions are applied to it once. Since the prescribed values change in time, we store the data needed
# to update the right hand side with `get_rhs_data` before `apply!` modifies the matrix.
Δt = 0.01
T = 1.0

A = M + Δt * K
rhsdata = get_rhs_data(ch, A)
apply!(A, ch)
A_factorized = cholesky(Symmetric(A))

# The solution is exported to a VTK file in every time step, and the files are collected in a ParaView collection (`.pvd`),
# which can be opened in ParaView to animate the solution. The `.vtu` files contain the Bézier elements, so the geometry
# and the solution are represented exactly (use a "Nonlinear Subdivision Level" larger than zero in ParaView to view them smoothly).
u = zeros(ndofs(dh))
pvd = paraview_collection("tspline_heat_equation")

VTKIGAFile("tspline_heat_equation-0", grid) do vtk
    write_solution(vtk, dh, u)
    pvd[0.0] = vtk
end

t = 0.0
step = 0
while t < T - Δt/2
    global t += Δt
    global step += 1
    ## Update the prescribed values to the new time, and solve for the new temperature
    update!(ch, t)
    b = M * u
    apply_rhs!(rhsdata, b, ch)
    u .= A_factorized \ b

    VTKIGAFile("tspline_heat_equation-$step", grid) do vtk
        write_solution(vtk, dh, u)
        pvd[t] = vtk
    end
end
close(pvd);

# As a simple check, the temperature should stay between the minimum and maximum of the prescribed temperatures.
# This holds approximately, since the spline basis is non-negative.
extrema(u)

#md # ## [Plain program](@id tspline_heat_equation-plain-program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`tspline_heat_equation.jl`](tspline_heat_equation.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```
