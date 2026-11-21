```@setup export
using FerriteIGA
grid = BezierGrid(generate_nurbs_patch(:plate_with_hole, (8, 4), 2))
ip = IGAInterpolation{RefQuadrilateral, 2}()
qr = QuadratureRule{RefQuadrilateral}(3)
cellvalues = BezierCellValues(qr, ip)
dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
u = rand(ndofs(dh)) # some data
```

# [Export](@id export_topicguide)

When the problem is solved, and the solution vector `u` is known, we typically want to
visualize it. As in Ferrite, the simplest way to do this is to write the solution to a VTK
file, which can be viewed in e.g. [ParaView](https://www.paraview.org/). For a
[`BezierGrid`](@ref), Ferrite's `VTKGridFile` is replaced with a [`VTKIGAFile`](@ref), which is
used in exactly the same way:

```@example export
VTKIGAFile("my_solution.vtu", grid) do vtk
    write_solution(vtk, dh, u)
end
nothing # hide
```

Instead of using the `do`-block, it is also possible to do

```@example export
vtk = VTKIGAFile("my_solution.vtu", grid)
write_solution(vtk, dh, u)
# etc.
close(vtk)
nothing # hide
```

The file is written to disk when it is closed. The following functions can be used to add
data to the file:

 - `write_solution(vtk, dh, u)`: the fields of the `DofHandler` `dh`, for the solution
   vector `u`,
 - `write_projection(vtk, proj, vals, name)`: the result of an `L2Projector`,
 - `write_node_data(vtk, nodedata, name)`: data given per node (control point),
 - `write_cell_data(vtk, celldata, name)`: data given per cell.

They have the same meaning as for a `VTKGridFile`, see the
[Ferrite export topic guide](https://ferrite-fem.github.io/Ferrite.jl/stable/topics/export/).

## How the grid is written

VTK supports *Bézier cells* (`VTK_BEZIER_CURVE`, `VTK_BEZIER_QUADRILATERAL` and
`VTK_BEZIER_HEXAHEDRON`), which are rational Bézier cells defined by their Bézier control
points and weights. These are exactly the cells obtained with Bézier extraction (see the
[BezierGrid](@ref beziergrid_topicguide) topic guide), so `VTKIGAFile` writes each cell of the
`BezierGrid` as a Bézier cell with its Bézier control points as the points of the cell. The
Bézier weights are written as the point data array `RationalWeights`, which ParaView uses to
evaluate the rational geometry.

Since the Bézier control points are local to each cell, the VTK grid is *discontinuous*:
every cell has its own copy of its points, and neighboring cells are not connected in the VTK
file. The VTK file therefore has more points than the `BezierGrid` has nodes.

!!! tip "Visualization of curved cells in ParaView"
    ParaView renders Bézier cells by subdividing them. The number of subdivisions is set by
    the *Nonlinear Subdivision Level* property of the representation (in the advanced
    properties). Increase it to see the curved geometry and the variation of the fields
    within the cells.

## Solution fields and projections

The spline solution is not interpolatory: the degrees of freedom are coefficients of the
basis functions, not values of the solution at the control points. `write_solution` and
`write_projection` therefore *evaluate* the fields in each cell, at the points of the
reference cell that correspond to the Bernstein polynomials of the
[`IGAInterpolation`](@ref) (equidistant points, including the vertices), and attach these
values to the points of the cell in the VTK file. ParaView then interpolates the values over
the cell, so the visualized field is an approximation of the spline field which converges as
the mesh is refined.

Quantities computed in the quadrature points, such as stresses or fluxes, can be projected to
the spline space with Ferrite's `L2Projector` and then written with `write_projection`. As an
example, we project the gradient of `u`:

```@example export
gradients = [Vec{2, Float64}[] for _ in 1:getncells(grid)]
for cellid in 1:getncells(grid)
    reinit!(cellvalues, getcoordinates(grid, cellid))
    ue = u[celldofs(dh, cellid)]
    for q_point in 1:getnquadpoints(cellvalues)
        push!(gradients[cellid], function_gradient(cellvalues, q_point, ue))
    end
end

projector = L2Projector(ip, grid)
∇u_projected = project(projector, gradients, qr)

VTKIGAFile("my_solution.vtu", grid) do vtk
    write_solution(vtk, dh, u)
    write_projection(vtk, projector, ∇u_projected, "∇u")
end
nothing # hide
```

## Node and cell data

`write_node_data` writes data given per node of the `BezierGrid`, i.e. per control point.
The value of a node is assigned to all copies of that node in the discontinuous VTK grid.
Note that the points of the VTK grid are the Bézier control points, not the NURBS control
points, so node data is mainly useful for data that is constant or for debugging. Cell data
is written with `write_cell_data`, exactly as in Ferrite:

```@example export
VTKIGAFile("my_solution.vtu", grid) do vtk
    write_cell_data(vtk, collect(1.0:getncells(grid)), "cell number")
end
nothing # hide
```

## Time series

For simulations with multiple time steps, one VTK file is written for each step, and the
files are collected in a ParaView data file (`.pvd`) with `paraview_collection` from
[WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl). The underlying WriteVTK file of a
`VTKIGAFile` is stored in its `vtk` field:

```@example export
using WriteVTK
pvd = paraview_collection("my_results")
for (step, t) in enumerate(range(0, 1, 5))
    # Do calculations to update u
    VTKIGAFile("my_results_$step.vtu", grid) do vtk
        write_solution(vtk, dh, u)
        pvd[t] = vtk.vtk
    end
end
vtk_save(pvd)
nothing # hide
```
