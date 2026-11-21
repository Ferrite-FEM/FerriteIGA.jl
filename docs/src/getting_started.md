# [Getting started](@id getting-started)

FerriteIGA is designed such that an isogeometric program looks, as far as possible, like a
standard Ferrite program. This page lists the (few) places where an isogeometric program
differs from a Ferrite program. It assumes that you are familiar with Ferrite, for example
through the first Ferrite tutorial,
[Heat equation](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/heat_equation/).

## From Ferrite to FerriteIGA

The table below summarizes the differences between a Ferrite program and a FerriteIGA
program. Everything that is not listed, such as the `DofHandler`, the `ConstraintHandler`,
the assembler, the element routine and the linear solve, is exactly the same.

| Step                 | Ferrite                                   | FerriteIGA |
|:---------------------|:------------------------------------------|:-----------|
| Load the package     | `using Ferrite`                           | `using FerriteIGA` (re-exports Ferrite) |
| Geometry             | `generate_grid(Quadrilateral, nels)`      | `generate_nurbs_patch(:name, nels, order)`, a [`NURBSMesh`](@ref) |
| Grid                 | a `Grid`                                  | `BezierGrid(mesh)`, a [`BezierGrid`](@ref) |
| Interpolation        | `Lagrange{RefQuadrilateral, order}()`     | [`IGAInterpolation{RefQuadrilateral, order}()`](@ref IGAInterpolation) |
| Cell values          | `CellValues(qr, ip)`                      | [`BezierCellValues(qr, ip)`](@ref BezierCellValues) |
| Facet values         | `FacetValues(fqr, ip)`                    | [`BezierFacetValues(fqr, ip)`](@ref BezierFacetValues) |
| Cell coordinates     | `getcoordinates(grid, cellid)` returns a vector of nodes | returns a [`BezierCoords`](@ref) |
| Cell loop            | `for cell in CellIterator(dh)`            | `for cellid in 1:getncells(grid)` or an [`IGACellCache`](@ref) |
| Export               | `VTKGridFile(name, grid)`                 | [`VTKIGAFile(name, grid)`](@ref VTKIGAFile) |

The main conceptual difference is in the first two steps. In standard finite elements the
mesh is an approximation of the geometry, and the nodes lie on it. In isogeometric analysis
the geometry is a spline *patch*, described by knot vectors, control points and weights. The
patch is converted to a [`BezierGrid`](@ref), where the nodes are the control points and the
cells are the knot spans. Each cell also stores a *Bézier extraction operator*, which
[`BezierCellValues`](@ref) uses to evaluate the spline basis functions. All of this happens
behind the scenes, and is explained in the [topic guides](topics/index.md).

Continue with the first tutorial, [Infinite plate with hole](@ref tutorial-plate-with-hole),
which solves a linear elasticity problem on the same geometry.
