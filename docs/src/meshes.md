# Meshes

Functions that generate meshes. A generated patch is a [`NURBSMesh`](@ref) (see the [NURBSMesh](@ref) page for refinement) or a [`TSplineMesh`](@ref), and is converted to a grid for the analysis with `BezierGrid(mesh)`.

```@docs
generate_nurbs_patch
```

T-spline meshes are generated with [`generate_tspline_patch`](@ref), which is documented on the [T-splines](@ref) page.

## Patches of `generate_nurbs_patch`

```julia
generate_nurbs_patch(name::Symbol, nels::NTuple{dim,Int}, order; kwargs...)
```

- `nels`: number of elements in each parametric direction. Its length gives the parametric dimension.
- `order`: polynomial degree, an `Int` or a tuple with one degree per direction. (A `BezierGrid` requires the same degree in all directions.)

Common keyword arguments, where they apply:

- `multiplicity`: multiplicity of the interior knots in each direction (default 1, at most the degree). A higher multiplicity lowers the continuity between the elements.
- `sdim`: spatial dimension of the control points, for a patch that is embedded in a higher dimension.

### Boxes

| `name` | Dimension | Keyword arguments |
|:-------|:----------|:------------------|
| `:line` | 1D | `size` (required), `cornerpos = (0.0,)`, `multiplicity`, `sdim = 1` |
| `:rectangle` | 2D | `size` (required), `cornerpos = (0.0, 0.0)`, `multiplicity`, `sdim = 2` |
| `:cube` | 3D | `size = (2.0, 2.0, 2.0)`, `cornerpos = (-1.0, -1.0, -1.0)`, `multiplicity` |
| `:hypercube` | 1D, 2D or 3D | as `:line`, `:rectangle` or `:cube`, depending on the length of `nels` |

`cornerpos` is the corner with the smallest coordinates and `size` the side lengths.

```julia
generate_nurbs_patch(:rectangle, (8, 4), 2; cornerpos = (0.0, 0.0), size = (2.0, 1.0))
generate_nurbs_patch(:cube, (4, 4, 4), 2; cornerpos = (-1.0, -1.0, -1.0), size = (2.0, 2.0, 2.0))
```

### Curved patches

All keyword arguments in this table are required, except `multiplicity`.

| `name` | Geometry | `nels` | Keyword arguments |
|:-------|:---------|:-------|:------------------|
| `:singly_curved` | Arc with a thickness, in 2D | 2 | `α` (angle), `R` (radius), `thickness`, `multiplicity` |
| `:singly_curved` | Arc with a thickness and a width, in 3D | 3 | `α`, `R`, `thickness`, `width`, `multiplicity` |
| `:singly_curved_shell` | Cylindrical surface in 3D | 2 | `α`, `R`, `width`, `multiplicity` |
| `:singly_curved_beam` | Arc (a curve) in 2D | 1 | `α`, `R`, `multiplicity` |
| `:doubly_curved` | Surface with two radii in 3D | 2 | `r1`, `r2`, `α1`, `α2`, `multiplicity` |

### Circles and cylinders

These patches use rational (NURBS) basis functions of degree 2, so the circular shapes are exact.

| `name` | Geometry | Restrictions | Keyword arguments |
|:-------|:---------|:-------------|:------------------|
| `:plate_with_hole` | Quarter of a square plate with side 4 and a circular hole with radius 1 | `order = 2`, `nels[1]` even and at least 2 | none |
| `:ring` | Circular ring | `order = 2`, `nels = (4, 1)` | `ri`, `ro` (inner and outer radius, required) |
| `:cylinder_sector` | Solid half cylinder | `order = 2`, `nels[2]` even and greater than 2, `nels[3]` greater than 1 | `L` (length), `r` (radius), `x0 = -L/2` (start of the axis) |

```julia
generate_nurbs_patch(:plate_with_hole, (20, 10), 2)
generate_nurbs_patch(:ring, (4, 1), 2; ri = 1.0, ro = 2.0)
generate_nurbs_patch(:cylinder_sector, (4, 4, 2), 2; L = 2.0, r = 1.0)
```

## `generate_grid` with a `BezierCell`

Ferrite's `generate_grid` has methods for `BezierCell`s. They create a box with `generate_nurbs_patch`, convert it to a `BezierGrid`, and add facetsets for the sides.

| Call | Facetsets |
|:-----|:----------|
| `generate_grid(BezierCell{RefLine, p}, (nx,), left::Vec{1}, right::Vec{1})` | `"left"`, `"right"` |
| `generate_grid(BezierCell{RefQuadrilateral, p}, (nx, ny), LL, LR, UR, UL)` | `"left"`, `"right"`, `"bottom"`, `"top"` |
| `generate_grid(BezierCell{RefHexahedron, p}, (nx, ny, nz), left::Vec{3}, right::Vec{3})` | `"left"`, `"right"`, `"front"`, `"back"`, `"bottom"`, `"top"` |

For the quadrilateral, only the corners `LL` (lower left) and `UR` (upper right) are used; the grid is always an axis-aligned rectangle.

```julia
generate_grid(BezierCell{RefLine, 2}, (10,), Vec(0.0), Vec(1.0))
generate_grid(BezierCell{RefQuadrilateral, 2}, (8, 4), Vec(0.0, 0.0), Vec(2.0, 0.0), Vec(2.0, 1.0), Vec(0.0, 1.0))
generate_grid(BezierCell{RefHexahedron, 2}, (4, 4, 4), Vec(-1.0, -1.0, -1.0), Vec(1.0, 1.0, 1.0))
```
