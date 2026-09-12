# Meshes

Geometry starts as a `NURBSMesh`. Convert it with `BezierGrid(mesh)` before assembly. Straight boxes can also use `generate_grid` with a `BezierCell`. Refinement is on the [Splines](@ref) page.

## `generate_nurbs_patch`

```julia
patch = generate_nurbs_patch(name, nels, order; kwargs...)
```

`name` is a symbol. `nels` is how many elements in each direction. `order` is the degree, or a tuple if the degree differs by direction.

### Boxes

`:line`, `:rectangle`, `:cube`, and `:hypercube` are straight patches. `:line` is 1D, `:rectangle` is 2D, `:cube` is 3D. `:hypercube` covers all three.

`cornerpos` is the lower corner and `size` is the side lengths. `multiplicity` repeats interior knots. `sdim` puts a 1D or 2D patch into 2D or 3D.

```julia
generate_nurbs_patch(:rectangle, (8, 4), 2; cornerpos=(0.0, 0.0), size=(2.0, 1.0))
generate_nurbs_patch(:cube, (4, 4, 4), 2; cornerpos=(-1.0, -1.0, -1.0), size=(2.0, 2.0, 2.0))
```

The same boxes with `generate_grid`:

```julia
generate_grid(BezierCell{RefLine, 2}, (10,), Vec(0.0), Vec(1.0))
generate_grid(BezierCell{RefQuadrilateral, 2}, (8, 4), Vec(0.0, 0.0), Vec(2.0, 0.0), Vec(2.0, 1.0), Vec(0.0, 1.0))
generate_grid(BezierCell{RefHexahedron, 2}, (4, 4, 4), Vec(-1.0, -1.0, -1.0), Vec(1.0, 1.0, 1.0))
```

### Curved patches

`:singly_curved` is an arc with thickness. In 2D the keywords are `α` (angle), `R` (radius), and `thickness`. In 3D add `width`.

`:singly_curved_shell` is a 3D surface (no thickness) with `α`, `R`, and `width`.

`:singly_curved_beam` is a 1D arc in the plane, with `α` and `R`.

`:doubly_curved` is a 3D surface with two radii and two angles, `r1`, `r2`, `α1`, `α2`.

`:doubly_curved_nurbs` is a quadratic NURBS surface (`r1`, `r2`, `α2`). The degree is fixed.

### Circles and cylinders

`:plate_with_hole` is a quarter plate with a circular hole of radius 1, outer size 4. The degree is 2. The first entry of `nels` must be even and at least 2.

```julia
generate_nurbs_patch(:plate_with_hole, (20, 10), 2)
```

`:ring` is a circular ring. The element count is fixed at `(4, 1)` and the degree at 2. Set inner and outer radii with `ri` and `ro`.

```julia
generate_nurbs_patch(:ring, (4, 1), 2; ri=1.0, ro=2.0)
```

`:cylinder_sector` is a half-cylinder solid. The degree is 2. `L` is the length, `r` the radius, and `x0` the start of the axis (default `-L/2`). `nels[2]` (around the circle) must be even and greater than 2. `nels[3]` must be greater than 1.

```julia
generate_nurbs_patch(:cylinder_sector, (4, 4, 2), 2; L=2.0, r=1.0)
```

## Parameter maps

`eval_parametric_coordinate(mesh, ξ)` maps a parameter `ξ` to a physical point. `parent_to_parametric_map(mesh, cellid, xi)` maps a point on the reference cell into the knot vector.
