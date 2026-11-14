# NURBSMesh

`NURBSMesh` represents one tensor-product B-spline or NURBS patch. It stores
the knot vectors, polynomial orders, control points, weights, and the
connectivity needed to visit each non-empty knot span. See [B-splines and NURBS](@ref) for an introduction.

## Refinement of a patch

A patch can be refined in three ways. None of them changes the geometry, and each acts on one parametric direction `dir`.

| Refinement | Function | Effect |
|:-----------|:---------|:-------|
| h-refinement (knot insertion) | [`knotinsertion`](@ref) | Inserts the knot `ξ` once. If `ξ` is a new value, there is one more element, joined with ``C^{p-1}`` continuity. If `ξ` is already a knot, its multiplicity ``m`` increases and the continuity there drops to ``C^{p-m}``. |
| p-refinement (order elevation) | [`orderelevation`](@ref) | Raises the degree by one. The multiplicity of all knots increases by one, so the continuity at the knots is kept. |
| k-refinement | [`smoothnesselevation`](@ref) | Raises the degree by one and then inserts `new_knots`. The new knots have multiplicity one at the elevated degree, so the continuity there is ``C^{p}`` (with ``p`` the original degree). |

These functions return a new mesh and leave the input unchanged:

```julia
mesh = generate_nurbs_patch(:hypercube, (1, 1), (2, 2); cornerpos = (0.0, 0.0), size = (2.0, 3.0))

mesh = knotinsertion(mesh, 0.0; dir = 1)     # h-refinement: insert ξ = 0 in the first direction
mesh = orderelevation(mesh; dir = 1)         # p-refinement: raise the degree in the first direction
grid = BezierGrid(mesh)
```

The versions ending in `!` work in place on the knot vectors, control points and weights, which avoids building a mesh after every step. Copy the arrays first, and build a new `NURBSMesh` at the end so that the element connectivity is up to date. `orderelevation!` and `smoothnesselevation!` return the new orders.

```julia
kv = map(copy, mesh.knot_vectors)
cp = copy(mesh.control_points)
w  = copy(mesh.weights)
orders = mesh.orders

knotinsertion!(kv, orders, cp, w, 0.5; dir = 1)
orders = orderelevation!(kv, orders, cp, w; dir = 2)
orders = smoothnesselevation!(kv, orders, cp, w, [0.25, 0.75]; dir = 1)

mesh = NURBSMesh(kv, orders, cp, w)
```

## API

```@docs
NURBSMesh
parent_to_parametric_map
eval_parametric_coordinate
knotinsertion
orderelevation
smoothnesselevation
knotinsertion!
orderelevation!
smoothnesselevation!
```
