# NURBSMesh

`NURBSMesh` represents one tensor-product B-spline or NURBS patch. It stores
the knot vectors, polynomial orders, control points, weights, and the
connectivity needed to visit each non-empty knot span.

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