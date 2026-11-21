```@meta
CurrentModule = FerriteIGA
DocTestSetup = :(using FerriteIGA)
```

# NURBSMesh

A [`NURBSMesh`](@ref) is a single tensor-product B-spline or NURBS patch. For more details
about the data structure and how to refine it, see the
[NURBSMesh topic guide](@ref nurbsmesh_topicguide), and for the underlying spline theory the
[Splines topic guide](@ref splines_topicguide).

## Main type

```@docs
NURBSMesh
```

## Evaluation

```@docs
eval_parametric_coordinate
parent_to_parametric_map
```

## Refinement

The refinement functions come in two variants. The non-mutating functions take a
[`NURBSMesh`](@ref) and return a new, refined mesh. The mutating functions (ending with `!`)
work directly on the arrays that define a patch, which is useful when several refinement
steps are applied before the mesh is (re)built.

```@docs
knotinsertion
orderelevation
smoothnesselevation
knotinsertion!
orderelevation!
smoothnesselevation!
```
