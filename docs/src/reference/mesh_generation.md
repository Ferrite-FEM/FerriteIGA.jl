```@meta
CurrentModule = FerriteIGA
DocTestSetup = :(using FerriteIGA)
```

# Mesh generation

FerriteIGA can generate a number of simple, single-patch geometries. The result is a
[`NURBSMesh`](@ref), which is converted to a [`BezierGrid`](@ref) before analysis. See the
[NURBSMesh topic guide](@ref nurbsmesh_topicguide) for how a `NURBSMesh` is refined and
converted.

Straight lines, rectangles and boxes can also be generated directly as a [`BezierGrid`](@ref)
with Ferrite's `generate_grid` and a [`BezierCell`](@ref), e.g.
`generate_grid(BezierCell{RefQuadrilateral, 2}, (4, 2))`.

```@docs
generate_nurbs_patch
```
