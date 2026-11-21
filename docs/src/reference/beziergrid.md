```@meta
CurrentModule = FerriteIGA
DocTestSetup = :(using FerriteIGA)
```

# BezierGrid

A [`BezierGrid`](@ref) is a Ferrite grid (`Ferrite.AbstractGrid`) in which the nodes are the
control points of a spline patch, and every cell additionally stores its Bézier extraction
operator. For more details see the [BezierGrid topic guide](@ref beziergrid_topicguide).

## Grid and cells

```@docs
BezierGrid
BezierCell
BezierExtractionOperator
```

## Cell coordinates

The function `getcoordinates(grid::BezierGrid, cellid)` returns a [`BezierCoords`](@ref),
which holds everything that is needed to [`reinit!`](@ref reinit_bezier) a
[`BezierCellValues`](@ref). The functions below give access to the individual parts.

```@docs
BezierCoords
get_nurbs_coordinates
get_nurbs_weights
get_nurbs_weights!
get_bezier_coordinates
get_bezier_coordinates!
get_extraction_operator
```

## Iteration caches

```@docs
IGACellCache
IGAFaceCache
```
