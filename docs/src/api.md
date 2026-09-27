# API

API of FerriteIGA

After `reinit!` on cell or facet values, assembly follows the
[Ferrite manual](https://ferrite-fem.github.io/Ferrite.jl/stable/).

## Interpolation and cells

```@docs
IGAInterpolation
BezierCell
BezierCoords
BezierGrid
BezierExtractionOperator
```

## Grid helpers

Coordinates, weights, and extraction operators on a [`BezierGrid`](@ref).

```@docs
getweights!
get_nurbs_weights
get_nurbs_coordinates
get_bezier_coordinates
get_bezier_coordinates!
get_extraction_operator
```

## Cell and facet values

```@docs
BezierCellValues
BezierFacetValues
set_bezier_operator!
```

## Iteration caches

Ferrite-style caches for looping over cells and facets on a [`BezierGrid`](@ref).

```@docs
IGACellCache
IGAFaceCache
```

## I/O

```@docs
VTKIGAFile
```

## NURBS patches

Mesh construction and evaluation. Refinement APIs live on the [NURBSMesh](@ref) page.

```@docs
generate_nurbs_patch
```
