```@meta
CurrentModule = FerriteIGA
DocTestSetup = :(using FerriteIGA)
```

# BezierCellValues

[`BezierCellValues`](@ref) and [`BezierFacetValues`](@ref) are the isogeometric
counterparts of Ferrite's `CellValues` and `FacetValues`. After a call to
[`reinit!`](@ref reinit_bezier) they behave exactly like their Ferrite counterparts, so all of
Ferrite's functions for `CellValues` and `FacetValues` (`getnquadpoints`, `getdetJdV`,
`shape_value`, `shape_gradient`, `shape_hessian`, `function_value`, `function_gradient`,
`spatial_coordinate`, `getnormal`, ...) can be used, see the
[Ferrite FEValues reference](https://ferrite-fem.github.io/Ferrite.jl/stable/reference/fevalues/).
For details about how the values are computed, see the
[BezierCellValues topic guide](@ref beziercellvalues_topicguide).

## Interpolation

```@docs
IGAInterpolation
```

## Cell and facet values

```@docs
BezierCellValues
BezierFacetValues
```

## [Reinitialization](@id reinit_bezier)

The values are updated for a cell with Ferrite's `reinit!` function, which has the following
methods for [`BezierCellValues`](@ref) and [`BezierFacetValues`](@ref):

```julia
reinit!(cv::BezierCellValues, bc::BezierCoords)
reinit!(cv::BezierCellValues, (xb, wb))
reinit!(fv::BezierFacetValues, bc::BezierCoords, facetnr::Int)
reinit!(fv::BezierFacetValues, (xb, wb), facetnr::Int)
```

The first form, with a [`BezierCoords`](@ref) obtained from `getcoordinates(grid, cellid)`,
is the recommended one. The second form takes the Bézier control points `xb` and Bézier
weights `wb` directly, and requires that the extraction operator and NURBS weights of the
cell have been set beforehand with [`set_bezier_operator!`](@ref).

```@docs
set_bezier_operator!
```
