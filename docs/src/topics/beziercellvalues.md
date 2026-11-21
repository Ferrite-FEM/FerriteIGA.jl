```@setup cellvalues
using FerriteIGA
```

# [BezierCellValues](@id beziercellvalues_topicguide)

In Ferrite, `CellValues` and `FacetValues` are used inside the element routines to query
the integration weights, shape function values and gradients, and much more (see the
[Ferrite FEValues topic guide](https://ferrite-fem.github.io/Ferrite.jl/stable/topics/FEValues/)).
FerriteIGA provides the isogeometric counterparts, [`BezierCellValues`](@ref) and
[`BezierFacetValues`](@ref). They have exactly the same interface: after a call to `reinit!`
for a cell, the functions `getdetJdV`, `shape_value`, `shape_gradient`, `shape_hessian`,
`function_value`, `function_gradient`, `spatial_coordinate`, `getnormal`, ... return the
values of the B-spline or NURBS basis functions of the patch on that cell. Element routines
written for Ferrite can therefore be used without modification.

This topic guide describes how the values are computed, i.e. what happens in `reinit!`, and
how the different variants of `reinit!` are used. It builds on the Bézier extraction
described in the [BezierGrid](@ref beziergrid_topicguide) topic guide.

## The reference basis: `IGAInterpolation`

In Ferrite, the reference shape functions are given by an interpolation, e.g.
`Lagrange{RefQuadrilateral, 2}()`. In FerriteIGA, the interpolation is an
[`IGAInterpolation`](@ref), whose reference shape functions are the Bernstein polynomials of
the given degree on the reference cell ``[-1, 1]^d``:

```@example cellvalues
ip = IGAInterpolation{RefQuadrilateral, 2}()
getnbasefunctions(ip)
```

The Bernstein polynomials are ordered like the nodes of a [`BezierCell`](@ref), with the
vertices first. As for Lagrange interpolations, vector-valued interpolations are created with
`ip^dim`. The order of the interpolation must match the order of the cells of the grid.

## Construction

A [`BezierCellValues`](@ref) is constructed from a quadrature rule and an interpolation,
exactly like a `CellValues`. As in Ferrite, the reference values of the Bernstein polynomials
(and their derivatives) at the quadrature points are computed once, at construction:

```@example cellvalues
qr = QuadratureRule{RefQuadrilateral}(3)
cv = BezierCellValues(qr, ip)
```

The keyword arguments `update_gradients`, `update_hessians` and `update_detJdV` control
which values are computed in `reinit!`, as for `CellValues`. For example, the
[Kirchhoff plate](@ref kirchhoff_plate) tutorial needs the second derivatives of the basis
functions, and therefore uses `update_hessians = true`. A Ferrite `CellValues` created with an
`IGAInterpolation` is also a `BezierCellValues`:

```@example cellvalues
CellValues(qr, ip) isa BezierCellValues
```

!!! note "Quadrature"
    NURBS basis functions are rational functions, so Gauss quadrature is not exact for them
    (for B-splines on straight-sided patches it is, as for polynomials). In practice, a Gauss
    rule with ``p + 1`` points per direction is the usual choice, as for standard finite
    elements of degree ``p``.

## [Reinitialization: from Bernstein polynomials to NURBS](@id beziercellvalues-reinit)

Calling `reinit!(cv, getcoordinates(grid, cellid))` updates the values for cell `cellid`. The
[`BezierCoords`](@ref) returned by `getcoordinates` provides the extraction operator
``\boldsymbol{C}^e``, the NURBS weights ``w_A``, the Bézier control points
``\boldsymbol{X}^b_a`` and the Bézier weights ``w^b_a`` of the cell. With these, `reinit!`
performs the following steps in each quadrature point ``\boldsymbol{\xi}_q``.

**1. Bézier extraction.** The pre-computed Bernstein values and derivatives are transformed
with the extraction operator and multiplied with the NURBS weights, giving the weighted
B-spline functions

```math
\tilde{N}_A(\boldsymbol{\xi}) = w_A \sum_{a} C^e_{Aa} B_a(\boldsymbol{\xi}), \qquad
\frac{\partial \tilde{N}_A}{\partial \boldsymbol{\xi}} = w_A \sum_{a} C^e_{Aa} \frac{\partial B_a}{\partial \boldsymbol{\xi}}, \qquad
\frac{\partial^2 \tilde{N}_A}{\partial \boldsymbol{\xi}^2} = w_A \sum_{a} C^e_{Aa} \frac{\partial^2 B_a}{\partial \boldsymbol{\xi}^2}.
```

**2. Rational basis.** The weight function and its derivatives are computed from the Bézier
weights,

```math
W = \sum_{a} w^b_a B_a, \qquad
\frac{\partial W}{\partial \boldsymbol{\xi}} = \sum_{a} w^b_a \frac{\partial B_a}{\partial \boldsymbol{\xi}}, \qquad
\frac{\partial^2 W}{\partial \boldsymbol{\xi}^2} = \sum_{a} w^b_a \frac{\partial^2 B_a}{\partial \boldsymbol{\xi}^2},
```

and the NURBS basis functions ``R_A = \tilde{N}_A / W`` and their derivatives with respect to
the reference coordinates follow from the quotient rule,

```math
\begin{align*}
\frac{\partial R_A}{\partial \boldsymbol{\xi}} &= \frac{1}{W} \frac{\partial \tilde{N}_A}{\partial \boldsymbol{\xi}} - \frac{\tilde{N}_A}{W^2} \frac{\partial W}{\partial \boldsymbol{\xi}}, \\
\frac{\partial^2 R_A}{\partial \boldsymbol{\xi}^2} &= \frac{1}{W} \frac{\partial^2 \tilde{N}_A}{\partial \boldsymbol{\xi}^2}
- \frac{1}{W^2} \left[ \frac{\partial \tilde{N}_A}{\partial \boldsymbol{\xi}} \otimes \frac{\partial W}{\partial \boldsymbol{\xi}} + \frac{\partial W}{\partial \boldsymbol{\xi}} \otimes \frac{\partial \tilde{N}_A}{\partial \boldsymbol{\xi}} + \tilde{N}_A \frac{\partial^2 W}{\partial \boldsymbol{\xi}^2} \right]
+ \frac{2 \tilde{N}_A}{W^3} \frac{\partial W}{\partial \boldsymbol{\xi}} \otimes \frac{\partial W}{\partial \boldsymbol{\xi}}.
\end{align*}
```

**3. Geometry mapping.** The geometry of the cell is the rational Bézier map
``\boldsymbol{x}(\boldsymbol{\xi}) = \sum_a \hat{R}_a(\boldsymbol{\xi}) \boldsymbol{X}^b_a``, with
``\hat{R}_a = w^b_a B_a / W`` (see the [BezierGrid](@ref beziergrid_topicguide) topic guide).
Its Jacobian ``\boldsymbol{J}`` (and, if hessians are requested, its second derivative
``\boldsymbol{\mathcal{H}}``) is computed from the Bézier control points and weights,

```math
\boldsymbol{J} = \frac{\partial \boldsymbol{x}}{\partial \boldsymbol{\xi}} = \sum_a \boldsymbol{X}^b_a \otimes \frac{\partial \hat{R}_a}{\partial \boldsymbol{\xi}},
\qquad
\boldsymbol{\mathcal{H}} = \frac{\partial^2 \boldsymbol{x}}{\partial \boldsymbol{\xi}^2} = \sum_a \boldsymbol{X}^b_a \otimes \frac{\partial^2 \hat{R}_a}{\partial \boldsymbol{\xi}^2},
```

where the derivatives of ``\hat{R}_a`` are given by the quotient rule as above. The volume
element is ``\mathrm{d}\Omega = \det(\boldsymbol{J}) \, w_q`` for the quadrature weight
``w_q``, which is what `getdetJdV` returns.

**4. Mapping to the physical cell.** Finally, the derivatives are mapped to the physical
coordinates exactly as for standard finite elements in Ferrite (the identity mapping, see the
[Ferrite FEValues topic guide](https://ferrite-fem.github.io/Ferrite.jl/stable/topics/FEValues/)),

```math
\frac{\partial R_A}{\partial \boldsymbol{x}} = \frac{\partial R_A}{\partial \boldsymbol{\xi}} \cdot \boldsymbol{J}^{-1},
\qquad
\frac{\partial^2 R_A}{\partial \boldsymbol{x}^2} = \boldsymbol{J}^{-\mathrm{T}} \cdot \left[ \frac{\partial^2 R_A}{\partial \boldsymbol{\xi}^2} - \frac{\partial R_A}{\partial \boldsymbol{x}} \cdot \boldsymbol{\mathcal{H}} \right] \cdot \boldsymbol{J}^{-1}.
```

!!! details "Why the Bézier weights?"
    The weight function can be computed either from the NURBS weights and the B-spline
    functions, ``W = \sum_A w_A N_A``, or from the Bézier weights and the Bernstein
    polynomials, ``W = \sum_a w^b_a B_a``. The two are identical, since
    ``\sum_A w_A N_A = \sum_A w_A \sum_a C^e_{Aa} B_a = \sum_a (\boldsymbol{C}^{e\mathrm{T}} \boldsymbol{w})_a B_a``.
    Similarly, the geometry is the same whether it is computed from the NURBS control points
    and the NURBS basis, or from the Bézier control points and the rational Bernstein basis.
    Using the Bézier form for the geometry makes each cell self-contained: it only requires
    the Bernstein polynomials and the Bézier control points and weights of the cell, which is
    also how the cells are exported to VTK (see the [Export](@ref export_topicguide) topic
    guide).

### Example

To illustrate, we reinitialize a `BezierCellValues` for a cell of the plate with a hole, and
check some properties of the NURBS basis:

```@example cellvalues
grid = BezierGrid(generate_nurbs_patch(:plate_with_hole, (4, 2), 2))
coords = getcoordinates(grid, 1)
reinit!(cv, coords)
nothing # hide
```

The basis functions form a partition of unity, so the values sum to one and the gradients sum
to zero in every quadrature point:

```@example cellvalues
q_point = 1
sum(i -> shape_value(cv, q_point, i), 1:getnbasefunctions(cv))
```

```@example cellvalues
sum(i -> shape_gradient(cv, q_point, i), 1:getnbasefunctions(cv))
```

The physical position of a quadrature point is obtained with `spatial_coordinate`, which
takes the `BezierCoords` of the cell:

```@example cellvalues
x = spatial_coordinate(cv, q_point, coords)
```

Summing the volume elements over all cells gives the area of the domain,
``16 - \pi/4``, up to the quadrature error for the rational functions:

```@example cellvalues
area = 0.0
for cellid in 1:getncells(grid)
    reinit!(cv, getcoordinates(grid, cellid))
    for q_point in 1:getnquadpoints(cv)
        global area += getdetJdV(cv, q_point)
    end
end
area - (16 - π / 4)
```

## Variants of `reinit!`

There are two ways of reinitializing a `BezierCellValues` for a cell:

| Call | Description |
|:-----|:------------|
| `reinit!(cv, coords::BezierCoords)` | Recommended. `coords` is obtained with `getcoordinates(grid, cellid)`, or from an [`IGACellCache`](@ref). |
| `reinit!(cv, (xb, wb))` | Takes the Bézier control points and weights directly. The extraction operator and NURBS weights of the cell must be set before, with [`set_bezier_operator!`](@ref). |

The second form is useful in hand-written loops that use pre-allocated arrays, as in the
[Infinite plate with hole](@ref tutorial-plate-with-hole) tutorial:

```@example cellvalues
n = Ferrite.nnodes_per_cell(grid, 1)
xb, wb = zeros(Vec{2, Float64}, n), zeros(n)
x, w = zeros(Vec{2, Float64}, n), zeros(n)
for cellid in 1:getncells(grid)
    get_bezier_coordinates!(xb, wb, x, w, grid, cellid)
    set_bezier_operator!(cv, get_extraction_operator(grid, cellid), w)
    reinit!(cv, (xb, wb))
    # ... element routine
end
```

!!! warning
    `reinit!` requires the cell to have an extraction operator. For cells without one (e.g.
    standard Ferrite cells in a mixed [`BezierGrid`](@ref)), use a standard `CellValues`
    with the node coordinates, `reinit!(cv, coords.x)`.

## `BezierFacetValues`

[`BezierFacetValues`](@ref) is the counterpart of Ferrite's `FacetValues`, and is
constructed from a `FacetQuadratureRule`. It is reinitialized for a facet with
`reinit!(fv, coords, facetnr)`. In addition to the values above, `getnormal` returns the
outward unit normal of the (curved) facet in each quadrature point.

As an example, we compute the length of the boundary of the hole in the plate, which should
be a quarter of the circumference of the unit circle. The hole is facet 1 of the cells in
the first row (see the [BezierGrid](@ref beziergrid_topicguide) topic guide):

```@example cellvalues
fqr = FacetQuadratureRule{RefQuadrilateral}(3)
fv = BezierFacetValues(fqr, ip)
hole = [FacetIndex(cellid, 1) for cellid in 1:4]
len = 0.0
for (cellid, facetnr) in hole
    reinit!(fv, getcoordinates(grid, cellid), facetnr)
    for q_point in 1:getnquadpoints(fv)
        global len += getdetJdV(fv, q_point)
    end
end
len - π / 2
```

Since the geometry is exact, the remaining (small) difference is the quadrature error, which
decreases rapidly with the number of quadrature points or elements.

## Iterating over cells and facets

Ferrite's `CellIterator` and `FacetIterator` do not support `BezierCellValues`, since they
provide the node coordinates instead of a [`BezierCoords`](@ref). The simplest alternative is
to loop over the cell indices and call `getcoordinates(grid, cellid)`. To avoid the
allocations of `getcoordinates`, use an [`IGACellCache`](@ref) (or [`IGAFaceCache`](@ref)
for facets), which stores pre-allocated data for one cell and is used like Ferrite's
`CellCache`:

```julia
cc = IGACellCache(dh)
for cellid in 1:getncells(grid)
    reinit!(cc, cellid)
    reinit!(cv, getcoordinates(cc))
    assemble_element!(ke, cv)
    assemble!(assembler, celldofs(cc), ke)
end
```

```julia
fc = IGAFaceCache(dh)
for facet in getfacetset(grid, "left")
    reinit!(fc, facet)
    reinit!(fv, getcoordinates(fc), facet[2])
    # ...
end
```
