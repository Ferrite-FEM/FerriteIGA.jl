# Bezier extraction

Bézier extraction is the technique that lets a finite element code such as Ferrite compute with splines. This page explains the idea step by step, starting with a curve. It assumes that you have read [B-splines and NURBS](@ref).

All figures on this page are computed with FerriteIGA (see `docs/generate_figures.jl`).

## Why extraction?

A finite element code is organized around the element. Every element has the same small set of shape functions, defined once on a reference element, and the code loops over the elements, evaluates these shape functions at the quadrature points, and assembles.

B-splines do not fit this pattern directly:

- A basis function of degree ``p`` stretches over ``p + 1`` elements, and it looks different on each of them.
- The pieces of the basis functions on an element depend on the knots around it, so they differ from element to element.

We would have to evaluate the Cox–de Boor recursion at every quadrature point of every element. Bézier extraction avoids this. It rewrites the spline basis on each element in terms of one fixed set of polynomials, the same on all elements. The rewriting is done by a small matrix per element, the **extraction operator**, which is computed once from the knot vector. After that, a spline element looks like an ordinary finite element to the rest of the code.

## Bernstein polynomials

The fixed set of polynomials is the **Bernstein basis**. On an element with the local coordinate ``t`` running from 0 to 1, the ``p + 1`` Bernstein polynomials of degree ``p`` are

```math
B_{a}(t) = \binom{p}{a-1} \, t^{a-1} (1 - t)^{p-a+1}, \qquad a = 1, \ldots, p+1 .
```

![Bernstein polynomials of degree 1, 2 and 3](figures/bezier_bernstein.svg)

*The Bernstein polynomials of degree 1, 2 and 3 on one element.*

They have much in common with B-splines: they are non-negative, they sum to one, and the first and the last one are interpolatory at the ends. In fact, they *are* B-splines, for the knot vector of a single element with both end knots repeated ``p + 1`` times, for example ``\{0, 0, 0, 1, 1, 1\}`` for ``p = 2``. And any polynomial of degree ``p`` on the element can be written as a combination of them.

That last property is all we need. On a single element, every B-spline is just a polynomial of degree ``p``, so it can be written as a combination of the Bernstein polynomials of that element.

## A spline on one element

Take the quadratic basis of the knot vector ``\{0, 0, 0, 1, 2, 3, 3, 3\}``, with five functions and three elements.

![Quadratic B-splines and the Bernstein polynomials of each element](figures/nurbs_bezier.svg)

*Top: the five quadratic B-splines. Bottom: the Bernstein polynomials, which are the same on every element.*

On the middle element (shaded), three B-splines are non-zero: ``N_2``, ``N_3`` and ``N_4``. Each of them is a combination of the three Bernstein polynomials of that element:

![The B-splines on the middle element as combinations of Bernstein polynomials](figures/bezier_decomposition.svg)

*The pieces of* ``N_2``*,* ``N_3`` *and* ``N_4`` *on the middle element (black), and the scaled Bernstein polynomials that add up to them (coloured). For* ``N_2`` *and* ``N_4`` *there is only one term, which coincides with the black curve.*

Written as a matrix,

```math
\begin{bmatrix} N_2 \\ N_3 \\ N_4 \end{bmatrix}
=
\begin{bmatrix}
\tfrac12 & 0 & 0 \\
\tfrac12 & 1 & \tfrac12 \\
0 & 0 & \tfrac12
\end{bmatrix}
\begin{bmatrix} B_1 \\ B_2 \\ B_3 \end{bmatrix}
\qquad \text{or} \qquad
\boldsymbol{N}^e = \boldsymbol{C}^e \boldsymbol{B} .
```

The matrix ``\boldsymbol{C}^e`` is the **Bézier extraction operator** of element ``e``. It has one row for each spline function that is non-zero on the element, and one column for each Bernstein polynomial. Three things about it are worth noting:

- It contains only numbers. It does not depend on where in the element we evaluate, so the same matrix is used at all quadrature points.
- It depends only on the knot vector, not on the control points or the geometry.
- Each column sums to one, because both the B-splines and the Bernstein polynomials sum to one.

## Where the numbers come from

The extraction operators are computed with knot insertion. Recall from [Smoothness and repeated knots](@ref) that a knot with multiplicity ``p`` makes the basis ``C^0`` and interpolatory there, and from [Refining without changing the shape](@ref) that inserting a knot does not change the spline: every old basis function is a combination of the new ones.

So insert every interior knot until it has multiplicity ``p``:

![Inserting knots until the basis is the Bernstein basis of each element](figures/bezier_knot_insertion.svg)

*Quadratic basis functions while the interior knots are repeated. In the last row, every interior knot appears twice, and the basis consists of the Bernstein polynomials of the three elements (neighbouring elements share the function that peaks at the knot between them).*

After the last insertion, the basis on each element is the Bernstein basis. The coefficients that express the old functions in the new ones are exactly the entries of the extraction operators. This is how `BezierGrid` computes them, with an algorithm that goes through the knot vector once and never forms the refined basis explicitly (Borden et al., 2011). In two and three dimensions, the operator of an element is the tensor (Kronecker) product of the one-dimensional operators of its knot spans.

In FerriteIGA, the extraction operators are stored in the `BezierGrid`, and [`get_extraction_operator`](@ref) returns the operator of a cell:

```@example bezier
using FerriteIGA

knots  = [0, 0, 0, 1, 2, 3, 3, 3.0]
points = [Vec(0.0, 0.0), Vec(0.5, 1.5), Vec(2.0, 2.0), Vec(3.5, 0.5), Vec(4.5, 1.5)]
curve  = NURBSMesh((knots,), (2,), points)
grid   = BezierGrid(curve)

C = get_extraction_operator(grid, 2)       # the middle element
natural = [1, 3, 2]                        # Ferrite numbers the two ends of a line first, then the interior
Matrix(reduce(hcat, C.C)')[natural, natural]
```

This is the matrix above. The operator is stored as one sparse row per spline function. Note the ordering: the nodes and shape functions of a cell follow Ferrite's convention (the ends or corners first, then edges, then the interior), so the rows and columns are permuted compared with the left-to-right order used in the text.

## Bézier control points

The extraction operator also gives a second description of the geometry. A spline curve is ``\boldsymbol{C}(\xi) = \sum_A N_A(\xi) \boldsymbol{P}_A``. On element ``e``, inserting ``\boldsymbol{N}^e = \boldsymbol{C}^e \boldsymbol{B}`` gives

```math
\boldsymbol{C}(\xi) = \sum_{a=1}^{p+1} B_a \, \boldsymbol{P}^b_a,
\qquad
\boldsymbol{P}^b_a = \sum_{A} C^e_{Aa} \, \boldsymbol{P}_A .
```

The points ``\boldsymbol{P}^b_a`` are the **Bézier control points** of the element. In matrix form they are ``(\boldsymbol{C}^e)^T`` times the spline control points of the element.

![A spline curve with its control points, and with the Bézier control points of each element](figures/bezier_curve.svg)

*Left: a quadratic spline curve with its 5 control points. Each element uses 3 of them, and neighbouring elements share 2. Right: the same curve described by 3 Bézier control points per element. The first and the last Bézier control point of an element lie on the curve, at the ends of the element.*

The two descriptions give exactly the same curve. The spline description is compact and smooth by construction. The Bézier description is local: each element is a self-contained polynomial piece, which is what a finite element code wants. The smoothness between the elements is not lost. It is built into the extraction operators, which make sure that the Bézier pieces fit together smoothly.

[`get_bezier_coordinates`](@ref) returns both sets of control points of a cell:

```@example bezier
xb, wb, x, w = get_bezier_coordinates(grid, 2)
x[natural]            # spline control points of the middle element: P₂, P₃, P₄
```

```@example bezier
xb[natural]           # its Bézier control points
```

The first Bézier control point is the midpoint between ``\boldsymbol{P}_2`` and ``\boldsymbol{P}_3``, as the first column ``(\tfrac12, \tfrac12, 0)`` of the extraction operator says.

## Surfaces

For a surface, the Bernstein basis of an element consists of the ``(p+1)^2`` products ``B_a(s) B_b(t)``, and each element has a ``(p+1) \times (p+1)`` net of Bézier control points.

![The control net of a NURBS patch and the Bézier control nets of its elements](figures/bezier_surface.svg)

*The quarter ring from the NURBS page, with 6 biquadratic elements. Left: the control net of the patch. Right: the Bézier control net of each element. The corners of each Bézier net lie on the corners of the element, and the nets of neighbouring elements share the points along their common edge.*

The picture to the right looks like a mesh of standard 9-node finite elements, and it can be used like one: the Bézier control points are the "nodes" that define the geometry of each element. This is also how FerriteIGA writes VTK files, since visualization programs know Bézier cells but not NURBS patches.

## Rational splines and extraction

NURBS have a weight ``w_A`` for each control point. The weights are extracted in the same way as the control points, which gives a Bézier weight for each Bézier control point,

```math
w^b_a = \sum_{A} C^e_{Aa} \, w_A .
```

The weight function in the denominator of the NURBS basis can then be evaluated with either basis, ``W = \sum_A w_A N_A = \sum_a w^b_a B_a``, and the rational basis functions on the element are

```math
R_A = \frac{w_A \, N_A}{W} = \frac{w_A \sum_a C^e_{Aa} B_a}{\sum_a w^b_a B_a} .
```

The Bézier control points of a rational spline are computed with the weights included,

```math
\boldsymbol{P}^b_a = \frac{1}{w^b_a} \sum_{A} C^e_{Aa} \, w_A \, \boldsymbol{P}_A ,
```

and the geometry is ``\sum_A R_A \boldsymbol{P}_A = \sum_a w^b_a B_a \boldsymbol{P}^b_a / W`` in both descriptions.

## T-splines and extraction

Nothing in the idea requires a tensor product patch. On a Bézier element of a T-spline, each basis function that is non-zero there is a polynomial, and can be written in the Bernstein basis of the element. The extraction operator again has one row per basis function. For T-splines the rows are computed one basis function at a time, from the local knot vectors of the function (Scott et al., 2011). This is what makes it possible to use T-splines with the same element routines as NURBS; see [T-splines](@ref).

## Extraction in FerriteIGA

The pieces map to the code as follows.

| Concept | In FerriteIGA |
|:--------|:--------------|
| Bernstein basis of an element | `IGAInterpolation{shape, p}()`, with `shape` one of `RefLine`, `RefQuadrilateral`, `RefHexahedron` |
| Extraction operators, one per cell | computed by `BezierGrid(mesh)`, returned by [`get_extraction_operator`](@ref) |
| Spline and Bézier control points and weights of a cell | `getcoordinates(grid, cellid)`, a [`BezierCoords`](@ref) with the fields `x`, `w` (spline) and `xb`, `wb` (Bézier), and the operator `beo` |
| Applying the operator to the shape functions | `reinit!(cv, coords)` for [`BezierCellValues`](@ref), `reinit!(fv, coords, facetid)` for [`BezierFacetValues`](@ref) |

`BezierCellValues` evaluate the Bernstein polynomials once, at the quadrature points of the reference element, exactly as Ferrite's `CellValues` do for Lagrange polynomials. `reinit!` then multiplies them with the extraction operator and the weights of the cell, which gives the values and gradients of the NURBS basis functions:

```@example bezier
patch = generate_nurbs_patch(:plate_with_hole, (20, 10), 2)
grid = BezierGrid(patch)

ip = IGAInterpolation{RefQuadrilateral, 2}()          # the Bernstein basis
cv = BezierCellValues(QuadratureRule{RefQuadrilateral}(4), ip)

coords = getcoordinates(grid, 1)                      # control points, weights and extraction operator of cell 1
reinit!(cv, coords)
sum(shape_value(cv, 1, i) for i in 1:getnbasefunctions(cv))     # the NURBS functions sum to one
```

After `reinit!`, `shape_value`, `shape_gradient` and `getdetJdV` refer to the NURBS basis and the exact geometry. The degrees of freedom of the cell are the coefficients of the spline functions, i.e. they belong to the spline control points, not to the Bézier control points. From here on, the `DofHandler`, the `ConstraintHandler` and the assembly work as described in the [Ferrite manual](https://ferrite-fem.github.io/Ferrite.jl/stable/). A complete example is [Infinite plate with hole](@ref).

An ordinary Ferrite grid can be wrapped with `BezierGrid(grid)`. Its cells have no extraction operator, and are treated as standard finite elements.

## Summary

- On each element, the spline functions are combinations of the Bernstein polynomials of that element: ``\boldsymbol{N}^e = \boldsymbol{C}^e \boldsymbol{B}``.
- The extraction operator ``\boldsymbol{C}^e`` is a small matrix of numbers. It is computed once from the knot vectors, by repeating the interior knots.
- The transpose of the operator turns the spline control points (and weights) of an element into its Bézier control points (and weights).
- With the Bernstein polynomials as shape functions and the extraction operator applied in `reinit!`, a spline element is used like a standard finite element.

## Further reading

- M. J. Borden, M. A. Scott, J. A. Evans, and T. J. R. Hughes. Isogeometric finite element data structures based on Bézier extraction of NURBS. *Int. J. Numer. Meth. Engng.*, 87:15–47, 2011. [doi:10.1002/nme.2968](https://doi.org/10.1002/nme.2968)
- M. A. Scott, M. J. Borden, C. V. Verhoosel, T. W. Sederberg, and T. J. R. Hughes. Isogeometric finite element data structures based on Bézier extraction of T-splines. *Int. J. Numer. Meth. Engng.*, 88:126–156, 2011. [doi:10.1002/nme.3167](https://doi.org/10.1002/nme.3167)
