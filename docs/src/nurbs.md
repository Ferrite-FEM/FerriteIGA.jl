# B-splines and NURBS

This page is a step-by-step introduction to the functions that isogeometric analysis is built on. It starts with a single curve in one parameter and works up to surfaces. No prior knowledge of splines is needed, only some familiarity with finite elements. The [Splines](@ref) page is the short reference version of the same material, and the [T-splines](@ref) tutorial continues where this page ends.

All figures on this page are computed with FerriteIGA (see `docs/generate_figures.jl`).

## Why splines?

In a standard finite element method, the solution on each element is a polynomial, and neighbouring elements are joined so that the solution is continuous. Its derivatives, however, jump from one element to the next.

A *spline* is also a piecewise polynomial, but the pieces are joined more smoothly: not only the function, but also its first derivatives (and possibly higher ones) are continuous across the element boundaries. Splines are what CAD programs use to describe curves and surfaces. Isogeometric analysis uses the same functions for the geometry and for the unknown field, so the analysis works on the exact CAD geometry and gets a smooth basis as a bonus.

Three ingredients define a spline curve:

1. a polynomial **degree** ``p``,
2. a **knot vector**, which says where the polynomial pieces begin and end,
3. a list of **control points**, which gives the curve its shape.

We look at them one at a time.

## The knot vector

A knot vector is a list of non-decreasing numbers,

```math
\Xi = \{\xi_1, \xi_2, \ldots, \xi_{n+p+1}\}, \qquad \xi_1 \le \xi_2 \le \cdots \le \xi_{n+p+1}.
```

The numbers ``\xi_i`` are called *knots*. They are coordinates on the *parameter* axis ``\xi``, the axis along which we travel when we trace out the curve. The knots cut the parameter axis into intervals, called *knot spans*. Each knot span with a positive length is one **element**: on it, the spline is a single polynomial.

Two things are worth remembering from the start:

- A knot vector with ``n + p + 1`` knots gives ``n`` basis functions of degree ``p``. Read the other way: the number of basis functions is the number of knots minus ``p + 1``.
- Knots may be repeated. A repeated knot does not create a new element (the span between two equal knots has zero length), but it changes how smoothly the pieces are joined. We come back to this in [Smoothness and repeated knots](@ref).

## Building the basis functions

The basis functions of a spline are called **B-splines**. They are built degree by degree with the *Cox–de Boor recursion*. The recursion starts with degree zero, where each function is simply one on a single knot span and zero everywhere else,

```math
N_{i,0}(\xi) =
\begin{cases}
1 & \text{if } \xi_i \leq \xi < \xi_{i+1}, \\
0 & \text{otherwise}.
\end{cases}
```

A function of degree ``p`` is then a blend of two neighbouring functions of degree ``p-1``,

```math
N_{i,p}(\xi) = \frac{\xi - \xi_i}{\xi_{i+p} - \xi_i} \, N_{i,p-1}(\xi) + \frac{\xi_{i+p+1} - \xi}{\xi_{i+p+1} - \xi_{i+1}} \, N_{i+1,p-1}(\xi).
```

The first fraction ramps up from zero to one over the support of ``N_{i,p-1}``, and the second ramps down from one to zero over the support of ``N_{i+1,p-1}``. (If a denominator is zero, which happens for repeated knots, that term is dropped.)

The figure below shows what the recursion does for the knot vector ``\{0, 1, 2, 3, 4, 5\}``. Two box functions blend into a hat function, and two hat functions blend into a smooth quadratic bump.

![The Cox-de Boor recursion for degrees 0, 1 and 2](figures/nurbs_cox_de_boor.svg)

*Each row is built from the row above it. Every step up in degree makes the functions one knot span wider and one derivative smoother, and there is one function fewer.*

Reading the figure from top to bottom shows the main properties of B-splines:

- **Local support.** ``N_{i,p}`` is non-zero only between the knots ``\xi_i`` and ``\xi_{i+p+1}``, that is, on ``p + 1`` knot spans.
- **Non-negative.** The functions are never negative.
- **Smooth.** With distinct knots, a function of degree ``p`` has ``p - 1`` continuous derivatives at the knots (it is ``C^{p-1}``). The hats are continuous (``C^0``), the quadratic bumps also have a continuous slope (``C^1``).
- **Overlap.** On each element, exactly ``p + 1`` functions are non-zero. In the bottom row, the element between ``\xi_3`` and ``\xi_4`` is the only one that is covered by all three functions, and this is also the only element where they sum to one.

That last observation is a problem: near the ends of the knot vector there are too few functions. The standard cure is the open knot vector.

## The open knot vector

A knot vector is **open** if the first and the last knot are repeated ``p + 1`` times. This is the standard choice in CAD and in isogeometric analysis, and it is what `generate_nurbs_patch` creates. The figure shows the quadratic basis for the open knot vector

```math
\Xi = \{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5\}.
```

There are 11 knots, so there are ``11 - (2 + 1) = 8`` basis functions.

![Quadratic B-spline basis for an open knot vector](figures/nurbs_open_basis.svg)

*The 8 quadratic basis functions of the open knot vector* ``\{0,0,0,1,2,3,4,4,5,5,5\}``*. There are 5 elements: the knot spans between 0, 1, 2, 3, 4 and 5.*

With an open knot vector:

- The functions form a **partition of unity** on the whole domain: at every ``\xi``, they sum to one. This is needed for the functions to represent constants (rigid body motions, for example).
- The first and the last function equal one at the ends, and all other functions are zero there. The basis is *interpolatory* at the ends, like a finite element basis. This makes it easy to apply boundary conditions.

In the interior, the functions are not interpolatory: ``N_3`` never reaches one, and at its peak ``N_2`` and ``N_4`` are also non-zero. This is the price of smoothness, and the reason why the coefficients of a spline are called *control points* and not nodal values.

The figure also contains a surprise at ``\xi = 4``, where ``N_6`` has a sharp peak. The knot 4 appears twice in the knot vector, and that is the topic of the next section.

## Smoothness and repeated knots

Every time a knot is repeated, the basis loses one continuous derivative at that knot. If a knot appears ``m`` times (it has *multiplicity* ``m``), a spline of degree ``p`` is ``C^{p-m}`` there.

![Effect of the knot multiplicity on the continuity](figures/nurbs_multiplicity.svg)

*Cubic basis functions (``p = 3``) with the knot 2 inserted once, twice and three times. The functions that are non-zero at ``\xi = 2`` are blue.*

For a cubic spline, a simple knot gives ``C^2`` continuity. A double knot gives ``C^1``. A triple knot gives ``C^0``: the basis has a kink, one function reaches the value one, and the spline is interpolatory at that knot, exactly as at an element boundary in a standard finite element mesh. With multiplicity ``p + 1``, the spline is discontinuous; this is what happens at the ends of an open knot vector.

So the knot vector controls two things at once: *where* the elements are, and *how smoothly* they are connected. A spline with all interior knots repeated ``p`` times is an ordinary ``C^0`` finite element basis. That basis is called the Bernstein basis, and it is the link between splines and Ferrite (see [Bezier extraction](@ref)).

## Local knot vectors

Look again at the recursion formula: to compute ``N_{i,p}``, only the knots ``\xi_i, \ldots, \xi_{i+p+1}`` are used. Each basis function therefore has its own small **local knot vector** of ``p + 2`` knots, and that is all it needs. The rest of the global knot vector is irrelevant to it.

![A cubic basis function and its local knot vector](figures/nurbs_local_support.svg)

*The cubic function ``N_4`` of the knot vector ``\{0,0,0,0,1,2,3,4,5,5,5,5\}`` is defined by five knots only: ``\{0, 1, 2, 3, 4\}``. It is non-zero between the first and the last of them.*

For a B-spline, this is just a different way of looking at the same thing. It becomes the central idea for [T-splines](@ref), where each basis function gets its own local knot vector that is not cut out of one global knot vector.

## From basis functions to curves

A B-spline curve is a sum of control points ``\boldsymbol{P}_i`` multiplied by the basis functions,

```math
\boldsymbol{C}(\xi) = \sum_{i=1}^{n} N_{i,p}(\xi) \, \boldsymbol{P}_i .
```

As ``\xi`` runs through the knot vector, the point ``\boldsymbol{C}(\xi)`` traces out the curve. At every ``\xi``, the curve point is a weighted average of the control points, with the basis function values as weights.

![A quadratic B-spline curve with its control polygon](figures/nurbs_curve.svg)

*A quadratic B-spline curve (blue) built from the basis in the section [The open knot vector](@ref) and 8 control points. The grey line connecting the control points is the control polygon. The black dots are the images of the knots, i.e. the boundaries of the 5 elements.*

Compare the curve with the basis functions above:

- The curve starts at ``\boldsymbol{P}_1`` and ends at ``\boldsymbol{P}_8``, because the basis is interpolatory at the ends.
- In between, the curve does not pass through the control points. It follows the control polygon smoothly, like a rope pulled towards the control points. It is tangent to the control polygon at the ends.
- At ``\xi = 4``, the double knot makes the basis ``C^0``. The curve passes through ``\boldsymbol{P}_6`` and has a corner there.
- The curve always lies inside the *convex hull* of its control points, because the basis functions are non-negative and sum to one.
- The **elements** are the pieces of the curve between the black dots. The control points are not the nodes of the elements: there are 8 control points and 5 elements, and each element is influenced by ``p + 1 = 3`` control points.

The last point is the most important difference to a standard finite element mesh. In isogeometric analysis there are two meshes: the *control mesh* (control points) and the *physical mesh* (elements), and they are not the same thing.

### Local control

Since each basis function is non-zero on only ``p + 1`` knot spans, a control point only influences that part of the curve.

![Moving one control point changes the curve only locally](figures/nurbs_local_control.svg)

*Moving ``\boldsymbol{P}_4`` only changes the curve for ``1 \le \xi \le 4`` (orange), which is where ``N_4`` is non-zero. The rest of the curve (blue) stays where it was (grey: the original curve).*

For analysis this means sparse matrices, as for finite elements: two basis functions only couple if their supports overlap.

## Refining without changing the shape

To get a more accurate analysis, we need more basis functions. For splines, the basic tool is **knot insertion**: add a new knot to the knot vector, and compute new control points such that the curve is exactly the same as before.

![Knot insertion keeps the curve but adds control points and elements](figures/nurbs_knot_insertion.svg)

*Inserting the knots 1.5 and 2.5 adds two control points and two elements. The curve itself is unchanged. The new control polygon (orange) lies closer to the curve than the old one (light grey).*

Why is this possible? Because every old basis function can be written exactly as a combination of the new, narrower basis functions:

![A basis function written as a sum of refined basis functions](figures/nurbs_refinement_basis.svg)

*After inserting the knot 1.5, the old function ``N_4`` (dashed) is the sum of two new functions with the factors ½ and ⅚. The new space contains the old one, so nothing is lost by refining.*

In FerriteIGA, knot insertion is done with [`knotinsertion`](@ref):

```@example nurbs
using FerriteIGA

knots  = [0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5.0]
points = [Vec(0.0, 0.3), Vec(0.1, 1.9), Vec(1.6, 2.3), Vec(2.5, 0.2),
          Vec(3.7, 0.3), Vec(4.5, 2.2), Vec(5.3, 0.9), Vec(6.4, 1.3)]
curve = NURBSMesh((knots,), (2,), points)       # the curve in the figures above

refined = knotinsertion(curve, 1.5; dir = 1)
refined.knot_vectors[1]
```

```@example nurbs
length(curve.control_points), length(refined.control_points)
```

There are two more ways to refine a spline. *Order elevation* ([`orderelevation`](@ref)) raises the degree and keeps the continuity at the knots. Raising the degree first and inserting knots afterwards ([`smoothnesselevation`](@ref)) gives new elements that are joined with ``C^{p}`` continuity; this is called *k-refinement* and has no counterpart in standard finite elements. See the Refinement section of the [Splines](@ref) page.

## Rational splines

B-splines are piecewise polynomials, and polynomials cannot describe a circle exactly. Since circles, cylinders and spheres are everywhere in engineering, CAD uses a slightly more general function: each control point gets a positive **weight** ``w_i``, and the basis functions are divided by their weighted sum,

```math
R_{i}(\xi) = \frac{w_i \, N_{i,p}(\xi)}{\sum_{j=1}^{n} w_j \, N_{j,p}(\xi)},
\qquad
\boldsymbol{C}(\xi) = \sum_{i=1}^{n} R_{i}(\xi) \, \boldsymbol{P}_i .
```

The functions ``R_i`` are ratios of polynomials, hence *rational*. A curve built from them is a Non-Uniform Rational B-Spline, a **NURBS**. ("Non-uniform" only says that the knots need not be equally spaced.)

The weight of a control point tells how strongly it pulls the curve towards itself:

![The effect of a weight on a NURBS curve, and on the basis functions](figures/nurbs_circle.svg)

*Left: a quadratic curve with three control points and the weights ``(1, w_2, 1)``. A larger ``w_2`` pulls the curve towards ``\boldsymbol{P}_2``. With ``w_2 = 1`` the curve is an ordinary B-spline (a parabola), and with ``w_2 = \sqrt{2}/2`` it is exactly a quarter of a circle (the thick grey line underneath). Right: a large weight makes the basis function of that control point larger, at the expense of the others.*

NURBS keep all the good properties of B-splines: the rational functions are non-negative, have the same local support and the same continuity, and they sum to one. If all weights are equal, the denominator is constant and a NURBS is just a B-spline. Knot insertion works in the same way.

For the analysis, the weights are part of the geometry description. They are stored in the mesh (`mesh.weights`), and the unknown field uses the same rational basis functions as the geometry. This is the *isoparametric* idea of finite elements, and the origin of the name *isogeometric*.

## Surfaces: tensor products

A surface needs two parameters, ``\xi`` and ``\eta``. The basis functions of a spline surface are simply products of a B-spline in ``\xi`` and a B-spline in ``\eta``, each with its own knot vector and degree,

```math
N_{ij}(\xi, \eta) = N_{i}(\xi) \, M_{j}(\eta).
```

This is called a **tensor product**. With ``n`` functions in ``\xi`` and ``m`` functions in ``\eta`` there are ``n \cdot m`` basis functions, and equally many control points, arranged in a regular ``n \times m`` grid called the *control net*.

![A bivariate basis function as the product of two univariate functions](figures/nurbs_tensor_basis.svg)

*A biquadratic basis function for the knot vectors ``\{0,0,0,1,2,3,4,4,4\}`` in both directions. Its value at ``(\xi, \eta)`` is the product of the two curves at the top and at the right. It is non-zero on ``3 \times 3`` elements (orange frame).*

![The same basis function drawn as a surface](figures/nurbs_tensor_surface.svg)

*The same function as a surface over the parameter domain.*

The surface itself is defined as for curves, with a weight per control point for NURBS,

```math
\boldsymbol{S}(\xi, \eta) = \sum_{i=1}^{n} \sum_{j=1}^{m} R_{ij}(\xi, \eta) \, \boldsymbol{P}_{ij},
\qquad
R_{ij}(\xi, \eta) = \frac{w_{ij} \, N_{i}(\xi) \, M_{j}(\eta)}{\sum_{k,l} w_{kl} \, N_{k}(\xi) \, M_{l}(\eta)}.
```

A solid works in the same way with three parameters and products of three functions. One such block, with its knot vectors, degrees, control points and weights, is called a **patch**. In FerriteIGA a patch is a [`NURBSMesh`](@ref).

### Three ways to look at a patch

A patch can be drawn in three different spaces, and it helps to keep them apart.

![A NURBS patch in index space, parameter space and physical space](figures/nurbs_spaces.svg)

*A quadratic NURBS patch of a quarter ring with the knot vectors ``\{0,0,0,\tfrac13,\tfrac23,1,1,1\}`` and ``\{0,0,0,\tfrac12,1,1,1\}``. The 6 elements are blue in all three pictures.*

- In the **index space**, every knot is drawn as a line, with the same distance between all lines, whether or not the knots have the same value. Repeated knots are separate lines here. This space shows the *structure* of the patch: it is the natural place to count which knots a basis function uses. Only the cells between knots with different values are elements; the others have zero size.
- In the **parameter space**, the knots are drawn at their values. Repeated knots fall on top of each other, and what remains are the elements, as rectangles.
- The **physical space** is where the geometry lives. The elements are the images of the parameter rectangles under the map ``\boldsymbol{S}(\xi, \eta)``. The control net is also drawn here. Notice once more that the control points are not the corners of the elements, and that most of them do not lie on the surface.

The index space may look like an unnecessary detour here, but it is the space in which [T-splines](@ref) are defined.

## A patch in FerriteIGA

The quarter ring in the figure can be built directly from its definition. The coarsest description has a single element: three control points across the thickness, three around the arc, and the weight ``\sqrt{2}/2`` on the middle row, as for the quarter circle above. The control points are numbered with the first parameter direction running fastest. The arc is traversed clockwise, so that the two parameter directions form a right-handed system and the Jacobian of the map is positive.

```@example nurbs
s = sqrt(2) / 2
knots = [0, 0, 0, 1, 1, 1.0]
control_points = Vec{2,Float64}[]
weights = Float64[]
for r in (1.0, 1.5, 2.0)                                # second direction (η): the radius
    # first direction (ξ): along the arc, from the y-axis to the x-axis
    for (x, w) in ((Vec(0.0, 1.0), 1.0), (Vec(1.0, 1.0), s), (Vec(1.0, 0.0), 1.0))
        push!(control_points, r * x)
        push!(weights, w)
    end
end
ring = NURBSMesh((knots, copy(knots)), (2, 2), control_points, weights)
```

Knot insertion gives the mesh in the figure, without changing the geometry:

```@example nurbs
ring = knotinsertion(ring, 1/3; dir = 1)
ring = knotinsertion(ring, 2/3; dir = 1)
ring = knotinsertion(ring, 1/2; dir = 2)
ring.knot_vectors
```

For the analysis, the patch is converted to a [`BezierGrid`](@ref), which is a Ferrite grid with one cell per element. From here on, everything works as in Ferrite. As a check, we integrate the area of the ring, and verify that the basis functions sum to one at every quadrature point:

```@example nurbs
grid = BezierGrid(ring)
ip = IGAInterpolation{RefQuadrilateral, 2}()
cv = BezierCellValues(QuadratureRule{RefQuadrilateral}(4), ip)

area = 0.0
partition_of_unity = true
for cellid in 1:getncells(grid)
    reinit!(cv, getcoordinates(grid, cellid))
    for qp in 1:getnquadpoints(cv)
        global area += getdetJdV(cv, qp)
        global partition_of_unity &= sum(shape_value(cv, qp, i) for i in 1:getnbasefunctions(cv)) ≈ 1
    end
end
area, 3π/4, partition_of_unity
```

The computed area agrees with the exact area ``\tfrac{\pi}{4}(2^2 - 1^2) = \tfrac{3\pi}{4}`` of a quarter ring up to the accuracy of the quadrature, with only 6 elements. The geometry is exact; no finer mesh would describe the ring any better.

### What happens inside `BezierGrid`

Ferrite works element by element, with the same set of shape functions on every element. B-splines do not fit this pattern directly, since each basis function stretches over several elements. The solution is **Bézier extraction**: on each element, the ``p + 1`` B-splines that are non-zero there are written as combinations of the ``p + 1`` *Bernstein polynomials* of that element.

![B-splines and the Bernstein polynomials of each element](figures/nurbs_bezier.svg)

*Top: quadratic B-splines on three elements. Bottom: the Bernstein polynomials, which are the same on every element. On the shaded element, each of the B-splines ``N_2``, ``N_3`` and ``N_4`` is a combination of ``B_1``, ``B_2`` and ``B_3``. The coefficients form the extraction operator of that element.*

`IGAInterpolation` is the Bernstein basis, and `reinit!` multiplies it with the extraction operator of the cell (and the weights) to get the NURBS functions. The details are on the [Bezier extraction](@ref) page.

## The limitation of tensor products

The tensor product makes NURBS simple, but it has a price. A new knot in the ``\xi``-direction is a new knot for *all* basis functions, so it creates a line of new elements and control points across the whole patch.

![Tensor product refinement is global](figures/nurbs_global_refinement.svg)

*Left: to split the blue element of a NURBS patch, the new knot lines (orange) must run through the whole patch, which refines elements that did not need it. Right: the local refinement we would like to have. The new lines end in the middle of the mesh, in so-called T-junctions.*

In an analysis, we typically want small elements in a small region (a stress concentration, a boundary layer) and large elements elsewhere. With NURBS, every local refinement turns into rows and columns of unwanted elements. Removing this restriction is the idea behind T-splines, which are the subject of the next tutorial: [T-splines](@ref).

## Summary

- A B-spline basis is defined by a degree ``p`` and a knot vector. Knot spans are elements; ``n + p + 1`` knots give ``n`` basis functions.
- Each basis function is non-zero on ``p + 1`` knot spans and is defined by ``p + 2`` consecutive knots, its local knot vector. On each element, ``p + 1`` functions are non-zero.
- A knot with multiplicity ``m`` gives ``C^{p-m}`` continuity. Open knot vectors repeat the end knots ``p + 1`` times, which makes the basis interpolatory at the ends.
- Curves, surfaces and solids are sums of control points times basis functions. The control points are coefficients, not nodal values.
- Knot insertion adds basis functions without changing the geometry.
- NURBS add a weight per control point, which makes circles and other conic sections exact.
- Surfaces and solids are tensor products, which makes refinement global.

## Further reading

- L. Piegl and W. Tiller. *The NURBS Book*, 2nd edition. Springer, 1997.
- J. A. Cottrell, T. J. R. Hughes, and Y. Bazilevs. *Isogeometric Analysis. Toward Integration of CAD and FEA*. Wiley, 2009. [doi:10.1002/9780470749081](https://doi.org/10.1002/9780470749081)
- T. J. R. Hughes, J. A. Cottrell, and Y. Bazilevs. Isogeometric analysis: CAD, finite elements, NURBS, exact geometry and mesh refinement. *Comput. Methods Appl. Mech. Engrg.*, 194:4135–4195, 2005. [doi:10.1016/j.cma.2004.10.008](https://doi.org/10.1016/j.cma.2004.10.008)
- Y. Bazilevs et al. Isogeometric analysis using T-splines. *Comput. Methods Appl. Mech. Engrg.*, 199:229–263, 2010. [doi:10.1016/j.cma.2009.02.036](https://doi.org/10.1016/j.cma.2009.02.036) (Section 2 is a compact review of B-splines and NURBS.)
