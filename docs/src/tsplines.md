# T-splines

This page explains T-splines by starting from what is already known about NURBS, and then shows how they are used in FerriteIGA. It assumes that you have read the [B-splines and NURBS](@ref) page, in particular the sections on [Local knot vectors](@ref) and the index space ([Three ways to look at a patch](@ref)).

The first part of the page is an introduction. The sections from [Grids with T-spline cells](@ref) onwards are reference material. All figures are computed with FerriteIGA (see `docs/generate_figures.jl`).

## Why T-splines?

A NURBS patch is a tensor product: every knot line runs through the whole patch. Making one element smaller therefore makes a whole row and a whole column of elements smaller.

![Tensor product refinement is global](figures/nurbs_global_refinement.svg)

*Left: refining one element of a NURBS patch. Right: the local refinement we would like to have.*

T-splines (Sederberg et al., 2003) remove this restriction. A line of the mesh is allowed to end in the middle of the mesh. The point where it ends is called a **T-junction**, because the mesh looks like the letter T there. Everything else stays as close to NURBS as possible: the basis functions are still products of B-splines, there is still one control point per basis function, and a T-spline without T-junctions *is* a NURBS.

The question that the rest of this introduction answers is: if the knot lines no longer run through the whole patch, what are the basis functions?

## A NURBS patch seen as a T-spline

The answer builds on an observation from the NURBS page: a B-spline basis function of degree ``p`` does not need the whole knot vector. It only needs its own *local knot vector* of ``p + 2`` consecutive knots. So instead of saying "here is a global knot vector, and function number ``i`` uses the knots ``i`` to ``i + p + 1``", we can give each function a position in the mesh, called its **anchor**, and a rule for collecting the knots around the anchor. For a NURBS this gives nothing new. For a T-spline, it is the definition.

### One direction

Start with one parameter direction. In the index space, each knot is drawn as a line, with the same distance between all lines.

![Knot vectors, index lines and anchors in one direction](figures/tspline_anchors_1d.svg)

*Top: the basis functions of an open knot vector for ``p = 3`` and ``p = 2``. Bottom: the same functions in index space. Each function has an anchor (coloured dot). For odd degrees, the anchors sit on the lines; for even degrees, they sit between the lines. The green bracket shows the lines that* ``N_3`` *uses.*

The rule for the local knot vector is:

- **Odd degree** (``p = 1, 3, \ldots``): the anchor sits on a line. Take the knot of that line, and the knots of the ``(p+1)/2`` nearest lines on each side. For ``p = 3`` that is ``1 + 2 + 2 = 5 = p + 2`` knots.
- **Even degree** (``p = 2, 4, \ldots``): the anchor sits in the middle between two lines. Take the knots of the ``p/2 + 1`` nearest lines on each side. For ``p = 2`` that is ``2 + 2 = 4 = p + 2`` knots.
- **At the boundary**: if there are not enough lines on one side, repeat the knot of the last line.

Check the rule for the cubic function ``N_3`` in the figure. Its anchor is on line 3. Two lines to the left and two to the right gives the lines 1 to 5, with the knots ``\{0, 0, 1, 2, 3\}``. This is exactly the local knot vector of the third function of ``\{0,0,0,0,1,2,3,3,3,3\}``. For ``N_1``, with its anchor on line 1, there are no lines to the left, so the knot 0 is repeated: ``\{0, 0, 0, 0, 1\}``.

!!! note "Index lines and open knot vectors"
    In FerriteIGA, a T-spline is described by the knot value of each index line. Because the last knot is repeated automatically at the boundary, an open boundary needs fewer repeated knots than an open knot vector: the boundary value is given `p÷2 + 1` times instead of `p + 1` times. In the figure, the index lines have the knots `[0, 0, 1, 2, 3, 3]` both for the cubic knot vector `[0, 0, 0, 0, 1, 2, 3, 3, 3, 3]` and for the quadratic knot vector `[0, 0, 0, 1, 2, 3, 3, 3]`. There is one anchor per basis function.

### Two directions

In two dimensions the index space is a grid of lines. The anchors are the **vertices** of the grid for odd degrees, and the **centres of the cells** for even degrees. The local knot vectors are found by applying the rule above in each direction: start at the anchor, walk to the left and to the right, and collect the vertical lines that you cross; then walk down and up, and collect the horizontal lines.

![A tensor product mesh in index space with anchors and the lines used by one anchor](figures/tspline_index_nurbs.svg)

*A NURBS patch in index space, with one anchor (grey dot) per basis function. For the marked anchor, the thick lines are the paths walked from the anchor, and the ticks are the lines that are collected. The shaded rectangle between the outermost ticks is the support of the basis function.*

The basis function of an anchor is the product of the two B-splines defined by its two local knot vectors,

```math
N_A(\xi, \eta) = N[\Xi_A](\xi) \, N[\mathrm{H}_A](\eta),
```

where ``\Xi_A`` and ``\mathrm{H}_A`` are the local knot vectors in the two directions. For a tensor product mesh, this gives exactly the NURBS basis. In FerriteIGA, a NURBS patch can be converted with `TSplineMesh(nurbsmesh)`:

```@example tsplines
using FerriteIGA

nurbs = generate_nurbs_patch(:rectangle, (4, 4), (3, 3); size = (4.0, 4.0))
TSplineMesh(nurbs)
```

## The T-mesh

Now remove the tensor product restriction. A **T-mesh** is a division of the index space into rectangular cells. The cells must fit together without gaps and overlaps, but a cell may be several index units wide or high, so a line can end where a larger cell begins. The corners of the cells are the vertices of the T-mesh.

![A T-mesh in index space and in parameter space](figures/tspline_tmesh.svg)

*A T-mesh with two T-junctions (orange circles). The vertical line* ``i = 5`` *ends at* ``j = 4``*, and the horizontal line* ``j = 7`` *ends at* ``i = 7``*. Left: the index space, where the two lines with the knot value 0 (and 6) at the boundary are drawn separately. The grey dots are the anchors for an odd degree. Right: the same T-mesh in parameter space.*

In FerriteIGA, a T-mesh is given by the knot values of the index lines in each direction, and a list of cells `(i0, i1, j0, j1)`, where `i0:i1` and `j0:j1` are the index ranges of the cell. A tensor product mesh consists of unit cells:

```@example tsplines
knots = ([0.0, 0, 1, 2, 3, 4, 4], [0.0, 0, 1, 2, 3, 4, 4])     # knot value of each index line
cells = [(i, i + 1, j, j + 1) for j in 1:6 for i in 1:6]        # 6 × 6 unit cells
TSplineMesh(knots, 3, cells)
```

To let the vertical line ``i = 4`` end at ``j = 4``, we replace the cells on both sides of it with wider cells above ``j = 4``:

```@example tsplines
cells = NTuple{4,Int}[]
for j in 1:6, i in 1:6
    if j >= 4 && i == 3
        push!(cells, (3, 5, j, j + 1))      # a wide cell: there is no line i = 4 here
    elseif j >= 4 && i == 4
        continue                            # covered by the wide cell
    else
        push!(cells, (i, i + 1, j, j + 1))
    end
end
tmesh = TSplineMesh(knots, 3, cells)
```

The T-mesh has three vertices fewer than the tensor product mesh, and therefore three control points fewer. (The number of Bézier elements is explained [below](@ref "Bézier elements").)

A few rules apply to the T-mesh ([`TSplineMesh`](@ref) checks them):

- The cells must tile the index domain.
- The index lines with a repeated knot value at the boundary must be complete lines. T-junctions belong in the interior.
- If two vertices lie exactly opposite to each other on the two sides of a cell, they must be connected by a line. (Otherwise the mesh would have two T-junctions facing each other.)

## Basis functions of a T-spline

The rule for the local knot vectors is the same as for the NURBS patch above, with one clarification: when walking from the anchor, a line only counts **if it is actually there** on the path. Where a line has ended in a T-junction, there is nothing to cross, and the walk continues to the next line.

![Local index vectors of three anchors in a cubic T-spline](figures/tspline_rays_odd.svg)

*Three anchors of the cubic T-spline. Left: far from the T-junctions, the function uses five consecutive lines in each direction, as in a NURBS. Middle: the path to the right passes the end of the line* ``i = 5``*, which is not counted, so the function reaches further, to* ``i = 6``*. The same happens upwards with the line* ``j = 7``*. Right: the anchor at a T-junction. The line ends at the anchor, but the walk upwards continues through the cell above and collects the lines* ``j = 5`` *and* ``j = 6``*.*

For even degrees, the anchors are the centres of the cells, and the walk starts there:

![Local index vectors of three anchors in a quadratic T-spline](figures/tspline_rays_even.svg)

*The same T-mesh with ``p = 2``. Each cell has an anchor in its centre, also the large cells next to the T-junctions.*

The result of the walk is a list of line indices for each direction, the *local index vectors*. Replacing each index with its knot value gives the local knot vectors, and the basis function is the product of the two B-splines, exactly as for NURBS. The difference is that the supports of the functions are no longer arranged in a regular pattern:

![Two basis functions of a cubic T-spline](figures/tspline_function.svg)

*The basis functions of the middle and the right anchor of the cubic figure above, in parameter space. The frame marks the support. The left function skips the knot 2.5, since the line ends before it reaches the function. The right function, anchored at the T-junction, is narrower in* ``\xi`` *than its neighbours above.*

In the code, the local index and knot vectors are stored for each basis function. For the T-mesh that we built above, the function anchored at the vertex ``(3, 5)`` skips the line ``i = 4``:

```@example tsplines
A = findfirst(==((2 * 3, 2 * 5)), tmesh.anchors)     # the anchors are stored with doubled index coordinates
tmesh.local_indices[A]
```

```@example tsplines
tmesh.local_knots[A]
```

### Control points and the surface

There is one control point ``\boldsymbol{P}_A`` and one weight ``w_A`` for each anchor, and the T-spline surface is defined as for NURBS,

```math
\boldsymbol{S}(\xi, \eta) = \sum_{A} R_A(\xi, \eta) \, \boldsymbol{P}_A,
\qquad
R_A(\xi, \eta) = \frac{w_A \, N_A(\xi, \eta)}{\sum_{B} w_B \, N_B(\xi, \eta)} .
```

The control points are passed to [`TSplineMesh`](@ref) with the keyword `control_points`, numbered like the anchors: row by row from the bottom, from left to right in each row. If no control points are given, the *Greville points* of the basis functions are used ([`greville_points`](@ref)). The Greville point of a function is the mean of the interior knots of its local knot vectors. With these control points (and an analysis-suitable mesh, see below), the map ``\boldsymbol{S}`` is the identity, so the physical domain equals the parameter domain. This is the natural choice for rectangular domains.

## Bézier elements

For NURBS, the elements are the knot spans. For T-splines, there is a subtlety. A basis function is a piecewise polynomial, with the pieces separated by *its own* knot lines, and these lines cross the whole support of the function. So near a T-junction, the basis functions have knot lines where the T-mesh has no line.

The **Bézier elements** are the regions where all basis functions are polynomials. They are obtained by drawing the knot lines of all basis functions on top of the T-mesh:

![The Bézier elements of a T-spline](figures/tspline_elements.svg)

*The T-mesh (black) and the additional element boundaries (orange, dashed) in parameter space. The line that ended in a T-junction continues as an element boundary, further for* ``p = 3`` *than for* ``p = 2``*. The basis functions are smooth across the black lines to the same degree as across the orange ones; the orange lines are only needed to integrate correctly.*

This is why the T-mesh `tmesh` above has 16 Bézier elements, just like the tensor product mesh, although it has only 14 cells with a non-zero area.

`BezierGrid(mesh)` creates a grid with one cell per Bézier element. As for NURBS, each cell stores a Bézier extraction operator, which expresses the T-spline basis functions that are non-zero on the element in the Bernstein polynomials of the element (Scott et al., 2011). After that, the analysis works as for NURBS: see [Solving a problem](@ref) below.

## Analysis-suitable T-splines

So far, any T-mesh gives a set of basis functions. But not every T-mesh gives a *good* set of basis functions. The functions of a general T-spline may fail to sum to one, they can fail to reproduce polynomials (which destroys the accuracy of the analysis), and they can even be linearly dependent. The figure shows two T-meshes that differ only in the length of two lines:

![The sum of all basis functions for an analysis-suitable and a non analysis-suitable T-mesh](figures/tspline_pou.svg)

*The sum of all cubic basis functions minus one. Left: the sum is one everywhere, and there are 16 functions on every element, as for a bicubic NURBS. Right: both lines are one element longer. Around the two T-junctions, the functions sum to more than one, and some elements have 17 functions.*

The T-meshes that are free of these problems are called **analysis-suitable** (Li et al., 2012; Beirão da Veiga et al., 2013). Whether a T-mesh is analysis-suitable can be read off the mesh, with a simple geometric rule.

### T-junction extensions

Each T-junction gets an **extension**: a line segment through the T-junction, in the direction of the line that ended there. It has two parts.

- The **face extension** continues the line into the cell where it is missing.
- The **edge extension** goes the other way, back along the line.

The length of the extension is measured in *bays*: one bay is the distance to the next crossing line. For odd degrees, the face extension is ``(p+1)/2`` bays long and the edge extension ``(p-1)/2`` bays. For even degrees, both are ``p/2`` bays long. In total, the extension always spans ``p`` bays.

![Face and edge extensions of a T-junction](figures/tspline_extensions.svg)

*The extension of a T-junction for* ``p = 3`` *(face extension 2 bays, edge extension 1 bay) and* ``p = 2`` *(1 bay each).*

The extensions are not an arbitrary construction. Compare with the figure in [Bézier elements](@ref): the face extensions are exactly the additional element boundaries. They show how far the influence of a T-junction reaches.

### The rule

A T-junction is *horizontal* if the line that ends there is horizontal, and *vertical* otherwise.

> A T-mesh is **analysis-suitable** if no extension of a horizontal T-junction intersects an extension of a vertical T-junction.

![T-junction extensions of an analysis-suitable and a non analysis-suitable T-mesh](figures/tspline_as.svg)

*The two T-meshes of the figure above in index space, with the extensions of the horizontal T-junction (blue) and of the vertical T-junction (orange) for* ``p = 3``*. Left: the extensions do not meet. Right: with longer lines, the T-junctions come too close to each other and the extensions cross.*

The check is available as [`is_analysis_suitable`](@ref), and [`tjunction_extensions`](@ref) returns the extensions. The T-mesh that we built by hand has a single T-junction, so nothing can intersect:

```@example tsplines
is_analysis_suitable(tmesh), tjunction_extensions(tmesh)
```

Now we also let the horizontal line ``j = 5`` end at ``i = 3``, by replacing four cells to the left of it with two tall cells. The horizontal extension then crosses the vertical one:

```@example tsplines
bad_cells = filter(c -> !(c[1] <= 2 && c[3] in (4, 5)), cells)
append!(bad_cells, [(1, 2, 4, 6), (2, 3, 4, 6)])
bad = TSplineMesh(knots, 3, bad_cells)
is_analysis_suitable(bad), tjunction_extensions(bad)
```

```@example tsplines
extrema(length, bad.IEN)        # number of basis functions on the Bézier elements, 16 for a bicubic NURBS
```

Since the extensions get longer with the degree, a T-mesh that is analysis-suitable for one degree is also analysis-suitable for all lower degrees, but not necessarily for higher ones. The mesh to the right in the figure is analysis-suitable for ``p = 2``.

### What analysis-suitable T-splines guarantee

For an analysis-suitable T-spline of any degree:

- The basis functions are **linearly independent**, for any choice of knot values.
- They form a **partition of unity** and reproduce all polynomials of degree ``p``. A mesh with its control points at the Greville points therefore describes its parameter domain exactly, and patch tests are passed.
- Exactly ``(p+1)^2`` basis functions are non-zero on each Bézier element, as for a NURBS. All cells of the `BezierGrid` then have the same type.
- The approximation is optimal: the error in the ``L^2`` norm decreases as ``h^{p+1}``, as for NURBS.
- If the *extended* T-mesh (the T-mesh together with all T-junction extensions) of one T-spline is contained in the extended T-mesh of another, then the spline spaces are nested (Li and Scott, 2014). This is the basis for local refinement that keeps the geometry unchanged.

The mathematical reason is that analysis-suitable T-meshes are exactly the *dual-compatible* ones (Beirão da Veiga et al., 2013): any two basis functions share, in at least one direction, a common knot vector from which both local knot vectors are taken. In that direction they behave like two B-splines of the same knot vector, and this is enough to carry over the proofs from B-splines. [`is_dual_compatible`](@ref) checks this property directly on the local index vectors.

### Which T-meshes are analysis-suitable?

- **Tensor product meshes**, since they have no T-junctions.
- **Meshes where all T-junctions have the same direction.** If only vertical lines end, there are no horizontal extensions, and nothing can intersect. This is the case for layered meshes that get finer towards one side, which is what `generate_tspline_patch(:graded_rectangle, ...)` creates.
- **Meshes where horizontal and vertical T-junctions are far enough apart.** When a region is refined in both directions, some of the new lines usually have to be made longer than strictly necessary, so that the extensions of the T-junctions do not meet.

A typical mesh that is *not* analysis-suitable has a horizontal and a vertical T-junction close to each other. The simplest example is a single cell that is split into four: the four new T-junctions have extensions that cross in the middle of the cell.

For a T-mesh that is not analysis-suitable, none of the guarantees above hold. Elements can have more or fewer than ``(p+1)^2`` basis functions, and linear independence depends on the knot values ([`is_linearly_independent`](@ref) checks it numerically). FerriteIGA is designed for analysis-suitable T-meshes; see [Limitations](@ref) for what works otherwise.

## Graded meshes

Writing down the cells of a T-mesh by hand quickly becomes tedious. [`generate_tspline_patch`](@ref) generates the most common locally refined mesh: a rectangle that is divided into horizontal layers, each with its own number of elements.

![A graded T-spline mesh](figures/tspline_graded.svg)

*A cubic T-spline of the rectangle* ``[0, 2] \times [0, 1]`` *with three layers of 2 × 1, 4 × 2 and 8 × 4 elements, from bottom to top. Black: the T-mesh. Orange: additional boundaries of the Bézier elements. Blue dots: the control points (Greville points).*

The horizontal lines span the whole width, and the vertical lines of a layer end at the layer below. All T-junctions are vertical, so the mesh is analysis-suitable for any numbers of elements in the layers.

```@example tsplines
mesh = generate_tspline_patch(:graded_rectangle, [2, 4, 8], [1, 2, 4], 3; size = (2.0, 1.0))
```

```@example tsplines
is_analysis_suitable(mesh)
```

## Solving a problem

From here on, a T-spline is used like a NURBS patch. `BezierGrid(mesh)` creates the grid, `IGAInterpolation` is the Bernstein basis of the elements, and `reinit!` applies the Bézier extraction:

```@example tsplines
grid = BezierGrid(mesh)
ip = IGAInterpolation{RefQuadrilateral, 3}()
cv = BezierCellValues(QuadratureRule{RefQuadrilateral}(4), ip)

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, union(getfacetset(grid, "left"), getfacetset(grid, "bottom")), x -> 0.0))
close!(ch)

ndofs(dh)
```

The assembly loop is the usual one:

```julia
for cellid in 1:getncells(grid)
    coords = getcoordinates(grid, cellid)
    reinit!(cv, coords)
    n = getnbasefunctions(cv)
    dofs = celldofs(dh, cellid)
    # ... compute the element matrix with shape_value(cv, qp, i) and shape_gradient(cv, qp, i),
    #     and assemble as usual
end
```

A complete example is [Transient heat equation on a graded T-spline mesh](@ref).

## Grids with T-spline cells

`BezierGrid(mesh::TSplineMesh)` creates one [`TSplineCell`](@ref) per Bézier element. Its nodes are the control points whose basis functions are non-zero on the element. For analysis-suitable meshes, all cells have ``(p+1)^2`` nodes and the same type. If the number of nodes differs between the elements (which only happens for meshes that are not analysis-suitable), the grid has several cell types. In that case:

- `add!(dh, :u, ip)` creates one `SubDofHandler` per cell type and uses a [`TSplineInterpolation`](@ref) with the right number of functions for each.
- Use `ndofs_per_cell(dh, cellid)` or `getnbasefunctions(cv)` (after `reinit!`) instead of `ndofs_per_cell(dh)`. Use `CellIterator(sdh)` for each SubDofHandler instead of `CellIterator(dh)`.

The grid contains the facetsets and nodesets `"left"`, `"right"`, `"bottom"` and `"top"` of the parametric domain. Facet integrals (e.g. Neumann conditions) work with `BezierFacetValues` on these facetsets, as for NURBS.

## Periodic T-splines

The T-mesh may extend outside of the parametric domain (the `domain` keyword of [`TSplineMesh`](@ref)). Bézier elements are then only created inside the domain, and only the basis functions that are non-zero in the domain are kept. This is used for periodic splines: with `generate_tspline_patch(:graded_rectangle, ...; periodic = true)`, the knot lines of each layer are extended periodically, by `p` element widths, outside of the left and right sides. The basis functions are not interpolatory at these sides, and the `p` columns of basis functions that are non-zero on the left side are translations of those on the right side. Tying them with

```julia
add!(ch, PeriodicIGADirichlet(:u, getfacetset(grid, "left"), getfacetset(grid, "right")))
```

gives a periodic solution with ``C^{p-1}`` continuity across the sides, whereas tying the control points of a mesh with open knot vectors only gives ``C^0`` continuity.

## 3D T-splines

T-splines in 3D are defined in the same way: the T-mesh consists of boxes `(i0, i1, j0, j1, k0, k1)` in index space, the anchors are the vertices of the T-mesh (odd degrees) or the centres of the boxes (even degrees), and the local knot vectors are found by walking through the T-mesh in each of the three directions. Each basis function is a product of three univariate B-splines, and each Bézier element (a `TSplineCell{RefHexahedron}`) has a Bernstein basis with ``(p+1)^3`` functions.

The vertices of a 3D T-mesh are the points where planes of all three directions meet. These are the corners of the boxes, but also points where two edges of the T-mesh cross without being a corner of a box, for example where a plane `x = const` that ends at a layer interface from below meets a plane `y = const` that ends there from above.

`generate_tspline_patch(:graded_box, nels_x, nels_y, nels_z, p)` generates a box which is layered in the z-direction. Layer `k` has `nels_x[k] × nels_y[k] × nels_z[k]` elements; the planes `z = const` span the whole box, while the planes `x = const` and `y = const` end at the layer interfaces. The layers do not have to refine each other. With `periodic = true` (or a tuple `(px, py)`), the knot vectors are periodic in x (and/or y), as for [Periodic T-splines](@ref). The grid contains the sets `"left"`, `"right"` (x), `"front"`, `"back"` (y) and `"bottom"`, `"top"` (z).

```julia
# 2×2×1 elements at the bottom, 4×4×2 elements at the top, periodic in x
mesh = generate_tspline_patch(:graded_box, [2, 4], [2, 4], [1, 2], 3; periodic = (true, false))
is_analysis_suitable(mesh) # true
grid = BezierGrid(mesh)
```

There is no established definition of T-junction extensions in 3D. Instead, *dual-compatibility* is used as the criterion: [`is_dual_compatible`](@ref) checks that for every pair of basis functions with overlapping index ranges, there is a direction in which their local index vectors are taken from a common index vector and differ. Then the products of the univariate B-spline dual functionals form a dual basis, which gives linear independence and at most ``(p+1)^3`` functions per element. In 2D, dual-compatibility is equivalent to analysis-suitability (Beirão da Veiga et al., 2013).

## Boundary conditions

The degrees of freedom of a spline are the coefficients of the basis functions, not nodal values. The following constraints are available for `BezierGrid`s (NURBS and T-splines):

- `Dirichlet` on a facetset prescribes the control points whose basis functions are non-zero on the facets, with the value of `f` evaluated at the control points. This is exact for constant (and, with control points at the Greville points, linear) boundary data, but otherwise only an approximation. It supports time dependent data via `update!`.
- [`ProjectedIGADirichlet`](@ref) computes the boundary values by an L2 projection onto the spline functions on the boundary, which is exact whenever the data can be represented by the spline, and generally more accurate. The projection is computed once, when it is added to the `ConstraintHandler`, so the data can not depend on time.
- [`PeriodicIGADirichlet`](@ref) ties the control points on two identical sides of the domain (up to a translation) with affine constraints. See [Periodic T-splines](@ref) for smooth periodicity.

Ferrite's `ProjectedDirichlet` is not supported for spline fields.

## Limitations

- The degree is the same in all directions. Parametric dimensions 2 and 3 are supported.
- There is no local refinement algorithm that computes the control points of a refined T-spline from a coarser one. A T-spline is either generated with [`generate_tspline_patch`](@ref), converted from a NURBS patch, or given with its control points.
- `TSplineMesh(nurbsmesh)` requires open knot vectors, and throws an error otherwise.
- The facets of a `TSplineCell` can not be computed from its nodes, so functions that rely on the grid topology, such as `addfacetset!` and `ExclusiveTopology`, do not work; use the predefined sets instead. `L2Projector` and `apply_analytical_iga!` are not supported for T-spline grids.
- T-meshes that are not analysis-suitable are only partly supported. `BezierCellValues` require exactly ``(p+1)^{dim}`` basis functions on an element, and throw an error on other elements. A `BezierGrid` can be created if all elements have at least ``(p+1)^{dim}`` basis functions.
- The index vectors of the basis functions are computed by walking along the lines through the anchors (Bazilevs et al., 2010). For even degrees, Beirão da Veiga et al. (2013) instead use the lines that span the whole cell of the anchor; the two definitions coincide for analysis-suitable T-meshes.

## References

- T. W. Sederberg, J. Zheng, A. Bakenov, and A. Nasri. T-splines and T-NURCCs. *ACM Trans. Graph.*, 22:477–484, 2003.
- Y. Bazilevs, V. M. Calo, J. A. Cottrell, J. A. Evans, T. J. R. Hughes, S. Lipton, M. A. Scott, and T. W. Sederberg. Isogeometric analysis using T-splines. *Comput. Methods Appl. Mech. Engrg.*, 199:229–263, 2010. [doi:10.1016/j.cma.2009.02.036](https://doi.org/10.1016/j.cma.2009.02.036)
- M. A. Scott, M. J. Borden, C. V. Verhoosel, T. W. Sederberg, and T. J. R. Hughes. Isogeometric finite element data structures based on Bézier extraction of T-splines. *Int. J. Numer. Meth. Engng.*, 88:126–156, 2011. [doi:10.1002/nme.3167](https://doi.org/10.1002/nme.3167)
- X. Li, J. Zheng, T. W. Sederberg, T. J. R. Hughes, and M. A. Scott. On linear independence of T-spline blending functions. *Comput. Aided Geom. Design*, 29:63–76, 2012.
- L. Beirão da Veiga, A. Buffa, G. Sangalli, and R. Vázquez. Analysis-suitable T-splines of arbitrary degree: definition, linear independence and approximation properties. *Math. Models Methods Appl. Sci.*, 23:1979–2003, 2013. [doi:10.1142/S0218202513500231](https://doi.org/10.1142/S0218202513500231)
- X. Li and M. A. Scott. Analysis-suitable T-splines: characterization, refineability, and approximation. *Math. Models Methods Appl. Sci.*, 24:1141–1164, 2014. [arXiv:1211.5669](https://arxiv.org/abs/1211.5669)

## API

```@docs
TSplineMesh
generate_tspline_patch
greville_points
is_analysis_suitable
tjunction_extensions
is_dual_compatible
is_linearly_independent
TSplineCell
TSplineInterpolation
Ferrite.add!(::DofHandler{<:Any, <:BezierGrid}, ::Symbol, ::Interpolation)
Ferrite.add!(::ConstraintHandler{<:DofHandler{<:Any, <:BezierGrid}}, ::Dirichlet)
```
