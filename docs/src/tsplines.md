# T-splines

A NURBS patch is a tensor product, so refining one element inserts a knot line through the whole patch. T-splines (Sederberg, Zheng, Bakenov and Nasri, 2003) remove this restriction: lines in the control mesh may end in *T-junctions*, which allows local refinement, for example a mesh that is fine at the top and coarse at the bottom.

As for NURBS, a [`TSplineMesh`](@ref) describes the spline, and `BezierGrid(mesh)` converts it to a grid of Bézier elements with extraction operators ([Scott et al., 2011](https://doi.org/10.1002/nme.3167)). Assembly then works as for NURBS.

```julia
using FerriteIGA

# Coarse (2 columns) at the bottom, fine (8 columns) at the top
mesh = generate_tspline_patch(:graded_rectangle, [2, 4, 8], [1, 2, 4], 3)
is_analysis_suitable(mesh) # true

grid = BezierGrid(mesh)
ip = IGAInterpolation{RefQuadrilateral, 3}()
cv = BezierCellValues(QuadratureRule{RefQuadrilateral}(5), ip)

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, union(getfacetset(grid, "left"), getfacetset(grid, "bottom")), x -> 0.0))
close!(ch)

for cellid in 1:getncells(grid)
    coords = getcoordinates(grid, cellid)
    reinit!(cv, coords)
    n = getnbasefunctions(cv)     # may vary between cells, see below
    dofs = celldofs(dh, cellid)   # length n
    # ... assemble as usual
end
```

## The T-mesh

The T-mesh is defined in *index space*. The vertical lines have indices `i = 1:n1` and knot values `knots[1][i]`, and the horizontal lines have indices `j = 1:n2` and knot values `knots[2][j]`. The T-mesh is given as a list of cells, which are rectangles `(i0, i1, j0, j1)` that tile the index domain `[1,n1]×[1,n2]` (boxes `(i0, i1, j0, j1, k0, k1)` in 3D, see [3D T-splines](@ref)). A face boundary may contain vertices of neighbouring faces; these vertices are the T-junctions.

Consecutive knot values may be equal, which gives zero knot intervals. Repeat the boundary value `p÷2 + 1` times to get an open (interpolatory) boundary, as for open knot vectors in NURBS. For example, `[0, 0, 1, 2, 3, 3]` is the cubic counterpart of `[0, 0, 0, 0, 1, 2, 3, 3, 3, 3]`. The lines with zero knot interval at the boundary must span the whole domain.

A basis function is associated with each *anchor*. For odd degrees the anchors are the vertices of the T-mesh; for even degrees they are the centres of the faces. The local knot vectors of a basis function are found by marching from its anchor in each direction until `(p+1)÷2` (odd `p`) or `p÷2 + 1` (even `p`) orthogonal lines have been crossed. If the boundary is reached first, the last knot is repeated ([Bazilevs et al., 2010](https://doi.org/10.1016/j.cma.2009.02.036)). The control points are numbered as the anchors, with the x-index running fastest (row by row from bottom to top in 2D).

The Bézier elements are bounded by all knot lines of all basis functions (the *elemental T-mesh*). This means that an element can be smaller than a face of the T-mesh.

## Analysis-suitable T-splines

General T-splines can have linearly dependent basis functions, and they need not be a partition of unity or reproduce polynomials. *Analysis-suitable* (AS) T-splines are a mildly restricted class without these problems (Li, Zheng, Sederberg, Hughes and Scott, 2012; [Li and Scott, 2014](https://arxiv.org/abs/1211.5669)).

Each T-junction has an *extension* that continues its missing edge and spans `p` bays in total ([Beirão da Veiga, Buffa, Sangalli and Vázquez, 2013](https://doi.org/10.1142/S0218202513500231)):

- The **face extension** extends `⌊(p+1)/2⌋` bays into the face where the edge is missing.
- The **edge extension** extends `⌈(p-1)/2⌉` bays in the opposite direction.

That is, `(p+1)/2` and `(p-1)/2` bays for odd degrees, and `p/2` bays in both directions for even degrees.

A T-junction is horizontal if its missing edge is horizontal (`⊢`, `⊣`), otherwise vertical (`⊥`, `⊤`). A T-mesh is **analysis-suitable** if no horizontal T-junction extension intersects a vertical T-junction extension. Use [`is_analysis_suitable`](@ref) to check this, and [`tjunction_extensions`](@ref) to inspect the extensions.

For an analysis-suitable T-spline:
- the basis functions are linearly independent for any choice of knot values;
- they form a partition of unity and reproduce all polynomials of degree `p`, so a mesh with control points at the Greville points ([`greville_points`](@ref)) represents an affine geometry exactly, and patch tests are passed;
- they are locally linearly independent: exactly `(p+1)^2` basis functions are non-zero on each Bézier element. This follows from the equivalence of analysis-suitable and *dual-compatible* T-meshes (Beirão da Veiga et al., 2013);
- the approximation is optimal, `O(h^(p+1))` in the L2 norm, as for NURBS;
- nested T-meshes give nested spaces, which is the basis for local refinement.

### Which T-meshes are analysis-suitable

- **Tensor product meshes**, i.e. NURBS (`TSplineMesh(nurbsmesh)`), have no T-junctions.
- **Meshes where all T-junctions have the same orientation.** If, for example, every horizontal line spans the whole domain and only vertical lines end in T-junctions, there are no horizontal extensions, so the mesh is always analysis-suitable. This is the case for layered meshes that get finer towards one side, such as `generate_tspline_patch(:graded_rectangle, ...)`, for any number of elements in the layers.
- **Locally refined meshes where perpendicular T-junctions are far enough apart.** When a region is refined in both directions, the new lines in one direction can be extended beyond the refined region, so that the vertical extensions end before they reach the horizontal ones. Use `is_analysis_suitable` to check the result.

A typical mesh that is **not** analysis-suitable has a horizontal and a vertical T-junction close to each other, e.g. a single face split into four, or a refined block whose corner has T-junctions in both directions. The extensions of these T-junctions cross near the corner. For such meshes, the guarantees above are lost. Elements may have more or fewer than `(p+1)^2` basis functions (e.g. 17 or 15 for cubics), the basis is generally not a partition of unity, and linear independence depends on the knot values ([`is_linearly_independent`](@ref) checks it numerically for given knots). A `BezierGrid` can still be created, but `BezierCellValues` only support square extraction operators, and throw an error on elements that do not have exactly `(p+1)^2` basis functions.

The index vectors of the basis functions are computed by marching along the lines through the anchors (Bazilevs et al., 2010). For even degrees, Beirão da Veiga et al. (2013) instead use the lines spanning the whole anchor face; the two definitions coincide for analysis-suitable T-meshes.

## Grids with T-spline cells

`BezierGrid(mesh::TSplineMesh)` creates one [`TSplineCell`](@ref) per Bézier element. Its nodes are the control points whose basis functions are non-zero on the element. For analysis-suitable meshes, all cells have `(p+1)^2` nodes and the same type. If the number of nodes differs between elements (only for meshes that are not analysis-suitable), the grid has several cell types. In that case:

- `add!(dh, :u, ip)` creates one `SubDofHandler` per cell type and uses a [`TSplineInterpolation`](@ref) with the right number of functions for each.
- Use `ndofs_per_cell(dh, cellid)` or `getnbasefunctions(cv)` (after `reinit!`) instead of `ndofs_per_cell(dh)`. Use `CellIterator(sdh)` for each SubDofHandler instead of `CellIterator(dh)`.

The grid contains the facetsets and nodesets `"left"`, `"right"`, `"bottom"` and `"top"` of the parametric domain. Facet integrals (e.g. Neumann conditions) work with `BezierFacetValues` on these facetsets, as for NURBS.

## Periodic T-splines

The T-mesh may extend outside of the parametric domain (the `domain` keyword of [`TSplineMesh`](@ref)). Bézier elements are then only created inside the domain, and only the basis functions that are non-zero in the domain are kept. This is used for periodic splines: with `generate_tspline_patch(:graded_rectangle, ...; periodic = true)`, the knot lines of each layer are extended periodically, by `p` element widths, outside of the left and right sides. The basis functions are not interpolatory at these sides, and the `p` columns of basis functions that are non-zero on the left side are translations of those on the right side. Tying them with

```julia
add!(ch, PeriodicIGADirichlet(:u, getfacetset(grid, "left"), getfacetset(grid, "right")))
```

gives a periodic solution with `C^(p-1)` continuity across the sides, whereas tying the control points of a mesh with open knot vectors only gives `C^0` continuity.

## 3D T-splines

T-splines in 3D are defined in the same way: the T-mesh consists of boxes `(i0, i1, j0, j1, k0, k1)` in index space, the anchors are the vertices (odd degrees) or cell centres (even degrees), and the local knot vectors are found by marching through the T-mesh in each of the three directions. Each basis function is a product of three univariate B-splines, and each Bézier element (a `TSplineCell{RefHexahedron}`) has a Bernstein basis with `(p+1)^3` functions.

`generate_tspline_patch(:graded_box, nels_x, nels_y, nels_z, p)` generates a box which is layered in the z-direction. Layer `k` has `nels_x[k] × nels_y[k] × nels_z[k]` elements; the planes `z = const` span the whole box, while the planes `x = const` and `y = const` end at the layer interfaces. With `periodic = true` (or a tuple `(px, py)`), the knot vectors are periodic in x (and/or y), as for [Periodic T-splines](@ref). The grid contains the sets `"left"`, `"right"` (x), `"front"`, `"back"` (y) and `"bottom"`, `"top"` (z).

```julia
# 2×2×1 elements at the bottom, 4×4×2 elements at the top, periodic in x
mesh = generate_tspline_patch(:graded_box, [2, 4], [2, 4], [1, 2], 3; periodic = (true, false))
is_analysis_suitable(mesh) # true
grid = BezierGrid(mesh)
```

There is no established definition of T-junction extensions in 3D. Instead, *dual-compatibility* is used as the criterion for analysis-suitability: [`is_dual_compatible`](@ref) checks that for every pair of basis functions with overlapping index ranges, there is a direction in which their local index vectors overlap (are consecutive indices of a common global index vector) and differ. Then the products of the univariate B-spline dual functionals form a dual basis, which gives linear independence, local linear independence (exactly `(p+1)^3` functions per element), partition of unity and optimal approximation. In 2D, dual-compatibility is equivalent to analysis-suitability (Beirão da Veiga et al., 2013).

The graded box is dual-compatible: since the planes `z = const` are full, the index vectors of all basis functions overlap in z, and basis functions with the same z index vector lie in the same layer, where the T-mesh is a tensor product mesh in x and y.

## Boundary conditions

The degrees of freedom of a spline are the coefficients of the basis functions, not nodal values. The following constraints are available for `BezierGrid`s (NURBS and T-splines):

- `Dirichlet` on a facetset prescribes the control points whose basis functions are non-zero on the facets, with the value of `f` evaluated at the control points. This is exact for constant (and, with control points at the Greville points, linear) boundary data, but otherwise only an approximation. It supports time dependent data via `update!`.
- [`ProjectedIGADirichlet`](@ref) computes the boundary values by an L2 projection onto the spline functions on the boundary, which is exact whenever the data can be represented by the spline, and generally more accurate. The projection is computed once, when it is added to the `ConstraintHandler`, so the data can not depend on time.
- [`PeriodicIGADirichlet`](@ref) ties the control points on two identical sides of the domain (up to a translation) with affine constraints. See [Periodic T-splines](@ref) for smooth periodicity.

Ferrite's `ProjectedDirichlet` is not supported for spline fields.

The facets of a `TSplineCell` can not be computed from its nodes, so functions that rely on the grid topology, such as `addfacetset!` and `ExclusiveTopology`, do not work; use the predefined sets instead. `L2Projector` and `apply_analytical_iga!` are not supported for T-spline grids. Parametric dimensions 2 and 3 are supported.

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
