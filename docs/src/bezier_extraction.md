# Bezier extraction

A B-spline or NURBS function of degree $p$ overlaps several neighbouring knot spans. Ferrite evaluates shape functions one cell at a time. Bézier extraction ([Borden, Scott, Evans, and Hughes, 2011](https://doi.org/10.1002/nme.2968)) writes the spline basis on an element in terms of Bernstein polynomials. Those polynomials are $C^0$ and are local one cell, so the geometry can be integrated like a standard finite element mesh. The higher smoothness between cells is kept geometrically.

The extraction operator $\boldsymbol{C}^e$ is built from the knot vector once per cell. `BezierGrid` stores it. `IGAInterpolation` is the Bernstein basis. `BezierCellValues` applies $\boldsymbol{C}^e$ in `reinit!`.

## B-spline basis on an element

On element $e$,

```math
\boldsymbol{N}^e = \boldsymbol{C}^e \boldsymbol{B}^e,
```

where $\boldsymbol{N}^e$ are the B-spline functions on the element and $\boldsymbol{B}^e$ are the Bernstein polynomials of the same degree. $\boldsymbol{C}^e$ is sparse. In 2D and 3D it is built from the 1D operators in each direction (`compute_bezier_extraction_operators`).

## NURBS

NURBS place a weight $w_A$ at each control point. Let $\boldsymbol{W}^e$ be the diagonal matrix of those weights on the element. The rational basis is

```math
\boldsymbol{R}^e = \boldsymbol{W}^e \frac{\boldsymbol{N}^e}{W(\xi)} = \boldsymbol{W}^e \boldsymbol{C}^e \frac{\boldsymbol{B}^e}{W^b(\xi)},
```

with weight functions

```math
W(\xi) = \sum_A N_A(\xi)\, w_A, \qquad
W^b(\xi) = \sum_A B_A(\xi)\, w^b_A.
```

The two sums are equal, $W(\xi) = W^b(\xi)$. Applying $(\boldsymbol{C}^e)^{T}$ to the NURBS data gives the Bézier weights and control points,

```math
\boldsymbol{w}^b = (\boldsymbol{C}^e)^{T} \boldsymbol{w},
\qquad
w^b_A \boldsymbol{X}_{b,A} = \sum_B C^e_{BA}\, w_B \boldsymbol{X}_B.
```

The geometry can then be written in either basis,

```math
\boldsymbol{S}(\xi) = \sum_A R_A(\xi)\, \boldsymbol{X}_A = \sum_A B_A(\xi)\, \boldsymbol{X}_{b,A}.
```

Bernstein values are stored at quadrature points, as in Ferrite. Multiplication by $\boldsymbol{C}^e$ and the weights gives the NURBS values used in the weak form.

## Usage

`BezierGrid(patch)` turns a NURBS patch into a Ferrite grid. Each cell keeps its weights and the extraction operator $\boldsymbol{C}^e$. An ordinary Ferrite grid can be passed as `BezierGrid(grid)`.

`getcoordinates(grid, cellid)` returns a `BezierCoords` for that cell. `x` and `w` are the NURBS points and weights. `xb` and `wb` are the Bézier points and weights. `beo` is the operator.

`reinit!(cv, coords)` uses that data to fill the shape functions.

`IGAInterpolation{shape, p}()` is the Bernstein basis of degree $p$. `shape` is `RefLine`, `RefQuadrilateral`, or `RefHexahedron`. For a vector unknown write `ip^2`.

`BezierCellValues(qr, ip)` stores the basis at the quadrature points in `qr`. `reinit!` applies $\boldsymbol{C}^e`. On faces use `BezierFacetValues` and `reinit!(fv, coords, faceid)`.

```julia
using Ferrite, FerriteIGA

order = 2 # second order NURBS
nels = (20, 10) # number of elements
patch = generate_nurbs_patch(:plate_with_hole, nels, order)

# Convert the NURBS patch to a grid with Bézier extraction operators
grid = BezierGrid(patch)

# Interpolation and shape values (Bernstein polynomials)
ip = IGAInterpolation{RefQuadrilateral, order}()
qr = QuadratureRule{RefQuadrilateral}(4)
cv = BezierCellValues(qr, ip)

# Update cell values
coords = getcoordinates(grid, 1)
reinit!(cv, coords)
```

After `reinit!`, `shape_value` and `shape_gradient` give NURBS values. Assembly then follows Ferrite (`DofHandler`, `ConstraintHandler`, `start_assemble`). See the [Ferrite manual](https://ferrite-fem.github.io/Ferrite.jl/stable/).

Displacement fields on a `BezierGrid` are written with `VTKIGAFile`.
```julia
VTKIGAFile("plate_with_hole.vtu", grid) do vtk
    write_solution(vtk, dh, u)
end
```
