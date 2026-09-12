# FerriteIGA.jl

Isogeometric analysis in the [Ferrite](https://github.com/Ferrite-FEM/Ferrite.jl) finite element library.

The basis is a B-spline or NURBS space, the same functions used for CAD geometry. Degrees of freedom, constraints, and assembly use the Ferrite interface. Bézier extraction converts the splines into Bernstein polynomials so the geometry works with a standard finite element code.

The [Splines](@ref), [Meshes](@ref), and [Bezier extraction](@ref) pages expand the topics below. Examples include the [Infinite plate with hole](@ref) and [structural vibrations](@ref structural_vibrations) of an elastic rod.

## Isogeometric analysis

Hughes, Cottrell, and Bazilevs (2005) introduced isogeometric analysis as a Galerkin method whose basis is taken from CAD, consisting of B-splines and eventually non-uniform rational B-splines (Cottrell, Hughes, and Bazilevs, *Isogeometric Analysis* (Wiley, 2009)). A NURBS patch used for design is then available as an analysis mesh. NURBS can represent curved and conic sections such as cylinders, spheres, and circles exactly.

Continuity is built into the knot vector, and the increased smoothness of IGA is a significant difference from standard finite elements. 

In linear elasticity the stress comes from derivatives of the displacement. In a standard $C^0$ mesh these derivatives jump at element edges, which affects stress concentrations and bending curvature. A higher-order smoothness patch keeps strains and stresses continuous inside the patch. Thin shells that use second derivatives of the displacement need at least continuity $C^1$. The [Infinite plate with hole](@ref) is a plane-stress problem with a circular hole that NURBS represent exactly.

The same smoothness provides another advantage in the eigenvalue spectrum and vibration problems. Cottrell et al. (2006) showed that the spline discretization produces an eigenvalue/frequency spectrum closer to the exact spectrum of the problem than a $C^0$ finite element space of the same degree, which drifts once the mode number exceeds about half the number of degrees of freedom. The [structural vibrations](@ref structural_vibrations) example shows the vibration spectrum of IGA.

The mesh can undergo knot insertion (h-refinement), order elevation (p-refinement), or a combination of both to increase smoothness (k-refinement).

## Bézier extraction

Spline functions overlap several neighbouring elements. A finite element code expects shape functions on a single cell. Bézier extraction (Borden, Scott, Evans, and Hughes, 2011) converts the B-spline or NURBS basis on each element into Bernstein polynomials. Those polynomials are local to the cell and $C^0$, so the isogeometric geometry can be used in Ferrite.

## Installation

The package is unregistered. From the Pkg REPL,

```
pkg> add https://github.com/Ferrite-FEM/FerriteIGA.jl
```

Ferrite is installed as a dependency.

## Quick start

Generate a NURBS patch, convert it to a Bézier mesh, and build cell values.

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

`reinit!` applies the extraction. After that call, `DofHandler`, `ConstraintHandler`, and assembly follow the [Ferrite documentation](https://ferrite-fem.github.io/Ferrite.jl/stable/).

A uniform box can also be built with `generate_grid` and a `BezierCell`.

```julia
grid = generate_grid(BezierCell{RefLine, 2}, (10,), Vec(0.0), Vec(1.0))
```

`generate_nurbs_patch` also provides lines, rectangles, cubes, rings, cylindrical sectors, a plate with a circular hole, and several singly and doubly curved shells. The full list is on the [Meshes](@ref) page.

## References

- T. J. R. Hughes, J. A. Cottrell, and Y. Bazilevs. Isogeometric analysis. CAD, finite elements, NURBS, exact geometry and mesh refinement. *Comput. Methods Appl. Mech. Engrg.*, 194:4135–4195, 2005. [doi:10.1016/j.cma.2004.10.008](https://doi.org/10.1016/j.cma.2004.10.008)
- J. A. Cottrell, A. Reali, Y. Bazilevs, and T. J. R. Hughes. Isogeometric analysis of structural vibrations. *Comput. Methods Appl. Mech. Engrg.*, 195:5257–5296, 2006. [doi:10.1016/j.cma.2005.09.027](https://doi.org/10.1016/j.cma.2005.09.027)
- M. J. Borden, M. A. Scott, J. A. Evans, and T. J. R. Hughes. Isogeometric finite element data structures based on Bézier extraction of NURBS. *Int. J. Numer. Meth. Engng.*, 87:15–47, 2011. [doi:10.1002/nme.2968](https://doi.org/10.1002/nme.2968)
- J. A. Cottrell, T. J. R. Hughes, and Y. Bazilevs. *Isogeometric Analysis. Toward Integration of CAD and FEA*. Wiley, 2009. [doi:10.1002/9780470749081](https://doi.org/10.1002/9780470749081)
