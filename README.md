# FerriteIGA.jl

Isogeometric analysis in the [Ferrite](https://github.com/Ferrite-FEM/Ferrite.jl) finite element library.

## Documentation

[![][docs-dev-img]][docs-dev-url]

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://ferrite-fem.github.io/FerriteIGA.jl/dev/

## Isogeometric analysis

Hughes, Cottrell, and Bazilevs (2005) introduced isogeometric analysis as a Galerkin method whose basis is taken from CAD, consisting of B-splines and non-uniform rational B-splines (Cottrell, Hughes, and Bazilevs, *Isogeometric Analysis* (Wiley, 2009)). A NURBS patch from CAD design can then be used as the analysis mesh. NURBS can represent curved and conic sections such as cylinders, spheres, and circles exactly. The increased smoothness of IGA is a significant difference from traditional finite elements.

Spline functions overlap several neighbouring elements. A finite element code expects shape functions on a single cell. Bézier extraction (Borden, Scott, Evans, and Hughes, 2011) converts the B-spline or NURBS basis on each element into Bernstein polynomials. Those polynomials are local to the cell so the CAD geometry can be used in Ferrite. 

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

`reinit!` applies the extraction. After that call the cell values are used as in Ferrite.

A uniform box can also be built with `generate_grid` and a `BezierCell`.

```julia
grid = generate_grid(BezierCell{RefLine, 2}, (10,), Vec(0.0), Vec(1.0))
```

`generate_nurbs_patch` also provides lines, rectangles, cubes, rings, cylindrical sectors, a plate with a circular hole, and several singly and doubly curved shells.

## References

- T. J. R. Hughes, J. A. Cottrell, and Y. Bazilevs. Isogeometric analysis. CAD, finite elements, NURBS, exact geometry and mesh refinement. *Comput. Methods Appl. Mech. Engrg.*, 194:4135–4195, 2005. [doi:10.1016/j.cma.2004.10.008](https://doi.org/10.1016/j.cma.2004.10.008)
- J. A. Cottrell, A. Reali, Y. Bazilevs, and T. J. R. Hughes. Isogeometric analysis of structural vibrations. *Comput. Methods Appl. Mech. Engrg.*, 195:5257–5296, 2006. [doi:10.1016/j.cma.2005.09.027](https://doi.org/10.1016/j.cma.2005.09.027)
- M. J. Borden, M. A. Scott, J. A. Evans, and T. J. R. Hughes. Isogeometric finite element data structures based on Bézier extraction of NURBS. *Int. J. Numer. Meth. Engng.*, 87:15–47, 2011. [doi:10.1002/nme.2968](https://doi.org/10.1002/nme.2968)
- J. A. Cottrell, T. J. R. Hughes, and Y. Bazilevs. *Isogeometric Analysis. Toward Integration of CAD and FEA*. Wiley, 2009. [doi:10.1002/9780470749081](https://doi.org/10.1002/9780470749081)
