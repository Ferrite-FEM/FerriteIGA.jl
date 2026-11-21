# Tutorials

On this page you find an overview of FerriteIGA tutorials. The tutorials explain and show how
FerriteIGA can be used to solve a range of problems with isogeometric analysis.

The tutorials all follow roughly the same structure:
 - **Introduction** introduces the problem to be solved and discusses the learning outcomes
   of the tutorial.
 - **Commented program** is the code for solving the problem with explanations and comments.
 - **Plain program** is the raw source code of the program.

When studying the tutorials it is a good idea to obtain a local copy of the code and run it
on your own machine as you read along. The tutorials assume that you are familiar with the
basics of Ferrite, as covered by the first
[Ferrite tutorial](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/heat_equation/),
and with the [Getting started](@ref getting-started) page, which summarizes how an isogeometric
program differs from a Ferrite program.

### Tutorial index

The tutorials are listed in roughly increasing order of complexity.

---

```@raw html
<div class="example-row">
<a class="example-thumb" href="plate_with_hole/">
  <img src="plate_with_hole.png" alt="Infinite plate with hole">
</a>
<div class="example-description">
```

#### [Tutorial 1: Infinite plate with hole](@ref tutorial-plate-with-hole)

This tutorial solves a linear elasticity problem: a plate with a circular hole subjected to a
uniaxial traction. The circular hole is
represented exactly by a NURBS patch. The tutorial shows the full isogeometric workflow:
generating a NURBS patch, converting it to a `BezierGrid`, assembling with
`BezierCellValues` and `BezierFacetValues`, and exporting the displacements and the projected
stresses to VTK.

**Keywords**: linear elasticity, NURBS, exact geometry, Neumann boundary conditions,
L2 projection, VTK export.

```@raw html
</div>
</div>
```

---

#### [Tutorial 2: Structural vibrations](@ref structural_vibrations)

This tutorial compares the natural frequencies of an elastic rod computed with quadratic
isogeometric analysis and with quadratic finite elements, reproducing a result from
[Cottrell2006](@citet). The smooth spline basis gives an accurate approximation of the whole
frequency spectrum, whereas the higher frequencies of the ``C^0`` finite element
approximation diverge from the exact ones.

**Keywords**: structural dynamics, generalized eigenvalue problem, spectral accuracy,
comparison with finite elements.

---

#### [Tutorial 3: Kirchhoff plate](@ref kirchhoff_plate)

This tutorial solves the bending of a simply supported Kirchhoff plate under a uniform
pressure. The weak form contains second derivatives of the deflection, which requires a
``C^1``-continuous basis. A quadratic B-spline basis provides this directly, which makes this
problem a natural application of isogeometric analysis. The tutorial shows how to compute
second derivatives with `BezierCellValues`.

**Keywords**: plate bending, fourth-order problem, ``C^1`` continuity, shape function hessians.
