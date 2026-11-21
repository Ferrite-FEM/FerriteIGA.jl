# FerriteIGA.jl

Welcome to the documentation for FerriteIGA.jl! FerriteIGA is an extension of the
[Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) finite element toolbox for
*isogeometric analysis* (IGA). It provides spline (B-spline and NURBS) geometries and
basis functions, and makes them available through the same interface as standard finite
elements in Ferrite. Programs written with FerriteIGA therefore look almost exactly like
Ferrite programs: degrees of freedom, boundary conditions, assembly and post-processing all
work in the same way.

!!! note
    Please help improve this documentation -- if something confuses you, chances are you're
    not alone. It's easy to do as you read along: just click on the "Edit on GitHub" link at
    the top of each page, and then [edit the files directly in your browser](https://help.github.com/articles/editing-files-in-another-user-s-repository/).
    Your changes will be vetted by developers before becoming permanent, so don't worry
    about whether you might say something wrong.

## How the documentation is organized

This high level view of the documentation structure will help you find what you are looking
for. The document is organized as follows[^1]:

 - [**Getting started**](getting_started.md) summarizes how an isogeometric program differs
   from a Ferrite program, and is the best place to start.
 - [**Tutorials**](tutorials/index.md) are thoroughly documented examples which guide you
   through the process of solving partial differential equations using FerriteIGA.
 - [**Topic guides**](topics/index.md) contain more in-depth explanations and discussions
   about isogeometric analysis concepts, and specifically how these are realized in
   FerriteIGA.
 - [**API reference**](reference/index.md) contains the technical API reference of functions
   and types (e.g. the documentation strings).
 - [**How-to guides**](howto/index.md) will guide you through the steps involved in
   addressing common tasks and use-cases. These usually build on top of the tutorials and
   thus assume basic knowledge of how FerriteIGA works.

[^1]: The organization of the document follows the [Diátaxis Framework](https://diataxis.fr),
    in the same way as the [Ferrite documentation](https://ferrite-fem.github.io/Ferrite.jl/stable/).

In addition, the document contains [**Developer documentation**](devdocs/index.md), which
documents the internals of FerriteIGA, and a list of [**Cited literature**](cited-literature.md).

Since FerriteIGA builds on Ferrite, the
[Ferrite documentation](https://ferrite-fem.github.io/Ferrite.jl/stable/) is a useful
companion. In particular, it is recommended to be familiar with the first Ferrite tutorial,
[Heat equation](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/heat_equation/),
before reading this documentation.

## Isogeometric analysis

Isogeometric analysis was introduced by [Hughes2005](@citet) as a Galerkin method whose
basis functions are taken from computer aided design (CAD): B-splines and non-uniform
rational B-splines (NURBS) [Cottrell2009](@cite). The idea is to use the *same* functions to
describe the geometry and to approximate the solution, which gives the method its name.
This has a number of consequences:

 - **Exact geometry.** A geometry designed in CAD is described exactly by the analysis
   model, already on the coarsest mesh. NURBS represent conic sections such as circles,
   cylinders and spheres exactly, so for example the circular hole in the
   [Infinite plate with hole](@ref tutorial-plate-with-hole) tutorial is exact, and the
   geometry is unchanged when the mesh is refined.
 - **Higher continuity.** The smoothness of a spline basis is set by its knot vector. A
   spline of degree ``p`` is, in general, ``C^{p-1}``-continuous between elements, whereas
   standard finite elements are only ``C^0``. Stresses and strains are therefore continuous
   inside a patch, and problems that need second derivatives of the unknown, such as the
   [Kirchhoff plate](@ref kirchhoff_plate), can be discretized directly with a quadratic
   basis.
 - **Better spectral properties.** For the same number of degrees of freedom, smooth spline
   spaces approximate the eigenvalues and eigenmodes of a problem considerably better than
   ``C^0`` finite elements, whose higher frequencies diverge from the exact ones
   [Cottrell2006](@cite). This is shown in the [Structural vibrations](@ref structural_vibrations)
   tutorial.
 - **Flexible refinement.** A mesh can be refined by inserting knots (h-refinement), by
   raising the degree (p-refinement), or by a combination that raises both the degree and the
   continuity (k-refinement), without changing the geometry.

A spline basis function is non-zero over several neighboring elements, which does not fit
directly into a finite element code that evaluates shape functions element by element.
FerriteIGA uses *Bézier extraction* [Borden2011](@cite) to overcome this: on each element the
spline basis is written as a linear combination of Bernstein polynomials, which are defined
on a single reference element exactly like standard finite element shape functions. This is
what allows Ferrite's machinery to be reused, and is explained in the
[BezierGrid](@ref beziergrid_topicguide) and
[BezierCellValues](@ref beziercellvalues_topicguide) topic guides.

## Installation

To use FerriteIGA you first need to install Julia, see <https://julialang.org/> for details.
Installing FerriteIGA can then be done from the Pkg REPL; press `]` at the `julia>` prompt
to enter `pkg>` mode:

```
pkg> add FerriteIGA
```

This will install FerriteIGA and all necessary dependencies, including Ferrite. Press
backspace to get back to the `julia>` prompt. Finally, to load FerriteIGA, use

```julia
using FerriteIGA
```

FerriteIGA re-exports Ferrite (and Tensors), so `using Ferrite` is not needed. You are now
all set to start using FerriteIGA!

### Getting help

If you have questions about FerriteIGA it is suggested to use the `#ferrite-fem` channel on
the [Julia Slack](https://julialang.org/slack/), or the `#Ferrite.jl` stream on
[Zulip](https://julialang.zulipchat.com/). Alternatively you can open an
[issue](https://github.com/Ferrite-FEM/FerriteIGA.jl/issues) on the GitHub repository.

## Contributing to FerriteIGA

FerriteIGA is still under active development. If you find a bug, or have ideas for
improvements, you are encouraged to interact with the developers on the
[FerriteIGA GitHub repository](https://github.com/Ferrite-FEM/FerriteIGA.jl).
