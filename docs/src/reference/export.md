```@meta
CurrentModule = FerriteIGA
DocTestSetup = :(using FerriteIGA)
```

# Export

FerriteIGA writes VTK files in which every cell is an exact rational Bézier cell, so that the
curved geometry and the spline solution are visualized without additional subdivision. See
the [Export topic guide](@ref export_topicguide) for more details.

```@docs
VTKIGAFile
```

The following Ferrite functions can be used with a [`VTKIGAFile`](@ref):

 - `write_solution(vtk::VTKIGAFile, dh, u)`
 - `write_projection(vtk::VTKIGAFile, proj::L2Projector, vals, name)`
 - `write_node_data(vtk::VTKIGAFile, nodedata, name)`
 - `write_cell_data(vtk::VTKIGAFile, celldata, name)`

See the
[Ferrite export reference](https://ferrite-fem.github.io/Ferrite.jl/stable/reference/export/)
for the documentation of these functions.
