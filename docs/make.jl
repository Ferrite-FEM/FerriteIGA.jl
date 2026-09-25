using Documenter
using FerriteIGA

include("generate.jl")

GENERATEDEXAMPLES = [joinpath("examples", f) for f in (
    "plate_with_hole.md",
    "structuralvibrations.md",
    "kirchhoff_plate.md",
    )]

makedocs(
    sitename = "FerriteIGA",
    format = Documenter.HTML(),
    doctest = false,
    modules = [FerriteIGA],
    checkdocs = :exports,
    warnonly = Documenter.except(:missing_docs),
    pages = Any[
        "Home" => "index.md",
        "Manual" => ["splines.md", "meshes.md", "bezier_extraction.md"],
        "Examples" => GENERATEDEXAMPLES,
        "API" => ["api.md", "nurbsmesh.md"],
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
deploydocs(
    repo = "github.com/Ferrite-FEM/FerriteIGA.jl.git",
    push_preview=true,
    devbranch = "master"
)
