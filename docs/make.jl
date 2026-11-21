using Documenter, DocumenterCitations
using FerriteIGA

const liveserver = "liveserver" in ARGS
# Draft mode skips running the code blocks (faster rebuilds). Pass "nodraft" to run them
# anyway, e.g. `push!(ARGS, "nodraft"); include("docs/liveserver.jl")`.
const draft = liveserver && !("nodraft" in ARGS)

# Plots is used in some tutorials and topic guides; make sure rendering is headless
haskey(ENV, "GKSwstype") || (ENV["GKSwstype"] = "100")
# Transparent plot backgrounds, so that the plots blend in with both the light and the dark
# themes (plots are inverted in the dark themes, see src/assets/custom.css)
using Plots: Plots
Plots.default(; background_color = :transparent)

# Generate tutorials
include("generate.jl")

bibtex_plugin = CitationBibliography(
    joinpath(@__DIR__, "src", "assets", "references.bib"),
    style = :numeric,
)

makedocs(
    sitename = "FerriteIGA.jl",
    format = Documenter.HTML(
        assets = ["assets/custom.css"],
        canonical = "https://ferrite-fem.github.io/FerriteIGA.jl/stable",
        collapselevel = 1,
    ),
    doctest = false,
    modules = [FerriteIGA],
    checkdocs = :exports,
    warnonly = Documenter.except(:missing_docs),
    draft = draft,
    pages = Any[
        "Home" => [
            "Introduction" => "index.md",
            "Getting started" => "getting_started.md",
        ],
        "Tutorials" => [
            "Tutorials overview" => "tutorials/index.md",
            "tutorials/plate_with_hole.md",
            "tutorials/structuralvibrations.md",
            "tutorials/kirchhoff_plate.md",
        ],
        "Topic guides" => [
            "Topic guide overview" => "topics/index.md",
            "topics/splines.md",
            "topics/nurbsmesh.md",
            "topics/beziergrid.md",
            "topics/beziercellvalues.md",
            "topics/export.md",
        ],
        "API reference" => [
            "Reference overview" => "reference/index.md",
            "reference/nurbsmesh.md",
            "reference/mesh_generation.md",
            "reference/beziergrid.md",
            "reference/beziercellvalues.md",
            "reference/export.md",
        ],
        "How-to guides" => [
            "How-to guide overview" => "howto/index.md",
        ],
        "devdocs/index.md",
        "cited-literature.md",
    ],
    plugins = [
        bibtex_plugin,
    ],
)

# Make sure there are no generated VTK files left around from the build
for (root, _, files) in walkdir(joinpath(@__DIR__, "build"))
    for file in files
        any(ext -> endswith(file, ext), (".vtu", ".pvd")) && rm(joinpath(root, file))
    end
end

# Insert a <br> in the side menu between the user-facing sections and the
# developer documentation (same as the Ferrite.jl docs)
for (root, _, files) in walkdir(joinpath(@__DIR__, "build")), file in joinpath.(root, files)
    endswith(file, ".html") || continue
    str = read(file, String)
    str = replace(str, r"""(<li(?: class="[^"]*")?><a class="tocitem" href="[^"]*devdocs/[^"]*">)""" => s"<br>\1")
    write(file, str)
end

# Deploy built documentation
if !liveserver
    deploydocs(
        repo = "github.com/Ferrite-FEM/FerriteIGA.jl.git",
        push_preview = true,
        devbranch = "master",
    )
end
