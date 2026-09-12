using FerriteIGA
using FerriteIGA: Bernstein, BSplineBasis
using Ferrite
using Test

include("test_meshes.jl")
include("test_values.jl")
include("test_bsplines.jl")
include("test_bertstein.jl")
include("test_bezier_extraction.jl")
include("test_nurbsmesh.jl")
include("test_l2proj.jl")
include("test_examples.jl")
include("test_utility_functions.jl")
include("test_iterators.jl")