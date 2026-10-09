

module TestPlateWithHoleExample
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/plate_with_hole.jl"))
        end
    end
end
module TestTSplineHeatEquationExample
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "../docs/src/literate/tspline_heat_equation.jl"))
        end
    end
end
