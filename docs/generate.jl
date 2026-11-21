# Generate the tutorials with Literate.jl
import Literate

TUTORIALS_IN = joinpath(@__DIR__, "src", "literate-tutorials")
TUTORIALS_OUT = joinpath(@__DIR__, "src", "tutorials")
mkpath(TUTORIALS_OUT)

for program in readdir(TUTORIALS_IN; join = true)
    name = basename(program)
    if endswith(program, ".jl")
        if draft
            code = "<< no script output when building as draft >>"
        else
            script = Literate.script(program, TUTORIALS_OUT)
            code = strip(read(script, String))
        end

        # Remove "hidden" lines and the Literate header from the plain program
        line_ending_symbol = occursin("\r\n", code) ? "\r\n" : "\n"
        code_clean = join(filter(x -> !endswith(x, "#hide"), split(code, r"\n|\r\n")), line_ending_symbol)
        code_clean = replace(code_clean, r"^# This file was generated .*$"m => "")
        code_clean = strip(code_clean)

        mdpost(str) = replace(str, "@__CODE__" => code_clean)
        Literate.markdown(program, TUTORIALS_OUT; postprocess = mdpost)
    elseif any(endswith.(program, [".png", ".jpg", ".gif", ".webp"]))
        cp(program, joinpath(TUTORIALS_OUT, name); force = true)
    else
        @warn "ignoring $program"
    end
end

# Remove generated VTK files; they should not be deployed with the docs
cd(TUTORIALS_OUT) do
    foreach(file -> any(ext -> endswith(file, ext), (".vtu", ".pvd")) && rm(file), readdir())
end
