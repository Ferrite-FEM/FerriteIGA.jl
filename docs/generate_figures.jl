# Generates the figures of the tutorials (docs/src/nurbs.md and docs/src/tsplines.md).
#
# The figures are computed with FerriteIGA itself (T-meshes, local knot vectors, T-junction extensions,
# Bézier elements, ...), so they always show what the implementation does.
#
# Run from make.jl, or standalone with `julia --project=docs docs/generate_figures.jl`.
# Set ENV["FERRITEIGA_FIGURES_PNG"] = "some/dir" to also write png-copies (handy for a quick look).

using Plots
using Plots: mm
using FerriteIGA, Ferrite

module TutorialFigures

using Plots
using Plots: mm
using FerriteIGA, Ferrite

const FIGDIR = joinpath(@__DIR__, "src", "figures")
const PNGDIR = get(ENV, "FERRITEIGA_FIGURES_PNG", "")

# # #
# Style
# # #

# Categorical colors (fixed order), a sequential blue ramp and a diverging blue-gray-red ramp
col(hex) = parse(Plots.Colorant, hex)
const SERIES  = col.(["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"])
const BLUE, ORANGE, AQUA = SERIES[1], SERIES[2], SERIES[3]
const INK     = col("#0b0b0b")   # primary text and mesh lines
const INK2    = col("#52514e")   # secondary text
const MUTED   = col("#a3a29c")   # de-emphasized marks
const FAINT   = col("#e4e3df")   # hairlines
const TINT    = col("#cde2fb")   # light fill
const SURFACE = col("#fcfcfb")
const SEQ     = cgrad(col.(["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]))
const DIV     = cgrad(col.(["#1c5cab", "#86b6ef", "#f0efec", "#ec9a9a", "#b83232"]))
const FONT    = "DejaVu Sans"    # has the subscripts and Greek letters used in the labels

function setup_style()
    gr()
    default(
        fontfamily = FONT, background_color = SURFACE, foreground_color = INK,
        foreground_color_axis = MUTED, foreground_color_border = MUTED, foreground_color_text = INK2,
        foreground_color_guide = INK2, guidefontsize = 9, tickfontsize = 8, titlefontsize = 9,
        titlefontcolor = INK, legendfontsize = 8, grid = false, legend = false, linewidth = 2,
        framestyle = :axes, tick_direction = :out, margin = 3mm, dpi = 200,
    )
end

function save(plt, name)
    mkpath(FIGDIR)
    savefig(plt, joinpath(FIGDIR, name * ".svg"))
    if !isempty(PNGDIR)
        mkpath(PNGDIR)
        savefig(plt, joinpath(PNGDIR, name * ".png"))
    end
    return nothing
end

label(str; color = INK, size = 8, align = :center) = text(str, Plots.font(family = FONT, pointsize = size, color = color, halign = align))

# Unicode subscripts, e.g. "N" * sub(12) == "N₁₂"
sub(i::Integer) = join(Char(0x2080 + d) for d in reverse(digits(i)))

# # #
# Univariate B-splines
# # #

# The B-spline function with the local knot vector `K` (p+2 knots), evaluated with the Cox-de Boor recursion.
function bspline(K::AbstractVector, x::Real)
    p = length(K) - 2
    (K[1] <= x <= K[end]) || return 0.0
    if p == 0
        return (K[1] <= x < K[2]) ? 1.0 : 0.0
    end
    a = K[p+1] > K[1] ? (x - K[1]) / (K[p+1] - K[1]) * bspline(view(K, 1:p+1), x) : 0.0
    b = K[p+2] > K[2] ? (K[p+2] - x) / (K[p+2] - K[2]) * bspline(view(K, 2:p+2), x) : 0.0
    return a + b
end

# Function number i of the knot vector `kv`. The last element is closed on the right.
function bspline(kv::AbstractVector, p::Int, i::Int, x::Real)
    x >= kv[end] && (x = prevfloat(float(kv[end])))
    return bspline(view(kv, i:i+p+1), x)
end

nfunctions(kv, p) = length(kv) - p - 1
allbasis(kv, p, x) = [bspline(kv, p, i, x) for i in 1:nfunctions(kv, p)]

# Sample points with the knots included, such that kinks are drawn sharply
function samples(kv; n = 60)
    u = unique(kv)
    return unique!(reduce(vcat, [collect(range(u[k], u[k+1], length = n)) for k in 1:length(u)-1]))
end

# Label of a knot value with its multiplicity, e.g. "0,0,0"
knotlabel(kv, v) = join(fill(isinteger(v) ? string(Int(v)) : string(v), count(==(v), kv)), ",")
knotticks(kv) = (unique(kv), [knotlabel(kv, v) for v in unique(kv)])

function plot_basis!(plt, kv, p; colors = i -> SERIES[mod1(i, 8)], labels = true, lw = 2, sp = 1, name = "N", labelsize = 8, subset = 1:nfunctions(kv, p))
    xs = samples(kv)
    for i in subset
        ys = [bspline(kv, p, i, x) for x in xs]
        plot!(plt, xs, ys; color = colors(i), lw, subplot = sp)
        if labels
            k = argmax(ys)
            annotate!(plt, xs[k], ys[k] + 0.07, label(name * sub(i); size = labelsize), subplot = sp)
        end
    end
    return plt
end

# # #
# Figures: B-splines and NURBS
# # #

function fig_cox_de_boor()
    kv = collect(0.0:5.0)
    plt = plot(layout = (3, 1), size = (680, 520), link = :x, left_margin = 6mm)
    for (sp, p) in enumerate(0:2)
        xs = samples(kv)
        for i in 1:nfunctions(kv, p)
            if p == 0 #Piecewise constants: draw the half-open intervals
                plot!(plt, [kv[i], kv[i+1]], [1, 1]; color = SERIES[i], subplot = sp)
                scatter!(plt, [kv[i]], [1]; color = SERIES[i], markerstrokecolor = SERIES[i], markersize = 4, subplot = sp)
                scatter!(plt, [kv[i+1]], [1]; color = SURFACE, markerstrokecolor = SERIES[i], markersize = 4, subplot = sp)
                annotate!(plt, (kv[i] + kv[i+1]) / 2, 1.13, label("N" * sub(i) * ",₀"), subplot = sp)
            else
                ys = [bspline(view(kv, i:i+p+1), x) for x in xs]
                plot!(plt, xs, ys; color = SERIES[i], subplot = sp)
                k = argmax(ys)
                annotate!(plt, xs[k], ys[k] + 0.13, label("N" * sub(i) * "," * sub(p)), subplot = sp)
            end
        end
        plot!(plt; subplot = sp, ylims = (0, 1.25), yticks = [0, 1], xlims = (-0.15, 5.15),
            xticks = (kv, ["ξ" * sub(i) for i in 1:6]),
            ylabel = ["p = 0", "p = 1", "p = 2"][sp])
    end
    save(plt, "nurbs_cox_de_boor")
end

const OPEN_KV = [0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5.0] #Quadratic example from Bazilevs et al. (2010), Fig. 12

function fig_open_basis()
    kv, p = OPEN_KV, 2
    plt = plot(size = (680, 300), ylims = (0, 1.15), yticks = [0, 0.5, 1], xticks = knotticks(kv),
        xlabel = "ξ  (the tick labels show the knots, repeated knots are written out)", bottom_margin = 5mm)
    plot!(plt, [0, 5], [1, 1]; color = FAINT, lw = 1)
    plot_basis!(plt, kv, p)
    save(plt, "nurbs_open_basis")
end

function fig_multiplicity()
    p = 3
    kvs = ([0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4.0], [0, 0, 0, 0, 1, 2, 2, 3, 4, 4, 4, 4.0], [0, 0, 0, 0, 1, 2, 2, 2, 3, 4, 4, 4, 4.0])
    titles = ["knot 2 once:  C² at ξ = 2", "knot 2 twice:  C¹ at ξ = 2", "knot 2 three times:  C⁰ at ξ = 2"]
    plt = plot(layout = (1, 3), size = (900, 260), bottom_margin = 6mm, left_margin = 4mm)
    for (sp, kv) in enumerate(kvs)
        vline!(plt, [2.0]; color = FAINT, lw = 1, subplot = sp)
        #Emphasize the functions which are non-zero at the repeated knot
        active = [i for i in 1:nfunctions(kv, p) if bspline(kv, p, i, 2.0) > 1e-12]
        plot_basis!(plt, kv, p; colors = i -> i in active ? BLUE : MUTED, labels = false, sp)
        plot!(plt; subplot = sp, title = titles[sp], ylims = (0, 1.08), yticks = [0, 1], xticks = knotticks(kv), xlabel = "ξ")
    end
    save(plt, "nurbs_multiplicity")
end

const CUBIC_KV = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5.0]

function fig_local_support()
    kv, p, A = CUBIC_KV, 3, 4
    plt = plot(size = (680, 300), ylims = (-0.16, 1.1), yticks = [0, 0.5, 1], xticks = knotticks(kv), xlabel = "ξ", bottom_margin = 5mm)
    plot_basis!(plt, kv, p; colors = i -> MUTED, labels = false, lw = 1.5)
    xs = samples(kv)
    ys = [bspline(kv, p, A, x) for x in xs]
    plot!(plt, xs, ys; color = BLUE, fillrange = 0, fillalpha = 0.15, lw = 2.5)
    annotate!(plt, 2.0, maximum(ys) + 0.07, label("N₄"))
    #The local knot vector
    K = kv[A:A+p+1]
    plot!(plt, [K[1], K[end]], [-0.08, -0.08]; color = BLUE, lw = 2)
    scatter!(plt, K, fill(-0.08, length(K)); color = BLUE, markerstrokecolor = SURFACE, markersize = 6)
    annotate!(plt, 4.12, -0.08, label("local knots"; align = :left, color = INK2))
    save(plt, "nurbs_local_support")
end

const CURVE_POINTS = [(0.0, 0.3), (0.1, 1.9), (1.6, 2.3), (2.5, 0.2), (3.7, 0.3), (4.5, 2.2), (5.3, 0.9), (6.4, 1.3)]

curve_point(kv, p, P, x; w = ones(length(P))) = (N = allbasis(kv, p, x) .* w; (sum(N .* first.(P)) / sum(N), sum(N .* last.(P)) / sum(N)))

function plot_polygon!(plt, P; color = MUTED, labels = true, sp = 1, offset = (0.0, 0.16), markersize = 5)
    plot!(plt, first.(P), last.(P); color, lw = 1, label = "", subplot = sp)
    scatter!(plt, first.(P), last.(P); color = SURFACE, markerstrokecolor = color, markerstrokewidth = 1.5, markersize, label = "", subplot = sp)
    if labels
        for (i, pt) in enumerate(P)
            annotate!(plt, pt[1] + offset[1], pt[2] + offset[2], label("P" * sub(i); color = INK2), subplot = sp)
        end
    end
end

function fig_curve()
    kv, p, P = OPEN_KV, 2, CURVE_POINTS
    plt = plot(size = (680, 330), aspect_ratio = :equal, framestyle = :none, xlims = (-0.4, 6.8), ylims = (-0.1, 2.7))
    plot_polygon!(plt, P)
    pts = [curve_point(kv, p, P, x) for x in samples(kv)]
    plot!(plt, first.(pts), last.(pts); color = BLUE, lw = 2.5)
    #Images of the knots: the element boundaries
    offsets = Dict(0 => (0.0, -0.17), 1 => (0.05, -0.17), 2 => (0.3, 0.05), 3 => (0.0, -0.17), 4 => (0.32, 0.02), 5 => (0.0, -0.17))
    for v in unique(kv)
        x, y = curve_point(kv, p, P, v)
        scatter!(plt, [x], [y]; color = INK, markerstrokecolor = SURFACE, markersize = 5)
        dx, dy = offsets[Int(v)]
        annotate!(plt, x + dx, y + dy, label("ξ=$(Int(v))"; color = INK2))
    end
    save(plt, "nurbs_curve")
end

function fig_local_control()
    kv, p, P = OPEN_KV, 2, CURVE_POINTS
    A = 4
    Pnew = copy(P)
    Pnew[A] = (P[A][1] + 0.5, P[A][2] - 0.75)
    plt = plot(size = (680, 360), aspect_ratio = :equal, framestyle = :none, xlims = (-0.4, 6.8), ylims = (-0.75, 2.7))
    plot_polygon!(plt, P; labels = false, color = FAINT)
    plot_polygon!(plt, Pnew; labels = false)
    xs = samples(kv)
    old = [curve_point(kv, p, P, x) for x in xs]
    new = [curve_point(kv, p, Pnew, x) for x in xs]
    plot!(plt, first.(old), last.(old); color = MUTED, lw = 2)
    #The curve only changes in the support of N_A
    a, b = kv[A], kv[A+p+1]
    for (range, color) in (((kv[1], a), BLUE), ((a, b), ORANGE), ((b, kv[end]), BLUE))
        idx = findall(x -> range[1] <= x <= range[2], xs)
        plot!(plt, first.(new[idx]), last.(new[idx]); color, lw = 2.5)
    end
    plot!(plt, [P[A][1], Pnew[A][1]], [P[A][2], Pnew[A][2]]; color = INK, lw = 1.5, arrow = true)
    annotate!(plt, P[A][1] - 0.05, P[A][2] + 0.17, label("P₄"; color = INK2))
    annotate!(plt, Pnew[A][1] + 0.15, Pnew[A][2] - 0.02, label("new position of P₄"; color = INK2, align = :left))
    for v in (a, b)
        x, y = curve_point(kv, p, Pnew, v)
        scatter!(plt, [x], [y]; color = INK, markerstrokecolor = SURFACE, markersize = 5)
        annotate!(plt, x + 0.02, y + 0.17, label("ξ=$(Int(v))"; color = INK2))
    end
    save(plt, "nurbs_local_control")
end

function fig_knot_insertion()
    kv, p, P = OPEN_KV, 2, CURVE_POINTS
    mesh = NURBSMesh((copy(kv),), (p,), [Vec(pt) for pt in P])
    refined = knotinsertion(knotinsertion(mesh, 1.5; dir = 1), 2.5; dir = 1)
    Pr = [Tuple(x) for x in refined.control_points]
    kvr = refined.knot_vectors[1]
    plt = plot(layout = (1, 2), size = (900, 250))
    for (sp, (kv_, P_, title)) in enumerate(((kv, P, "before: 8 control points, 5 elements"), (kvr, Pr, "after inserting ξ = 1.5 and 2.5: 10 control points, 7 elements")))
        sp == 2 && plot_polygon!(plt, P; labels = false, color = FAINT, sp)
        plot_polygon!(plt, P_; labels = false, sp, color = sp == 1 ? MUTED : ORANGE)
        pts = [curve_point(kv_, p, P_, x) for x in samples(kv_)]
        plot!(plt, first.(pts), last.(pts); color = BLUE, lw = 2.5, subplot = sp)
        for v in unique(kv_)
            x, y = curve_point(kv_, p, P_, v)
            scatter!(plt, [x], [y]; color = INK, markerstrokecolor = SURFACE, markersize = 4, subplot = sp)
        end
        plot!(plt; subplot = sp, title, aspect_ratio = :equal, framestyle = :none, xlims = (-0.4, 6.8), ylims = (-0.1, 2.7))
    end
    save(plt, "nurbs_knot_insertion")
end

function fig_refinement_basis()
    p = 3
    kv = CUBIC_KV
    kvr = sort([kv; 1.5])
    plt = plot(layout = (2, 1), size = (680, 430), link = :x, left_margin = 4mm)
    xs = samples(kvr)
    plot_basis!(plt, kv, p; colors = i -> i == 4 ? BLUE : MUTED, labels = false, sp = 1, lw = 1.5)
    plot!(plt, xs, [bspline(kv, p, 4, x) for x in xs]; color = BLUE, lw = 2.5, subplot = 1)
    annotate!(plt, 2.0, 0.76, label("N₄"), subplot = 1)
    plot_basis!(plt, kvr, p; colors = i -> MUTED, labels = false, sp = 2, lw = 1.5)
    #N_4 = 1/2 Nbar_4 + 5/6 Nbar_5 (the row of the knot insertion matrix)
    c4 = (1.5 - kv[4]) / (kv[4+p] - kv[4])
    c5 = (kv[5+p] - 1.5) / (kv[5+p] - kv[5])
    plot!(plt, xs, [bspline(kv, p, 4, x) for x in xs]; color = BLUE, lw = 1.5, ls = :dash, subplot = 2)
    plot!(plt, xs, [c4 * bspline(kvr, p, 4, x) for x in xs]; color = ORANGE, lw = 2.5, subplot = 2)
    plot!(plt, xs, [c5 * bspline(kvr, p, 5, x) for x in xs]; color = AQUA, lw = 2.5, subplot = 2)
    annotate!(plt, 1.5, 0.3, label("½ · N̄₄"), subplot = 2)
    annotate!(plt, 2.27, 0.46, label("⅚ · N̄₅"), subplot = 2)
    annotate!(plt, 2.0, 0.95, label("dashed:  N₄ = ½ · N̄₄ + ⅚ · N̄₅"), subplot = 2)
    plot!(plt; subplot = 1, ylims = (0, 1.08), yticks = [0, 1], xticks = knotticks(kv), ylabel = "original basis")
    plot!(plt; subplot = 2, ylims = (0, 1.08), yticks = [0, 1], xticks = knotticks(kvr), ylabel = "refined basis", xlabel = "ξ")
    save(plt, "nurbs_refinement_basis")
end

function fig_circle()
    kv, p = [0, 0, 0, 1, 1, 1.0], 2
    P = [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    weights = [0.25, sqrt(2) / 2, 1.0, 3.0]
    names = ["w₂ = 0.25", "w₂ = √2/2 ≈ 0.71  (exact circle)", "w₂ = 1  (B-spline: a parabola)", "w₂ = 3"]
    ramp = col.(["#86b6ef", "#2a78d6", "#1c5cab", "#0d366b"]) #One hue, ordered by the weight
    plt = plot(layout = grid(1, 2, widths = [0.5, 0.5]), size = (900, 400))
    #Left: the curves, on top of the exact quarter circle (wide grey band)
    θ = range(0, π / 2, length = 100)
    plot!(plt, cos.(θ), sin.(θ); color = FAINT, lw = 10, label = "", subplot = 1)
    plot_polygon!(plt, P; labels = false, sp = 1, markersize = 6)
    xs = range(0, 1, length = 200)
    for (k, w2) in enumerate(weights)
        pts = [curve_point(kv, p, P, x; w = [1.0, w2, 1.0]) for x in xs]
        plot!(plt, first.(pts), last.(pts); color = ramp[k], lw = 2, label = names[k], subplot = 1)
    end
    annotate!(plt, 1.04, 1.04, label("P₂ (weight w₂)"; align = :left, color = INK2), subplot = 1)
    annotate!(plt, 1.04, 0.0, label("P₁"; align = :left, color = INK2), subplot = 1)
    annotate!(plt, 0.0, 1.07, label("P₃"; color = INK2), subplot = 1)
    plot!(plt; subplot = 1, aspect_ratio = :equal, framestyle = :none, xlims = (-0.1, 1.5), ylims = (-0.08, 1.15),
        legend = :bottomleft, foreground_color_legend = nothing, background_color_legend = nothing)
    #Right: the rational basis functions for w₂ = 3 compared with the B-spline basis
    for i in 1:3
        plot!(plt, xs, [bspline(kv, p, i, x) for x in xs]; color = MUTED, lw = 1.5, subplot = 2)
        R = [(N = allbasis(kv, p, x) .* [1.0, 3.0, 1.0]; N[i] / sum(N)) for x in xs]
        plot!(plt, xs, R; color = SERIES[i], lw = 2.5, subplot = 2)
    end
    annotate!(plt, 0.5, 0.8, label("R₂"), subplot = 2)
    annotate!(plt, 0.09, 0.85, label("R₁"), subplot = 2)
    annotate!(plt, 0.91, 0.85, label("R₃"), subplot = 2)
    annotate!(plt, 0.5, 0.45, label("N₂"; color = INK2), subplot = 2)
    annotate!(plt, 0.2, 0.72, label("N₁"; color = INK2), subplot = 2)
    annotate!(plt, 0.8, 0.72, label("N₃"; color = INK2), subplot = 2)
    plot!(plt; subplot = 2, title = "w = (1, 3, 1):  rational basis R in colour,  B-splines N in grey",
        xlabel = "ξ", ylims = (0, 1.08), yticks = [0, 0.5, 1], xticks = [0, 0.5, 1], bottom_margin = 6mm, left_margin = 8mm)
    save(plt, "nurbs_circle")
end

const TP_KV = [0, 0, 0, 1, 2, 3, 4, 4, 4.0]

function fig_tensor_basis()
    kv, p, i, j = TP_KV, 2, 3, 4
    xs = collect(range(0, 4, length = 161))
    Ni = [bspline(kv, p, i, x) for x in xs]
    Nj = [bspline(kv, p, j, y) for y in xs]
    plt = plot(layout = grid(2, 2, heights = [0.28, 0.72], widths = [0.72, 0.28]), size = (640, 600))
    plot!(plt, xs, Ni; color = BLUE, subplot = 1, ylims = (0, 0.85), yticks = [0, 0.5], xticks = knotticks(kv), xlims = (0, 4),
        title = "N₃(ξ)", bottom_margin = 0mm)
    plot!(plt; subplot = 2, framestyle = :none)
    annotate!(plt, 0.5, 0.5, label("N(ξ, η) = N₃(ξ) · N₄(η)"; color = INK2), subplot = 2)
    heatmap!(plt, xs, xs, [Nj[b] * Ni[a] for b in eachindex(xs), a in eachindex(xs)]; color = SEQ, subplot = 3, colorbar = false,
        xlims = (0, 4), ylims = (0, 4), xticks = knotticks(kv), yticks = knotticks(kv), xlabel = "ξ", ylabel = "η", aspect_ratio = :equal)
    for v in unique(kv)
        plot!(plt, [v, v], [0, 4]; color = INK, lw = 0.7, subplot = 3)
        plot!(plt, [0, 4], [v, v]; color = INK, lw = 0.7, subplot = 3)
    end
    #The support
    plot!(plt, [kv[i], kv[i+p+1], kv[i+p+1], kv[i], kv[i]], [kv[j], kv[j], kv[j+p+1], kv[j+p+1], kv[j]]; color = ORANGE, lw = 2.5, subplot = 3)
    plot!(plt, Nj, xs; color = BLUE, subplot = 4, xlims = (0, 0.85), xticks = [0, 0.5], yticks = knotticks(kv), ylims = (0, 4),
        title = "N₄(η)", left_margin = 0mm)
    save(plt, "nurbs_tensor_basis")
end

function fig_tensor_surface()
    kv, p, i, j = TP_KV, 2, 3, 4
    xs = collect(range(0, 4, length = 81))
    Z = [bspline(kv, p, j, y) * bspline(kv, p, i, x) for y in xs, x in xs]
    plt = surface(xs, xs, Z; color = SEQ, colorbar = false, camera = (35, 40), size = (640, 420), xlabel = "ξ", ylabel = "η",
        zlims = (0, 0.6), zticks = [0, 0.25, 0.5], xticks = 0:4, yticks = 0:4, framestyle = :axes, linewidth = 0.3, linecolor = col("#1c5cab"))
    save(plt, "nurbs_tensor_surface")
end

# A quadratic NURBS patch of a quarter annulus, refined by knot insertion
function quarter_annulus()
    s = sqrt(2) / 2
    kv = [0, 0, 0, 1, 1, 1.0]
    ri, ro = 1.0, 2.0
    cps = Vec{2,Float64}[]
    ws = Float64[]
    for r in (ri, (ri + ro) / 2, ro), (c, w) in (((0.0, 1.0), 1.0), ((1.0, 1.0), s), ((1.0, 0.0), 1.0))
        push!(cps, Vec(r * c[1], r * c[2]))
        push!(ws, w)
    end
    mesh = NURBSMesh((copy(kv), copy(kv)), (2, 2), cps, ws)
    for ξ in (1 / 3, 2 / 3)
        mesh = knotinsertion(mesh, ξ; dir = 1)
    end
    return knotinsertion(mesh, 0.5; dir = 2)
end

function nurbs_point(mesh::NURBSMesh{2}, ξ, η)
    (kv1, kv2), (p, q) = mesh.knot_vectors, mesh.orders
    N1, N2 = allbasis(kv1, p, ξ), allbasis(kv2, q, η)
    x, W = zero(Vec{2,Float64}), 0.0
    for j in eachindex(N2), i in eachindex(N1)
        A = (j - 1) * length(N1) + i
        R = N1[i] * N2[j] * mesh.weights[A]
        x += R * mesh.control_points[A]
        W += R
    end
    return x / W
end

function fig_spaces()
    mesh = quarter_annulus()
    kv1, kv2 = mesh.knot_vectors
    n1, n2 = length(kv1), length(kv2)
    plt = plot(layout = grid(1, 3, widths = [0.38, 0.28, 0.34]), size = (1000, 380))
    #Index space: all knots are drawn with equal spacing. Cells between different knot values are elements.
    for j in 1:n2-1, i in 1:n1-1
        if kv1[i] < kv1[i+1] && kv2[j] < kv2[j+1]
            plot!(plt, Shape([i, i + 1, i + 1, i], [j, j, j + 1, j + 1]); color = TINT, linecolor = nothing, subplot = 1)
        end
    end
    for i in 1:n1
        plot!(plt, [i, i], [1, n2]; color = INK, lw = 1, subplot = 1)
    end
    for j in 1:n2
        plot!(plt, [1, n1], [j, j]; color = INK, lw = 1, subplot = 1)
    end
    frac(v) = v == 0 ? "0" : v == 1 ? "1" : v ≈ 1 / 3 ? "⅓" : v ≈ 2 / 3 ? "⅔" : "½"
    
    plot!(plt; subplot = 1, title = "index space", aspect_ratio = :equal, framestyle = :axes, xlims = (0.6, n1 + 0.4), ylims = (0.6, n2 + 0.4),
        xticks = (1:n1, ["$i\n$(frac(kv1[i]))" for i in 1:n1]), yticks = (1:n2, ["$j  ($(frac(kv2[j])))" for j in 1:n2]),
        xlabel = "knot index i  (second row: knot value ξ)", ylabel = "knot index j  (knot value η)", bottom_margin = 10mm, left_margin = 6mm)
    #Parameter space
    plot!(plt, Shape([0, 1, 1, 0], [0, 0, 1, 1]); color = TINT, linecolor = nothing, subplot = 2)
    for v in unique(kv1)
        plot!(plt, [v, v], [0, 1]; color = INK, lw = 1, subplot = 2)
    end
    for v in unique(kv2)
        plot!(plt, [0, 1], [v, v]; color = INK, lw = 1, subplot = 2)
    end
    plot!(plt; subplot = 2, title = "parameter space", aspect_ratio = :equal, xlims = (-0.05, 1.05), ylims = (-0.05, 1.05),
        xticks = (unique(kv1), frac.(unique(kv1))), yticks = (unique(kv2), frac.(unique(kv2))), xlabel = "ξ", ylabel = "η")
    #Physical space: images of the knot lines, and the control net
    ts = range(0, 1, length = 60)
    θ = range(0, π / 2, length = 60)
    plot!(plt, Shape([cos.(θ); 2 .* cos.(reverse(θ))], [sin.(θ); 2 .* sin.(reverse(θ))]); color = TINT, linecolor = nothing, subplot = 3)
    for v in unique(kv1)
        pts = [nurbs_point(mesh, v, t) for t in ts]
        plot!(plt, first.(pts), last.(pts); color = INK, lw = 1, subplot = 3)
    end
    for v in unique(kv2)
        pts = [nurbs_point(mesh, t, v) for t in ts]
        plot!(plt, first.(pts), last.(pts); color = INK, lw = 1, subplot = 3)
    end
    m1, m2 = n1 - 3, n2 - 3
    cp(i, j) = mesh.control_points[(j-1)*m1+i]
    for j in 1:m2
        plot!(plt, [cp(i, j)[1] for i in 1:m1], [cp(i, j)[2] for i in 1:m1]; color = ORANGE, lw = 1, subplot = 3)
    end
    for i in 1:m1
        plot!(plt, [cp(i, j)[1] for j in 1:m2], [cp(i, j)[2] for j in 1:m2]; color = ORANGE, lw = 1, subplot = 3)
    end
    scatter!(plt, first.(mesh.control_points), last.(mesh.control_points); color = SURFACE, markerstrokecolor = ORANGE, markerstrokewidth = 1.5, markersize = 4, subplot = 3)
    plot!(plt; subplot = 3, title = "physical space (control net in orange)", aspect_ratio = :equal, framestyle = :none, xlims = (-0.1, 2.15), ylims = (-0.1, 2.15))
    save(plt, "nurbs_spaces")
end

function fig_global_refinement()
    plt = plot(layout = (1, 2), size = (760, 360))
    n = 5
    for sp in 1:2
        plot!(plt, Shape([3, 4, 4, 3], [3, 3, 4, 4]); color = TINT, linecolor = nothing, subplot = sp)
        for k in 0:n
            plot!(plt, [k, k], [0, n]; color = INK, lw = 1, subplot = sp)
            plot!(plt, [0, n], [k, k]; color = INK, lw = 1, subplot = sp)
        end
        if sp == 1
            plot!(plt, [3.5, 3.5], [0, n]; color = ORANGE, lw = 2.5, subplot = sp)
            plot!(plt, [0, n], [3.5, 3.5]; color = ORANGE, lw = 2.5, subplot = sp)
        else
            plot!(plt, [3.5, 3.5], [3, 4]; color = ORANGE, lw = 2.5, subplot = sp)
            plot!(plt, [3, 4], [3.5, 3.5]; color = ORANGE, lw = 2.5, subplot = sp)
        end
        plot!(plt; subplot = sp, aspect_ratio = :equal, framestyle = :none, xlims = (-0.2, n + 0.2), ylims = (-0.2, n + 0.2),
            title = ["what a NURBS patch can do: new knot lines cross the whole patch", "what we would like: split only the blue element"][sp])
    end
    save(plt, "nurbs_global_refinement")
end

function fig_bezier()
    kv, p = [0, 0, 0, 1, 2, 3, 3, 3.0], 2
    plt = plot(layout = (2, 1), size = (680, 430), link = :x, left_margin = 4mm)
    for sp in 1:2
        plot!(plt, Shape([1, 2, 2, 1], [0, 0, 1.08, 1.08]); color = col("#eef4fd"), linecolor = nothing, subplot = sp)
    end
    plot_basis!(plt, kv, p; sp = 1)
    bern = [0, 0, 0, 1, 1, 1.0]
    for e in 1:3
        xs = range(e - 1, e, length = 60)
        for a in 1:3
            plot!(plt, xs, [bspline(bern, p, a, x - (e - 1)) for x in xs]; color = [INK2, MUTED, INK2][a], lw = 2, subplot = 2)
        end
        annotate!(plt, e - 0.5, 0.62, label("B₂"; color = INK2), subplot = 2)
        annotate!(plt, e - 0.87, 0.93, label("B₁"; color = INK2), subplot = 2)
        annotate!(plt, e - 0.13, 0.93, label("B₃"; color = INK2), subplot = 2)
    end
    plot!(plt; subplot = 1, ylims = (0, 1.15), yticks = [0, 1], xticks = knotticks(kv), ylabel = "B-splines")
    plot!(plt; subplot = 2, ylims = (0, 1.15), yticks = [0, 1], xticks = knotticks(kv), ylabel = "Bernstein polynomials", xlabel = "ξ")
    save(plt, "nurbs_bezier")
end

function nurbs_figures()
    fig_cox_de_boor()
    fig_open_basis()
    fig_multiplicity()
    fig_local_support()
    fig_curve()
    fig_local_control()
    fig_knot_insertion()
    fig_refinement_basis()
    fig_circle()
    fig_tensor_basis()
    fig_tensor_surface()
    fig_spaces()
    fig_global_refinement()
    fig_bezier()
end

include("generate_figures_tsplines.jl")

function generate()
    setup_style()
    nurbs_figures()
    tspline_figures()
    @info "Tutorial figures written to $FIGDIR"
end

end # module

TutorialFigures.generate()
