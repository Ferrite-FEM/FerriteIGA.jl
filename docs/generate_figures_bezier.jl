# Figures of the Bézier extraction page (included from generate_figures.jl)

const BEZIER_KV = [0, 0, 0, 1, 2, 3, 3, 3.0]
const BEZIER_POINTS = [(0.0, 0.0), (0.5, 1.5), (2.0, 2.0), (3.5, 0.5), (4.5, 1.5)]

# The Bernstein polynomial number a (1-based) of degree p on [0, 1]
bernstein(p, a, t) = binomial(p, a - 1) * t^(a - 1) * (1 - t)^(p - a + 1)

function fig_bezier_bernstein()
    plt = plot(layout = (1, 3), size = (900, 270), bottom_margin = 6mm, left_margin = 4mm)
    ts = range(0, 1, length = 100)
    for (sp, p) in enumerate(1:3)
        for a in 1:p+1
            ys = [bernstein(p, a, t) for t in ts]
            plot!(plt, ts, ys; color = SERIES[a], subplot = sp)
            k = argmax(ys)
            x = a == 1 ? 0.07 : a == p + 1 ? 0.93 : ts[k]
            y = (a == 1 || a == p + 1) ? 0.93 : ys[k] + 0.08
            annotate!(plt, x, y, label("B" * sub(a)), subplot = sp)
        end
        plot!(plt; subplot = sp, title = "p = $p:  $(p + 1) polynomials", ylims = (0, 1.08), yticks = [0, 1], xticks = [0, 1], xlabel = "position in the element")
    end
    save(plt, "bezier_bernstein")
end

# The pieces of the B-splines on the middle element, written in the Bernstein polynomials of that element
function fig_bezier_decomposition()
    kv, p, e = BEZIER_KV, 2, 2
    C, _ = FerriteIGA._compute_bezier_extraction_operators(p, kv)
    Ce = Matrix(C[e])
    functions = e:e+p #The B-splines which are non-zero on element e
    xs = range(e - 1, e, length = 80)
    frac(c) = c == 1 ? "" : c == 0.5 ? "½ · " : string(c) * " · "
    plt = plot(layout = (1, 3), size = (900, 290), bottom_margin = 6mm, left_margin = 4mm)
    for (sp, (row, A)) in enumerate(pairs(functions))
        xall = samples(kv)
        plot!(plt, xall, [bspline(kv, p, A, x) for x in xall]; color = MUTED, lw = 1.5, subplot = sp)
        plot!(plt, Shape([e - 1, e, e, e - 1], [0, 0, 1.08, 1.08]); color = col("#eef4fd"), linecolor = nothing, subplot = sp)
        terms = String[]
        for a in 1:p+1
            c = Ce[row, a]
            c == 0 && continue
            plot!(plt, xs, [c * bernstein(p, a, x - (e - 1)) for x in xs]; color = SERIES[a], lw = 1.5, subplot = sp)
            push!(terms, frac(c) * "B" * sub(a))
        end
        plot!(plt, xs, [bspline(kv, p, A, x) for x in xs]; color = INK, lw = 2.5, subplot = sp)
        plot!(plt; subplot = sp, title = "N" * sub(A) * " = " * join(terms, " + "),
            ylims = (0, 1.08), yticks = [0, 1], xticks = knotticks(kv), xlabel = "ξ")
    end
    save(plt, "bezier_decomposition")
end

# Repeating the interior knots turns the B-spline basis into the Bernstein basis of each element
function fig_bezier_knot_insertion()
    p = 2
    kvs = (BEZIER_KV, sort([BEZIER_KV; 1.0]), sort([BEZIER_KV; 1.0; 2.0]))
    names = ("original knot vector", "knot 1 inserted", "knots 1 and 2 inserted")
    plt = plot(layout = (3, 1), size = (680, 560), link = :x, left_margin = 4mm)
    for (sp, kv) in enumerate(kvs)
        plot_basis!(plt, kv, p; sp, labels = false)
        plot!(plt; subplot = sp, ylims = (0, 1.08), yticks = [0, 1], xticks = knotticks(kv),
            title = names[sp] * ":  {" * join(fmt.(kv), ", ") * "},  $(nfunctions(kv, p)) functions", xlabel = sp == 3 ? "ξ" : "")
    end
    save(plt, "bezier_knot_insertion")
end

function fig_bezier_curve()
    kv, p, P = BEZIER_KV, 2, BEZIER_POINTS
    mesh = NURBSMesh((copy(kv),), (p,), [Vec(pt) for pt in P])
    refined = knotinsertion(knotinsertion(mesh, 1.0; dir = 1), 2.0; dir = 1)
    Pb = [Tuple(x) for x in refined.control_points]
    plt = plot(layout = (1, 2), size = (900, 330))
    pts = [curve_point(kv, p, P, x) for x in samples(kv)]
    for sp in 1:2
        plot!(plt, first.(pts), last.(pts); color = sp == 1 ? BLUE : MUTED, lw = 2.5, subplot = sp)
    end
    plot_polygon!(plt, P; sp = 1, offset = (0.0, 0.2))
    #Each element has p+1 Bézier control points; neighbouring elements share one
    plot_polygon!(plt, P; sp = 2, labels = false, color = FAINT)
    for e in 1:3
        ids = (e-1)*p+1:e*p+1
        xs = range(e - 1, e, length = 50)
        epts = [curve_point(kv, p, P, x) for x in xs]
        plot!(plt, first.(epts), last.(epts); color = SERIES[e], lw = 2.5, subplot = 2)
        plot!(plt, first.(Pb[ids]), last.(Pb[ids]); color = SERIES[e], lw = 1, subplot = 2)
        scatter!(plt, first.(Pb[ids]), last.(Pb[ids]); color = SURFACE, markerstrokecolor = SERIES[e], markerstrokewidth = 1.5, markersize = 5, subplot = 2)
    end
    for sp in 1:2, v in unique(kv)
        x, y = curve_point(kv, p, P, v)
        scatter!(plt, [x], [y]; color = INK, markerstrokecolor = SURFACE, markersize = 5, subplot = sp)
    end
    for sp in 1:2
        plot!(plt; subplot = sp, aspect_ratio = :equal, framestyle = :none, xlims = (-0.3, 4.8), ylims = (-0.3, 2.4),
            title = ["the spline: 5 control points, shared by the elements", "the same curve: 3 Bézier control points per element"][sp])
    end
    save(plt, "bezier_curve")
end

function fig_bezier_surface()
    mesh = quarter_annulus()
    grid = BezierGrid(mesh)
    kv1, kv2 = mesh.knot_vectors
    ordering = FerriteIGA._bernstein_ordering(FerriteIGA.Bernstein{RefQuadrilateral,2}())
    plt = plot(layout = (1, 2), size = (900, 430))
    ts = range(0, 1, length = 60)
    θ = range(0, π / 2, length = 60)
    for sp in 1:2
        plot!(plt, Shape([cos.(θ); 2 .* cos.(reverse(θ))], [sin.(θ); 2 .* sin.(reverse(θ))]); color = TINT, linecolor = nothing, subplot = sp)
        for v in unique(kv1)
            pts = [nurbs_point(mesh, v, t) for t in ts]
            plot!(plt, first.(pts), last.(pts); color = INK, lw = 1, subplot = sp)
        end
        for v in unique(kv2)
            pts = [nurbs_point(mesh, t, v) for t in ts]
            plot!(plt, first.(pts), last.(pts); color = INK, lw = 1, subplot = sp)
        end
    end
    #Left: the control net of the spline
    m1, m2 = length(kv1) - 3, length(kv2) - 3
    cp(i, j) = mesh.control_points[(j-1)*m1+i]
    segs = [[((cp(i, j)[1], cp(i, j)[2]), (cp(i + 1, j)[1], cp(i + 1, j)[2])) for j in 1:m2 for i in 1:m1-1];
            [((cp(i, j)[1], cp(i, j)[2]), (cp(i, j + 1)[1], cp(i, j + 1)[2])) for j in 1:m2-1 for i in 1:m1]]
    segments!(plt, segs; sp = 1, color = ORANGE, lw = 1)
    scatter!(plt, first.(mesh.control_points), last.(mesh.control_points); subplot = 1, color = SURFACE, markerstrokecolor = ORANGE, markerstrokewidth = 1.5, markersize = 4)
    #Right: the 3 × 3 Bézier control points of each element
    for cellid in 1:getncells(grid)
        xb, _ = get_bezier_coordinates(grid, cellid)
        X = Matrix{Tuple{Float64,Float64}}(undef, 3, 3)
        for (k, tensor_index) in enumerate(ordering)
            X[tensor_index] = (xb[k][1], xb[k][2])
        end
        esegs = [[(X[i, j], X[i+1, j]) for j in 1:3 for i in 1:2]; [(X[i, j], X[i, j+1]) for j in 1:2 for i in 1:3]]
        segments!(plt, esegs; sp = 2, color = ORANGE, lw = 1)
        scatter!(plt, first.(vec(X)), last.(vec(X)); subplot = 2, color = SURFACE, markerstrokecolor = ORANGE, markerstrokewidth = 1.5, markersize = 4)
    end
    for sp in 1:2
        plot!(plt; subplot = sp, aspect_ratio = :equal, framestyle = :none, xlims = (-0.1, 2.15), ylims = (-0.1, 2.15),
            title = ["NURBS control net:  $(length(mesh.control_points)) control points for the patch", "Bézier control net:  9 control points for each of the $(getncells(grid)) elements"][sp])
    end
    save(plt, "bezier_surface")
end

function bezier_figures()
    fig_bezier_bernstein()
    fig_bezier_decomposition()
    fig_bezier_knot_insertion()
    fig_bezier_curve()
    fig_bezier_surface()
end
