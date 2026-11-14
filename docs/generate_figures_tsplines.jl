# Figures of the T-spline tutorial (included from generate_figures.jl)

# # #
# T-meshes used in the tutorial
# # #

# A T-mesh on [0,n]² with full lines at the integers, and additional vertical segments (x, y0, y1) and
# horizontal segments (y, x0, x1). The boundary knots are repeated, such that the boundary is open.
function tmesh_with_segments(vsegs, hsegs; p, n = 6)
    r = p ÷ 2
    line_values(segs) = (vals = sort(unique(Float64[0:n; [s[1] for s in segs]])); [fill(0.0, r); vals; fill(Float64(n), r)])
    k1, k2 = line_values(vsegs), line_values(hsegs)
    n1, n2 = length(k1), length(k2)
    vedge(i, j) = isinteger(k1[i]) || any(s -> s[1] == k1[i] && s[2] <= k2[j] && k2[j+1] <= s[3], vsegs)
    hedge(i, j) = isinteger(k2[j]) || any(s -> s[1] == k2[j] && s[2] <= k1[i] && k1[i+1] <= s[3], hsegs)
    owner = zeros(Int, n1 - 1, n2 - 1)
    cells = NTuple{4,Int}[]
    for j in 1:n2-1, i in 1:n1-1
        owner[i, j] != 0 && continue
        i1 = i + 1
        while i1 < n1 && !vedge(i1, j); i1 += 1; end
        j1 = j + 1
        while j1 < n2 && !hedge(i, j1); j1 += 1; end
        push!(cells, (i, i1, j, j1))
        owner[i:i1-1, j:j1-1] .= length(cells)
    end
    return TSplineMesh((k1, k2), p, cells)
end

# One vertical line ends at η = 2 and one horizontal line ends at ξ = 4: analysis-suitable for cubics
mesh_good(p) = tmesh_with_segments([(2.5, 0.0, 2.0)], [(4.5, 4.0, 6.0)]; p)
# The same two lines, one element longer: their T-junction extensions cross for cubics
mesh_bad(p) = tmesh_with_segments([(2.5, 0.0, 3.0)], [(4.5, 3.0, 6.0)]; p)

# # #
# Drawing helpers
# # #

function segments!(plt, segs; sp = 1, label = "", kwargs...)
    isempty(segs) && return plt
    xs, ys = Float64[], Float64[]
    for ((x1, y1), (x2, y2)) in segs
        append!(xs, (x1, x2, NaN))
        append!(ys, (y1, y2, NaN))
    end
    plot!(plt, xs, ys; subplot = sp, label, kwargs...)
end

index_edges(mesh::TSplineMesh{2}) = [[((i, j), (i, j + 1)) for (i, j) in mesh.walls[1]]; [((i, j), (i + 1, j)) for (i, j) in mesh.walls[2]]]

function parametric_edges(mesh::TSplineMesh{2})
    k1, k2 = mesh.knots
    segs = [((k1[a[1]], k2[a[2]]), (k1[b[1]], k2[b[2]])) for (a, b) in index_edges(mesh)]
    return unique!(filter!(s -> s[1] != s[2], segs))
end

fmt(v) = isinteger(v) ? string(Int(v)) : string(v)

# Axes of the index space: the index of each line, and its knot value on a second row
function index_axes!(plt, mesh; sp = 1, values = true, kwargs...)
    k1, k2 = mesh.knots
    n1, n2 = length(k1), length(k2)
    xt = values ? ["$i\n$(fmt(k1[i]))" for i in 1:n1] : string.(1:n1)
    yt = values ? ["$j  ($(fmt(k2[j])))" for j in 1:n2] : string.(1:n2)
    plot!(plt; subplot = sp, aspect_ratio = :equal, xlims = (0.5, n1 + 0.5), ylims = (0.5, n2 + 0.5),
        xticks = (1:n1, xt), yticks = (1:n2, yt), kwargs...)
end

draw_index_mesh!(plt, mesh; sp = 1, color = INK, lw = 1) = segments!(plt, index_edges(mesh); sp, color, lw)
draw_parametric_mesh!(plt, mesh; sp = 1, color = INK, lw = 1) = segments!(plt, parametric_edges(mesh); sp, color, lw)

function anchor_id(mesh, X)
    A = findfirst(==(Tuple(round.(Int, 2 .* X))), mesh.anchors)
    A === nothing && error("The T-mesh has no anchor at the index coordinates $X")
    return A
end
anchor_position(mesh, A) = mesh.anchors[A] ./ 2

function draw_anchors!(plt, mesh; sp = 1, color = MUTED, markersize = 3, parametric = false)
    pts = [anchor_position(mesh, A) for A in eachindex(mesh.anchors)]
    if parametric
        pts = greville_points(mesh)
    end
    scatter!(plt, first.(pts), last.(pts); subplot = sp, color, markerstrokecolor = color, markersize)
end

# The T-junctions (i, j, direction of the missing edge) of a 2D T-mesh
function tjunctions(mesh::TSplineMesh{2})
    vedges, hedges = mesh.walls
    n1, n2 = length.(mesh.knots)
    out = Tuple{Int,Int,Symbol}[]
    for i in 2:n1-1, j in 2:n2-1
        left, right, down, up = (i - 1, j) in hedges, (i, j) in hedges, (i, j - 1) in vedges, (i, j) in vedges
        left + right + down + up == 3 || continue
        push!(out, (i, j, !left ? :left : !right ? :right : !down ? :down : :up))
    end
    return out
end

function draw_tjunctions!(plt, mesh; sp = 1, markersize = 6)
    T = tjunctions(mesh)
    scatter!(plt, [t[1] for t in T], [t[2] for t in T]; subplot = sp, color = SURFACE, markerstrokecolor = ORANGE, markerstrokewidth = 2, markersize)
end

# The local index vectors of the anchor A: the marching paths, the lines that are met (ticks), and the support
function draw_local_index_vectors!(plt, mesh, A; sp = 1, color = BLUE, tick = 0.16)
    X = anchor_position(mesh, A)
    I1, I2 = mesh.local_indices[A]
    plot!(plt, Shape([I1[1], I1[end], I1[end], I1[1]], [I2[1], I2[1], I2[end], I2[end]]); subplot = sp, color, fillalpha = 0.13, linecolor = nothing)
    segs = [((I1[1], X[2]), (I1[end], X[2])), ((X[1], I2[1]), (X[1], I2[end]))]
    append!(segs, [((i, X[2] - tick), (i, X[2] + tick)) for i in unique(I1)])
    append!(segs, [((X[1] - tick, j), (X[1] + tick, j)) for j in unique(I2)])
    segments!(plt, segs; sp, color, lw = 3)
    scatter!(plt, [X[1]], [X[2]]; subplot = sp, color, markerstrokecolor = SURFACE, markerstrokewidth = 1.5, markersize = 8)
end

index_vector_string(I) = "(" * join(I, ", ") * ")"

tspline_value(mesh, A, x, y) = bspline(mesh.local_knots[A][1], x) * bspline(mesh.local_knots[A][2], y)

# # #
# Figures
# # #

# The relation between the knot vector of a B-spline and the index lines and anchors of a T-spline, in 1D
function fig_ts_anchors_1d()
    plt = plot(layout = grid(2, 2, heights = [0.62, 0.38]), size = (900, 430), left_margin = 4mm)
    for (c, p) in enumerate((3, 2))
        kv = p == 3 ? [0, 0, 0, 0, 1, 2, 3, 3, 3, 3.0] : [0, 0, 0, 1, 2, 3, 3, 3.0]
        s = cld(p, 2)
        lines = kv[s+1:end-s] #The index lines: the knot vector without the outermost knots
        n = nfunctions(kv, p)
        top, bottom = c, c + 2
        plot_basis!(plt, kv, p; sp = top)
        plot!(plt; subplot = top, ylims = (0, 1.15), yticks = [0, 1], xticks = knotticks(kv), xlabel = "ξ  (parameter space)",
            title = "p = $p:  knot vector {" * join(fmt.(kv), ", ") * "}")
        #Index space
        m = length(lines)
        segments!(plt, [((i, -0.3), (i, 0.3)) for i in 1:m]; sp = bottom, color = INK, lw = 1.5)
        segments!(plt, [((1, 0.0), (m, 0.0))]; sp = bottom, color = MUTED, lw = 1)
        anchors = isodd(p) ? collect(1.0:m) : collect(1.5:1:m-0.5)
        @assert length(anchors) == n
        for (A, x) in enumerate(anchors)
            scatter!(plt, [x], [0.0]; subplot = bottom, color = SERIES[A], markerstrokecolor = SURFACE, markersize = 8)
            annotate!(plt, x, 0.55, label("N" * sub(A)), subplot = bottom)
        end
        #Local index vector of the third function
        A = 3
        k = isodd(p) ? (p + 1) ÷ 2 : p ÷ 2 + 1
        lo = isodd(p) ? anchors[A] - k : anchors[A] - k + 0.5
        hi = isodd(p) ? anchors[A] + k : anchors[A] + k - 0.5
        lo, hi = max(lo, 1), min(hi, m)
        segments!(plt, [((lo, -0.62), (hi, -0.62)); [((i, -0.74), (i, -0.5)) for i in lo:hi]]; sp = bottom, color = SERIES[A], lw = 2.5)
        annotate!(plt, (lo + hi) / 2, -0.98, label("lines used by N₃"; color = INK2), subplot = bottom)
        plot!(plt; subplot = bottom, ylims = (-1.2, 0.8), xlims = (0.4, m + 0.6), yticks = false, yaxis = false,
            xticks = (1:m, ["$i\n$(fmt(lines[i]))" for i in 1:m]), xlabel = "index of the line, and its knot value  (index space)", bottom_margin = 8mm)
    end
    save(plt, "tspline_anchors_1d")
end

# A NURBS patch seen as a T-spline: tensor product mesh in index space
function fig_ts_index_nurbs()
    k = [0, 0, 1, 2, 3, 4, 4.0]
    cells = [(i, i + 1, j, j + 1) for j in 1:6 for i in 1:6]
    plt = plot(layout = (1, 2), size = (900, 460))
    for (sp, (p, X)) in enumerate(((3, (4, 4)), (2, (3.5, 3.5))))
        mesh = TSplineMesh((k, copy(k)), p, cells)
        A = anchor_id(mesh, X)
        draw_index_mesh!(plt, mesh; sp)
        draw_anchors!(plt, mesh; sp)
        draw_local_index_vectors!(plt, mesh, A; sp)
        I1, I2 = mesh.local_indices[A]
        index_axes!(plt, mesh; sp, xlabel = "index i  (second row: knot value ξ)", ylabel = "index j  (knot value η)",
            title = (p == 3 ? "p = 3: one anchor at every vertex" : "p = 2: one anchor in every cell") * "\nmarked anchor: i = $(index_vector_string(I1)),  j = $(index_vector_string(I2))",
            bottom_margin = 8mm, left_margin = 5mm)
    end
    save(plt, "tspline_index_nurbs")
end

function fig_ts_tmesh()
    mesh = mesh_good(3)
    plt = plot(layout = grid(1, 2, widths = [0.55, 0.45]), size = (940, 500))
    draw_index_mesh!(plt, mesh; sp = 1)
    draw_anchors!(plt, mesh; sp = 1)
    draw_tjunctions!(plt, mesh; sp = 1, markersize = 8)
    index_axes!(plt, mesh; sp = 1, xlabel = "index i  (second row: knot value ξ)", ylabel = "index j  (knot value η)",
        title = "index space:  $(length(mesh.cells)) cells,  $(getnnodes(mesh)) anchors,  T-junctions circled", bottom_margin = 8mm, left_margin = 5mm)
    draw_parametric_mesh!(plt, mesh; sp = 2)
    plot!(plt; subplot = 2, aspect_ratio = :equal, xlims = (-0.2, 6.2), ylims = (-0.2, 6.2), xticks = 0:6, yticks = 0:6,
        xlabel = "ξ", ylabel = "η", title = "parameter space")
    save(plt, "tspline_tmesh")
end

function fig_ts_rays_odd()
    mesh = mesh_good(3)
    anchors = ((8, 3), (3, 6), (5, 4))
    names = ("far from the T-junctions", "its paths pass the ends of two lines", "the T-junction itself")
    plt = plot(layout = (1, 3), size = (1000, 400))
    for (sp, X) in enumerate(anchors)
        A = anchor_id(mesh, X)
        draw_index_mesh!(plt, mesh; sp)
        draw_local_index_vectors!(plt, mesh, A; sp, color = SERIES[sp])
        I1, I2 = mesh.local_indices[A]
        index_axes!(plt, mesh; sp, values = false, xlabel = "i", ylabel = sp == 1 ? "j" : "",
            title = "anchor $(X): $(names[sp])\ni = $(index_vector_string(I1)),  j = $(index_vector_string(I2))")
    end
    save(plt, "tspline_rays_odd")
end

function fig_ts_rays_even()
    mesh = mesh_good(2)
    anchors = ((8.5, 2.5), (5.0, 4.5), (2.5, 7.0))
    plt = plot(layout = (1, 3), size = (1000, 400))
    for (sp, X) in enumerate(anchors)
        A = anchor_id(mesh, X)
        draw_index_mesh!(plt, mesh; sp)
        draw_local_index_vectors!(plt, mesh, A; sp, color = SERIES[sp])
        I1, I2 = mesh.local_indices[A]
        index_axes!(plt, mesh; sp, values = false, xlabel = "i", ylabel = sp == 1 ? "j" : "",
            title = "anchor $(X)\ni = $(index_vector_string(I1)),  j = $(index_vector_string(I2))")
    end
    save(plt, "tspline_rays_even")
end

function fig_ts_function()
    mesh = mesh_good(3)
    xs = collect(range(0, prevfloat(6.0), length = 241))
    plt = plot(layout = (1, 2), size = (900, 450))
    for (sp, X) in enumerate(((3, 6), (5, 4)))
        A = anchor_id(mesh, X)
        Z = [tspline_value(mesh, A, x, y) for y in xs, x in xs]
        heatmap!(plt, xs, xs, Z; subplot = sp, color = SEQ, colorbar = false)
        draw_parametric_mesh!(plt, mesh; sp, lw = 0.7)
        K1, K2 = mesh.local_knots[A]
        segments!(plt, [((K1[1], K2[1]), (K1[end], K2[1])), ((K1[end], K2[1]), (K1[end], K2[end])), ((K1[end], K2[end]), (K1[1], K2[end])), ((K1[1], K2[end]), (K1[1], K2[1]))];
            sp, color = SERIES[sp+1], lw = 2.5)
        plot!(plt; subplot = sp, aspect_ratio = :equal, xlims = (0, 6), ylims = (0, 6), xticks = 0:6, yticks = 0:6, xlabel = "ξ", ylabel = "η",
            title = "basis function of anchor $(X)\nξ-knots {" * join(fmt.(K1), ", ") * "},  η-knots {" * join(fmt.(K2), ", ") * "}")
    end
    save(plt, "tspline_function")
end

element_edges(mesh) = reduce(vcat, [[((a, c), (b, c)), ((b, c), (b, d)), ((b, d), (a, d)), ((a, d), (a, c))] for ((a, b), (c, d)) in mesh.elements])

function fig_ts_elements()
    plt = plot(layout = (1, 2), size = (900, 460))
    for (sp, p) in enumerate((3, 2))
        mesh = mesh_good(p)
        #All element boundaries; the ones which are edges of the T-mesh are then covered by the T-mesh
        segments!(plt, element_edges(mesh); sp, color = ORANGE, lw = 2, ls = :dash)
        draw_parametric_mesh!(plt, mesh; sp, lw = 1.5)
        ncells = count(c -> all(d -> mesh.knots[d][c[d][1]] < mesh.knots[d][c[d][2]], 1:2), mesh.cells)
        plot!(plt; subplot = sp, aspect_ratio = :equal, xlims = (-0.2, 6.2), ylims = (-0.2, 6.2), xticks = 0:6, yticks = 0:6, xlabel = "ξ", ylabel = "η",
            title = "p = $p:  $ncells cells of the T-mesh,  $(getncells(mesh)) Bézier elements")
    end
    save(plt, "tspline_elements")
end

# Split the extension of a T-junction into the face extension and the edge extension
function extension_parts(mesh, (i, j, missing))
    hext, vext = tjunction_extensions(mesh)
    if missing in (:left, :right)
        (i1, i2, _) = only(filter(e -> e[3] == j && e[1] <= i <= e[2], hext))
        face = missing == :left ? ((i1, j), (i, j)) : ((i, j), (i2, j))
        edge = missing == :left ? ((i, j), (i2, j)) : ((i1, j), (i, j))
    else
        (_, j1, j2) = only(filter(e -> e[1] == i && e[2] <= j <= e[3], vext))
        face = missing == :down ? ((i, j1), (i, j)) : ((i, j), (i, j2))
        edge = missing == :down ? ((i, j), (i, j2)) : ((i, j1), (i, j))
    end
    return face, edge
end

function draw_extensions!(plt, mesh; sp = 1, face_color = BLUE, edge_color = ORANGE, lw = 4)
    for T in tjunctions(mesh)
        face, edge = extension_parts(mesh, T)
        segments!(plt, [edge]; sp, color = edge_color, lw)
        segments!(plt, [face]; sp, color = face_color, lw)
    end
    draw_tjunctions!(plt, mesh; sp, markersize = 7)
end

function fig_ts_extensions()
    plt = plot(layout = (1, 2), size = (900, 480))
    for (sp, p) in enumerate((3, 2))
        mesh = tmesh_with_segments([(3.5, 0.0, 3.0)], NTuple{3,Float64}[]; p)
        n1, n2 = length.(mesh.knots)
        draw_index_mesh!(plt, mesh; sp)
        T = only(tjunctions(mesh))
        face, edge = extension_parts(mesh, T)
        nface, nedge = abs(face[2][2] - face[1][2]), abs(edge[2][2] - edge[1][2])
        bays(n) = n == 1 ? "1 bay" : "$n bays"
        segments!(plt, [face]; sp, color = BLUE, lw = 5)
        segments!(plt, [edge]; sp, color = ORANGE, lw = 5)
        scatter!(plt, [T[1]], [T[2]]; subplot = sp, color = SURFACE, markerstrokecolor = INK, markerstrokewidth = 1.5, markersize = 7)
        #Legend below the mesh
        segments!(plt, [((1.0, 0.2), (1.7, 0.2))]; sp, color = BLUE, lw = 5)
        annotate!(plt, 1.9, 0.2, label("face extension: $(bays(nface))"; align = :left), subplot = sp)
        segments!(plt, [((5.6, 0.2), (6.3, 0.2))]; sp, color = ORANGE, lw = 5)
        annotate!(plt, 6.5, 0.2, label("edge extension: $(bays(nedge))"; align = :left), subplot = sp)
        scatter!(plt, [1.35], [-0.45]; subplot = sp, color = SURFACE, markerstrokecolor = INK, markerstrokewidth = 1.5, markersize = 7)
        annotate!(plt, 1.9, -0.45, label("T-junction: the line coming from below ends here"; align = :left), subplot = sp)
        plot!(plt; subplot = sp, aspect_ratio = :equal, framestyle = :none, xlims = (0.6, n1 + 0.4), ylims = (-0.9, n2 + 0.4), title = "p = $p")
    end
    save(plt, "tspline_extensions")
end

function crossing_points(mesh)
    hext, vext = tjunction_extensions(mesh)
    return [(i, j) for (i1, i2, j) in hext for (i, j1, j2) in vext if i1 <= i <= i2 && j1 <= j <= j2]
end

function fig_ts_as()
    plt = plot(layout = (1, 2), size = (940, 500))
    for (sp, mesh) in enumerate((mesh_good(3), mesh_bad(3)))
        draw_index_mesh!(plt, mesh; sp)
        #Horizontal extensions in blue and vertical extensions in orange
        hext, vext = tjunction_extensions(mesh)
        segments!(plt, [((i1, j), (i2, j)) for (i1, i2, j) in hext]; sp, color = BLUE, lw = 4)
        segments!(plt, [((i, j1), (i, j2)) for (i, j1, j2) in vext]; sp, color = ORANGE, lw = 4)
        T = tjunctions(mesh)
        scatter!(plt, [t[1] for t in T], [t[2] for t in T]; subplot = sp, color = SURFACE, markerstrokecolor = INK, markerstrokewidth = 1.5, markersize = 6)
        X = crossing_points(mesh)
        if !isempty(X)
            scatter!(plt, first.(X), last.(X); subplot = sp, color = INK, markershape = :xcross, markersize = 9, markerstrokewidth = 3)
            annotate!(plt, X[1][1] + 0.25, X[1][2] + 0.35, label("the extensions cross"; align = :left), subplot = sp)
        end
        @assert isempty(X) == is_analysis_suitable(mesh)
        index_axes!(plt, mesh; sp, values = false, xlabel = "i", ylabel = "j",
            title = is_analysis_suitable(mesh) ? "analysis-suitable" : "not analysis-suitable")
    end
    save(plt, "tspline_as")
end

function fig_ts_pou()
    xs = collect(range(0, prevfloat(6.0), length = 181))
    plt = plot(layout = grid(1, 2, widths = [0.435, 0.565]), size = (960, 450))
    for (sp, mesh) in enumerate((mesh_good(3), mesh_bad(3)))
        Z = [sum(tspline_value(mesh, A, x, y) for A in 1:getnnodes(mesh)) - 1 for y in xs, x in xs]
        heatmap!(plt, xs, xs, Z; subplot = sp, color = DIV, clims = (-0.05, 0.05), colorbar = sp == 2, colorbar_title = sp == 2 ? "sum of all basis functions − 1" : "")
        segments!(plt, element_edges(mesh); sp, color = INK, lw = 0.6)
        draw_parametric_mesh!(plt, mesh; sp, lw = 1.2)
        nb = (mesh.order + 1)^2
        for (e, ((a, b), (c, d))) in enumerate(mesh.elements)
            n = length(mesh.IEN[e])
            n == nb || annotate!(plt, (a + b) / 2, (c + d) / 2, label(string(n); size = 7), subplot = sp)
        end
        lo, hi = extrema(length, mesh.IEN)
        plot!(plt; subplot = sp, aspect_ratio = :equal, xlims = (0, 6), ylims = (0, 6), xticks = 0:6, yticks = 0:6, xlabel = "ξ", ylabel = sp == 1 ? "η" : "",
            title = (is_analysis_suitable(mesh) ? "analysis-suitable" : "not analysis-suitable") *
                    (lo == hi ? ":  $nb functions on every element" : ":  $lo to $hi functions per element\n(the number is written in elements with more than $nb)"))
    end
    save(plt, "tspline_pou")
end

function fig_ts_graded()
    mesh = generate_tspline_patch(:graded_rectangle, [2, 4, 8], [1, 2, 4], 3; size = (2.0, 1.0))
    plt = plot(size = (900, 470), aspect_ratio = :equal, xlims = (-0.05, 2.05), ylims = (-0.05, 1.05), xlabel = "x", ylabel = "y",
        title = "generate_tspline_patch(:graded_rectangle, [2, 4, 8], [1, 2, 4], 3; size = (2.0, 1.0)):  " *
                "$(getncells(mesh)) Bézier elements,  $(getnnodes(mesh)) control points")
    segments!(plt, element_edges(mesh); color = ORANGE, lw = 1.5, ls = :dash)
    draw_parametric_mesh!(plt, mesh; lw = 1.5)
    g = mesh.control_points
    scatter!(plt, first.(g), last.(g); color = BLUE, markerstrokecolor = SURFACE, markersize = 5)
    save(plt, "tspline_graded")
end

function tspline_figures()
    fig_ts_anchors_1d()
    fig_ts_index_nurbs()
    fig_ts_tmesh()
    fig_ts_rays_odd()
    fig_ts_rays_even()
    fig_ts_function()
    fig_ts_elements()
    fig_ts_extensions()
    fig_ts_as()
    fig_ts_pou()
    fig_ts_graded()
end
