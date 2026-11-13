#!/usr/bin/env julia

# Multipatch IGA version of the Ferrite.jl logo.
#
# The polygonal faces are taken from the Ferrite.jl logo (`Ferrite.jl/docs/logo.jl`, which
# generates them with neper). Each face is meshed with one or more B-spline patches:
# four-sided faces become a single patch, while the five-sided faces are split into two
# four-sided subpatches by connecting one corner to the midpoint of the opposite edge
# (`pentagon = :split`), or meshed with a single patch with a kinked side (`pentagon = :kink`).
# The patches are not connected topologically (each patch has its own control points).
#
# The logo uses the `:kink` variant. Run from the repository root:
#   julia --project=docs docs/logo.jl
# to regenerate `docs/src/assets/logo.svg` (used in the docs sidebar) and
# `docs/src/assets/logo-name.svg` (logo with "FerriteIGA.jl" below, as Ferrite.jl's logo-name.svg).

using Ferrite, FerriteIGA

# Julia colors (same as the Ferrite.jl logo)
const julia_blue   = "rgb(25.1%,38.8%,84.7%)"
const julia_green  = "rgb(22%,59.6%,14.9%)"
const julia_purple = "rgb(58.4%,34.5%,69.8%)"
const julia_red    = "rgb(79.6%,23.5%,20%)"

# Faces of the Ferrite.jl logo in the unit square, corners in counter-clockwise order.
const FACES = [
    (julia_blue,   [(0.0, 0.0), (0.548, 0.0), (0.489, 0.224), (0.172, 0.368), (0.0, 0.325)]),
    (julia_red,    [(0.489, 0.224), (0.548, 0.0), (1.0, 0.0), (1.0, 0.38), (0.801, 0.455)]),
    (julia_green,  [(0.172, 0.368), (0.489, 0.224), (0.801, 0.455), (0.601, 0.775), (0.393, 0.672)]),
    (julia_purple, [(0.0, 0.325), (0.172, 0.368), (0.393, 0.672), (0.0, 1.0)]),
    (julia_purple, [(0.601, 0.775), (0.801, 0.455), (1.0, 0.38), (1.0, 1.0), (0.656, 1.0)]),
    (julia_red,    [(0.0, 1.0), (0.393, 0.672), (0.601, 0.775), (0.656, 1.0)]),
]

# Quality of a (convex) quad: smallest sine of its interior angles.
function quad_quality(q)
    minimum(1:4) do i
        a, b, c = q[mod1(i - 1, 4)], q[i], q[mod1(i + 1, 4)]
        u, v = a - b, c - b
        abs(u[1] * v[2] - u[2] * v[1]) / (hypot(u...) * hypot(v...))
    end
end

# Split a polygon into quads. A pentagon (P1,...,P5) is split from corner Pi to the midpoint
# M of the opposite edge, giving (Pi, Pi+1, Pi+2, M) and (Pi, M, Pi+3, Pi+4). The corner giving
# the best worst-case quad quality is chosen.
function split_into_quads(poly)
    P = [Vec{2}(p) for p in poly]
    length(P) == 4 && return [P]
    @assert length(P) == 5
    candidates = map(1:5) do i
        p(k) = P[mod1(i + k, 5)]
        M = (p(2) + p(3)) / 2
        [[p(0), p(1), p(2), M], [p(0), M, p(3), p(4)]]
    end
    return argmax(qs -> minimum(quad_quality, qs), candidates)
end

# B-spline patch with straight edges spanning the quad (P1, P2, P3, P4).
# The number of elements in each direction follows the target element size `h`.
function quad_patch(q; h, order)
    edgelen(a, b) = norm(q[b] - q[a])
    nel = (max(1, round(Int, (edgelen(1, 2) + edgelen(4, 3)) / 2h)),
           max(1, round(Int, (edgelen(1, 4) + edgelen(2, 3)) / 2h)))
    # Unit square patch, then map the control points with the bilinear map of the quad.
    # Since the control points sit at the Greville points this reproduces the map exactly.
    sq = generate_nurbs_patch(:hypercube, nel, (order, order); cornerpos = (0.0, 0.0), size = (1.0, 1.0))
    bilinear(ξ, η) = (1 - ξ) * (1 - η) * q[1] + ξ * (1 - η) * q[2] + ξ * η * q[3] + (1 - ξ) * η * q[4]
    cps = [bilinear(x[1], x[2]) for x in sq.control_points]
    return NURBSMesh(sq.knot_vectors, sq.orders, cps, sq.weights)
end

# Single B-spline patch for a pentagon, using the same trick as the `:plate_with_hole` patch:
# one side of the patch passes through a corner Pk of the pentagon via a doubled control
# point, giving a sharp corner in the geometry while the parametrization stays C1 (with a
# singular Jacobian at Pk). The corner with the largest interior angle is used as the kink.
#
#          Pk-1 ----- Pk ----- Pk+1       η = 1:  (Pk+1, Pk, Pk, Pk-1)
#           |          :         |
#          Pk+3 -------:------ Pk+2       η = -1: straight edge Pk+2 → Pk+3
function pentagon_patch(poly; h, order)
    P = [Vec{2}(p) for p in poly]
    @assert length(P) == 5
    p(i) = P[mod1(i, 5)]
    angle(i) = (u = p(i - 1) - p(i); v = p(i + 1) - p(i); acos(clamp(u ⋅ v / (norm(u) * norm(v)), -1, 1)))
    k = argmax(angle, 1:5)

    # Coarse quadratic patch: 4 x 3 control points, knot ξ = 0 at the kink.
    A, B = p(k + 2), p(k + 3)
    bottom = [A, 0.75A + 0.25B, 0.25A + 0.75B, B]   # Greville points of a straight edge
    top    = [p(k + 1), p(k), p(k), p(k - 1)]
    cps = [bottom; (bottom + top) / 2; top]
    knot_vectors = ([-1.0, -1, -1, 0, 1, 1, 1], [-1.0, -1, -1, 1, 1, 1])
    weights = ones(length(cps))

    orders = (2, 2)
    for _ in 3:order, dir in 1:2
        orders = orderelevation!(knot_vectors, orders, cps, weights; dir)
    end

    # Refine: each half of the ξ-direction (on either side of the kink) and the η-direction.
    halflen = (norm(p(k + 1) - p(k)) + norm(p(k) - p(k - 1)) + norm(B - A)) / 4
    heightlen = (norm(p(k + 1) - A) + norm(p(k - 1) - B)) / 2
    nξ = max(1, round(Int, halflen / h))
    nη = max(1, round(Int, heightlen / h))
    for ξ in [range(-1.0, 0.0; length = nξ + 1)[2:end-1]; range(0.0, 1.0; length = nξ + 1)[2:end-1]]
        knotinsertion!(knot_vectors, orders, cps, weights, ξ; dir = 1)
    end
    for η in range(-1.0, 1.0; length = nη + 1)[2:end-1]
        knotinsertion!(knot_vectors, orders, cps, weights, η; dir = 2)
    end
    return NURBSMesh(knot_vectors, orders, cps, weights)
end

# Patches for one face. `pentagon = :split` splits five-sided faces into two quad patches,
# `pentagon = :kink` uses a single patch with a doubled control point (see `pentagon_patch`).
function face_patches(poly; h, order, pentagon)
    if length(poly) == 5 && pentagon === :kink
        return [pentagon_patch(poly; h, order)]
    end
    return [quad_patch(q; h, order) for q in split_into_quads(poly)]
end

# Merge disconnected patches into a single BezierGrid, with one cellset per patch and per face.
function multipatch_grid(patches, patch_face)
    cells = eltype(BezierGrid(first(patches)).cells)[]
    nodes = Node{2, Float64}[]
    weights = Float64[]
    beo = FerriteIGA.Optional{BezierExtractionOperator{Float64}}[]
    cellsets = Dict{String, FerriteIGA.OrderedSet{Int}}()
    for (ip, patch) in enumerate(patches)
        g = BezierGrid(patch)
        noffset, coffset = length(nodes), length(cells)
        append!(cells, [typeof(c)(c.nodes .+ noffset) for c in g.cells])
        append!(nodes, g.nodes)
        append!(weights, g.weights)
        append!(beo, g.beo)
        cellrange = (coffset + 1):length(cells)
        cellsets["patch$ip"] = FerriteIGA.OrderedSet(cellrange)
        union!(get!(cellsets, "face$(patch_face[ip])", FerriteIGA.OrderedSet{Int}()), cellrange)
    end
    return BezierGrid(cells, nodes, weights, beo; cellsets)
end

function ferrite_iga_logo(; h = 0.1, order = 2, pentagon = :kink)
    @assert pentagon in (:split, :kink)
    patches = NURBSMesh{2, 2, Float64}[]
    patch_face = Int[]
    for (iface, (_, poly)) in enumerate(FACES)
        for patch in face_patches(poly; h, order, pentagon)
            push!(patches, patch)
            push!(patch_face, iface)
        end
    end
    return multipatch_grid(patches, patch_face), patches, patch_face
end

# ---------------------------------------------------------------------------------------
# SVG output, mimicking Ferrite.jl's logo.svg and logo-name.svg (PGFPlots/TikZ, scale = 5)
# ---------------------------------------------------------------------------------------

const L = 141.734715   # side of the logo square in pt (5 cm)
const MARGIN = 3.32    # gives the same 148.373 pt canvas as Ferrite.jl's logo.svg

# JuliaMono Regular glyphs (as in Ferrite.jl's logo-name.svg: 20.856 pt advance), in pt with
# the baseline at y = 0 and y pointing down. Extracted with fontTools from JuliaMono-Regular.ttf.
const GLYPH_ADVANCE = 20.856
const GLYPHS = Dict(
    'F' => "M3.441 -25.201H18.527V-22.594H6.5V-13.73H17.971V-11.193H6.5V0H3.441Z",
    'e' => "M2.52 -9.42Q2.52 -12.322 3.606 -14.608Q4.693 -16.893 6.587 -18.179Q8.481 -19.466 10.793 -19.466Q13.018 -19.466 14.703 -18.24Q16.389 -17.015 17.328 -14.599Q18.266 -12.183 18.266 -8.725H5.579Q5.579 -5.649 6.978 -4.076Q8.377 -2.503 11.036 -2.503Q12.774 -2.503 14.252 -3.076Q15.729 -3.65 16.563 -4.31L18.023 -2.016Q17.363 -1.443 16.25 -0.912Q15.138 -0.382 13.774 -0.052Q12.409 0.278 11.071 0.278Q8.551 0.278 6.604 -0.939Q4.658 -2.155 3.589 -4.362Q2.52 -6.57 2.52 -9.42ZM10.723 -16.72Q9.35 -16.72 8.186 -16.05Q7.022 -15.381 6.335 -14.156Q5.649 -12.931 5.649 -11.367H15.486Q15.486 -12.948 14.869 -14.165Q14.252 -15.381 13.165 -16.05Q12.079 -16.72 10.723 -16.72Z",
    'r' => "M6.309 -9.976Q6.309 -12.983 6.222 -14.191Q6.135 -15.399 5.927 -16.511H2.59V-19.118H8.255L8.534 -15.677Q9.333 -17.589 10.802 -18.527Q12.27 -19.466 13.852 -19.466Q15.19 -19.466 16.633 -18.962Q18.075 -18.458 18.857 -17.71L16.98 -15.381Q16.441 -15.955 15.494 -16.441Q14.547 -16.928 13.365 -16.928Q11.888 -16.928 11.002 -16.007Q10.115 -15.086 9.707 -13.183Q9.298 -11.28 9.298 -8.238V0H6.309Z",
    'i' => "M3.632 -19.118H11.957V-2.572H17.658V0H3.198V-2.572H8.916V-16.511H3.632ZM8.082 -25.236Q8.082 -26.487 8.603 -27.009Q9.124 -27.53 10.376 -27.53Q11.627 -27.53 12.149 -27.009Q12.67 -26.487 12.67 -25.236Q12.67 -23.984 12.149 -23.463Q11.627 -22.942 10.376 -22.942Q9.124 -22.942 8.603 -23.463Q8.082 -23.984 8.082 -25.236Z",
    't' => "M8.064 -4.171V-16.511H3.111V-19.118H8.064V-25.201H11.106V-19.118H16.528V-16.511H11.106V-5.214Q11.106 -4.067 11.384 -3.459Q11.662 -2.85 12.322 -2.607Q12.983 -2.364 14.234 -2.364Q14.79 -2.364 15.833 -2.546Q16.876 -2.729 17.745 -2.92V-0.278Q14.738 0.278 13.087 0.278Q10.567 0.278 9.316 -0.869Q8.064 -2.016 8.064 -4.171Z",
    'I' => "M3.893 -2.433H8.933V-22.768H4.554V-25.201H16.268V-22.768H11.992V-2.433H16.928V0H3.893Z",
    'G' => "M2.051 -12.653Q2.051 -18.84 4.076 -22.194Q6.1 -25.549 10.88 -25.549Q13.435 -25.549 15.46 -24.532Q17.484 -23.515 18.458 -21.343L16.059 -19.952Q15.121 -21.69 13.887 -22.299Q12.653 -22.907 11.158 -22.907Q8.638 -22.907 7.308 -21.569Q5.979 -20.23 5.544 -18.04Q5.11 -15.851 5.11 -12.635Q5.11 -9.872 5.483 -7.665Q5.857 -5.457 7.117 -3.884Q8.377 -2.312 10.828 -2.312Q13.678 -2.312 14.869 -4.51Q16.059 -6.709 16.059 -10.706H10.011V-13.139H18.736V0H16.129V-2.677H16.024Q15.381 -1.286 14.191 -0.469Q13 0.348 10.671 0.348Q5.857 0.348 3.954 -3.215Q2.051 -6.778 2.051 -12.653Z",
    'A' => "M8.69 -25.201H12.166L19.361 0H16.059L14.182 -7.056H6.604L4.797 0H1.564ZM13.556 -9.559 10.498 -22.003H10.358L7.334 -9.559Z",
    '.' => "M7.23 -2.659Q7.23 -5.857 10.428 -5.857Q13.626 -5.857 13.626 -2.659Q13.626 0.539 10.428 0.539Q7.23 0.539 7.23 -2.659Z",
    'j' => "M10.15 4.415Q11.21 4.415 11.966 3.806Q12.722 3.198 13.105 2.251Q13.487 1.303 13.487 0.313V-16.511H7.161V-19.118H16.476V0.243Q16.476 2.294 15.694 3.815Q14.912 5.336 13.496 6.144Q12.079 6.952 10.185 6.952Q8.273 6.952 6.804 6.126Q5.336 5.301 4.51 3.719Q3.685 2.138 3.685 -0.035H6.57Q6.57 2.086 7.517 3.25Q8.464 4.415 10.15 4.415ZM11.784 -25.236Q11.784 -26.487 12.305 -27.009Q12.826 -27.53 14.078 -27.53Q15.329 -27.53 15.851 -27.009Q16.372 -26.487 16.372 -25.236Q16.372 -23.984 15.851 -23.463Q15.329 -22.942 14.078 -22.942Q12.826 -22.942 12.305 -23.463Q11.784 -23.984 11.784 -25.236Z",
    'l' => "M8.116 -5.753V-24.975H2.103V-27.46H11.158V-6.274Q11.158 -4.189 12.07 -3.241Q12.983 -2.294 14.425 -2.294Q15.364 -2.294 16.163 -2.772Q16.963 -3.25 17.936 -4.206L19.674 -2.207Q17.276 0.417 14.304 0.417Q11.088 0.417 9.602 -1.043Q8.116 -2.503 8.116 -5.753Z",
)

fmt(x) = string(round(x; digits = 2))
svgpath(pts) = "M " * join((fmt(p[1]) * " " * fmt(p[2]) for p in pts), " L ")

# Logo drawing commands with the logo square at [x0, x0 + L] × [y0, y0 + L] (SVG coordinates).
function logo_svg(patches, patch_face; x0, y0, nsample = 6)
    io = IOBuffer()
    tosvg(x) = (x0 + L * x[1], y0 + L * (1 - x[2]))
    # Bounding box, as in the Ferrite.jl logo (drawn first so the patches cover its inner half)
    println(io, "<path style=\"fill:none;stroke:black;stroke-width:3.9851\" d=\"$(svgpath(tosvg.([(0, 0), (1, 0), (1, 1), (0, 1)]))) Z\"/>")
    for (ip, patch) in enumerate(patches)
        X(ξ, η) = tosvg(eval_parametric_coordinate(patch, Vec{2}((ξ, η))))
        knots = unique.(patch.knot_vectors)
        # Parameter values along a full edge, sampled between every pair of knots.
        sample(kv) = unique!(reduce(vcat, [range(a, b; length = nsample) for (a, b) in zip(kv[1:end-1], kv[2:end])]))
        tξ, tη = sample(knots[1]), sample(knots[2])
        boundary = [[X(s, -1.0) for s in tξ]; [X(1.0, s) for s in tη]; [X(s, 1.0) for s in reverse(tξ)]; [X(-1.0, s) for s in reverse(tη)]]
        println(io, "<path style=\"fill:$(FACES[patch_face[ip]][1]);stroke:none\" d=\"$(svgpath(boundary)) Z\"/>")
        lines = [[[X(ξ, s) for s in tη] for ξ in knots[1]]; [[X(s, η) for s in tξ] for η in knots[2]]]
        println(io, "<path style=\"fill:none;stroke:black;stroke-width:0.3985;stroke-linejoin:round\" d=\"",
            join(svgpath.(lines), " "), "\"/>")
    end
    return String(take!(io))
end

svgfile(w, h, x0, y0, body) = """
    <?xml version="1.0" encoding="UTF-8"?>
    <svg xmlns="http://www.w3.org/2000/svg" width="$(fmt(w))pt" height="$(fmt(h))pt" viewBox="$(fmt(x0)) $(fmt(y0)) $(fmt(w)) $(fmt(h))" version="1.1">
    $(body)</svg>
    """

function write_logos(dir; kwargs...)
    _, patches, patch_face = ferrite_iga_logo(; kwargs...)
    logo = logo_svg(patches, patch_face; x0 = 0.0, y0 = 0.0)

    # Logo only
    write(joinpath(dir, "logo.svg"), svgfile(L + 2MARGIN, L + 2MARGIN, -MARGIN, -MARGIN, logo))

    # Logo with name below. Placement as in Ferrite.jl's logo-name.svg: the baseline is
    # 44.2 pt below the logo and the text is centered 3.11 pt (0.11 cm) to the right of the
    # logo center, which puts the middle point between the vertical bars of "F" and "l"
    # below the center of the logo.
    name = "FerriteIGA.jl"
    textwidth = GLYPH_ADVANCE * length(name)
    tx, ty = L / 2 + 3.11 - textwidth / 2, L + 44.2
    text = join(("<path transform=\"translate($(fmt(tx + (i - 1) * GLYPH_ADVANCE)),$(fmt(ty)))\" d=\"$(GLYPHS[c])\"/>"
                 for (i, c) in enumerate(name)), "\n")
    x0, y0 = min(tx - 3.578, -MARGIN), -8.96
    w = max(tx + textwidth + 3.718, L + MARGIN) - x0
    h = ty + 10.907 - y0
    write(joinpath(dir, "logo-name.svg"), svgfile(w, h, x0, y0, logo * "<g style=\"fill:black;stroke:none\">\n" * text * "\n</g>\n"))
    return nothing
end

write_logos(joinpath(@__DIR__, "src", "assets"))
