using FerriteIGA: _bspline_bezier_coefficients, _compute_bezier_extraction_operators, compute_bezier_extraction_operators

#Faces of a T-mesh in index space, from predicates for the unit edges:
#hedge(i,j) for the edge (i,j)-(i+1,j), and vedge(i,j) for the edge (i,j)-(i,j+1).
function _faces_from_edges(n1, n2, hedge, vedge)
    owner = zeros(Int, n1-1, n2-1)
    faces = NTuple{4,Int}[]
    for j in 1:n2-1, i in 1:n1-1
        owner[i,j] != 0 && continue
        i1 = i+1; while i1 < n1 && !vedge(i1, j); i1 += 1; end
        j1 = j+1; while j1 < n2 && !hedge(i, j1); j1 += 1; end
        push!(faces, (i, i1, j, j1))
        owner[i:i1-1, j:j1-1] .= length(faces)
    end
    return faces
end

#Cubic T-spline on [0,n]², with full lines at the integers and additional vertical (x, y0, y1) and
#horizontal (y, x0, x1) segments.
function _tmesh_with_segments(vsegs, hsegs; p = 3, n = 8)
    vals = sort(unique(Float64[0:n; [s[1] for s in vsegs]; [s[1] for s in hsegs]]))
    r = p÷2
    knots = ([fill(0.0, r); vals; fill(Float64(n), r)], [fill(0.0, r); vals; fill(Float64(n), r)])
    N = length(knots[1])
    vedge(i,j) = (x = knots[1][i]; isinteger(x) || any(s -> s[1] == x && s[2] <= knots[2][j] && knots[2][j+1] <= s[3], vsegs))
    hedge(i,j) = (y = knots[2][j]; isinteger(y) || any(s -> s[1] == y && s[2] <= knots[1][i] && knots[1][i+1] <= s[3], hsegs))
    return TSplineMesh(knots, p, _faces_from_edges(N, N, hedge, vedge))
end

#Solve -Δu + c u = f with Dirichlet conditions u = uex on the facetsets `dbc_sets`. Returns the L2 error.
function _solve_tspline_poisson(grid, p, uex, f; dbc_sets = ("left", "right", "bottom", "top"), c = 0.0)
    ip = IGAInterpolation{RefQuadrilateral,p}()
    cv = BezierCellValues(QuadratureRule{RefQuadrilateral}(p+2), ip)
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)
    ch = ConstraintHandler(dh)
    for set in dbc_sets
        add!(ch, Dirichlet(:u, getfacetset(grid, set), x -> uex(x)))
    end
    close!(ch)
    update!(ch, 0.0)

    K = allocate_matrix(dh)
    F = zeros(ndofs(dh))
    assembler = start_assemble(K, F)
    for cellid in 1:getncells(grid)
        coords = getcoordinates(grid, cellid)
        reinit!(cv, coords)
        n = getnbasefunctions(cv)
        @assert n == ndofs_per_cell(dh, cellid)
        Ke = zeros(n, n)
        fe = zeros(n)
        for qp in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, qp)
            x = spatial_coordinate(cv, qp, coords)
            for i in 1:n
                δu = shape_value(cv, qp, i)
                ∇δu = shape_gradient(cv, qp, i)
                fe[i] += f(x) * δu * dΩ
                for j in 1:n
                    Ke[i,j] += (∇δu ⋅ shape_gradient(cv, qp, j) + c * δu * shape_value(cv, qp, j)) * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(dh, cellid), Ke, fe)
    end
    apply!(K, F, ch)
    a = K \ F
    apply!(a, ch)

    err = 0.0
    for cellid in 1:getncells(grid)
        coords = getcoordinates(grid, cellid)
        reinit!(cv, coords)
        ae = a[celldofs(dh, cellid)]
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            err += (function_value(cv, qp, ae) - uex(x))^2 * getdetJdV(cv, qp)
        end
    end
    return sqrt(err), dh, a
end

@testset "T-splines: univariate extraction" begin
    for (p, U) in ((2, [0,0,0,1,2,2,3,4,4,4.0]), (3, [0,0,0,0,1,2,2,3,4,4,4,4.0]), (3, [0,0,0,0,0.5,2,2,2,3,4,4,4,4.0]))
        C, _ = _compute_bezier_extraction_operators(p, U)
        spans = [(U[k], U[k+1]) for k in 1:length(U)-1 if U[k] < U[k+1]]
        for (e, (a, b)) in enumerate(spans)
            k = findlast(==(a), U)
            for (row, i) in enumerate((k-p):k)
                @test _bspline_bezier_coefficients(U[i:i+p+1], p, a, b) ≈ C[e][row, :] atol=1e-14
            end
        end
    end
    #A T-spline function where the element is a sub-interval of a knot span. 
    #The uniform cubic B-spline is (-3u³ + 3u² + 3u + 1)/6 for x = 1 + u ∈ [1,2].
    c = _bspline_bezier_coefficients([0.0, 1.0, 2.0, 3.0, 4.0], 3, 1.0, 1.5)
    for t in (0.0, 0.3, 0.7, 1.0)
        u = 0.5t
        @test sum(c[j+1] * binomial(3, j) * t^j * (1-t)^(3-j) for j in 0:3) ≈ (-3u^3 + 3u^2 + 3u + 1)/6
    end
end

@testset "T-splines: NURBS as T-splines" begin
    for (p, nurbs) in (
            (2, generate_nurbs_patch(:rectangle, (3, 4), 2; size = (2.0, 3.0))),
            (3, generate_nurbs_patch(:rectangle, (4, 3), 3; size = (2.0, 3.0))),
            (2, generate_nurbs_patch(:plate_with_hole, (4, 3), (2, 2))),
        )
        nurbs = knotinsertion(nurbs, nurbs.knot_vectors[1][p+2]; dir=1) #Repeated interior knot
        tmesh = TSplineMesh(nurbs)
        @test getnnodes(tmesh) == getnnodes(nurbs)
        @test getncells(tmesh) == getncells(nurbs)
        @test is_linearly_independent(tmesh)
        isodd(p) && @test is_analysis_suitable(tmesh)

        nurbsgrid = BezierGrid(nurbs)
        tgrid = BezierGrid(tmesh)
        @test tmesh.control_points == nurbs.control_points && tmesh.weights == nurbs.weights
        @test eltype(tgrid.cells) == TSplineCell{RefQuadrilateral, p, (p+1)^2}
        for cellid in 1:getncells(nurbsgrid)
            xb1, wb1 = get_bezier_coordinates(nurbsgrid, cellid)
            xb2, wb2 = get_bezier_coordinates(tgrid, cellid)
            @test xb1 ≈ xb2
            @test wb1 ≈ wb2
            #Same extraction operator rows for the same basis functions. The anchors are sorted row by row,
            #so the control points have the same numbering as in the NURBS patch.
            C1, C2 = nurbsgrid.beo[cellid].C, tgrid.beo[cellid].C
            nodes1, nodes2 = nurbsgrid.cells[cellid].nodes, tgrid.cells[cellid].nodes
            @test sort(collect(nodes1)) == collect(nodes2)
            for (a, A) in pairs(nodes2)
                b = findfirst(==(A), nodes1)
                @test Vector(C1[b]) ≈ Vector(C2[a])
            end
        end

        #Same discrete solution (reaction-diffusion with natural boundary conditions)
        uex(x) = x[1]^2 + x[2]
        f(x) = 1 + x[1]
        e1, _ = _solve_tspline_poisson(nurbsgrid, p, uex, f; dbc_sets = (), c = 1.0)
        e2, _ = _solve_tspline_poisson(tgrid, p, uex, f; dbc_sets = (), c = 1.0)
        @test e1 ≈ e2 rtol=1e-8
    end
end

@testset "T-splines: graded mesh" begin
    mesh = generate_tspline_patch(:graded_rectangle, [2, 4, 8], [1, 2, 4], 3; size = (2.0, 1.0))
    @test is_analysis_suitable(mesh)
    @test is_linearly_independent(mesh)
    hext, vext = tjunction_extensions(mesh)
    @test isempty(hext) && length(vext) > 0
    #Analysis-suitable T-splines are locally linearly independent: (p+1)^2 functions per element
    @test all(ien -> length(ien) == 16, mesh.IEN)

    #Polynomial partition of unity and non-negative coefficients
    beos = compute_bezier_extraction_operators(mesh)
    for beo in beos
        @test sum(Vector(r) for r in beo.C) ≈ ones(16)
        @test all(r -> all(>=(0), r.nzval), beo.C)
    end
    #Greville points give the identity geometry
    grid = BezierGrid(mesh)
    for cellid in 1:getncells(grid)
        (ξa, ξb, ηa, ηb) = mesh.elements[cellid]
        xb, _ = get_bezier_coordinates(grid, cellid)
        @test xb[1] ≈ Vec((ξa, ηa)) && xb[3] ≈ Vec((ξb, ηb))
    end
    @test length(getnodeset(grid, "bottom")) == 2 + 3 && length(getnodeset(grid, "top")) == 8 + 3
    for (set, nel) in (("bottom", 2), ("top", 8))
        facets = getfacetset(grid, set)
        @test length(facets) >= nel
        @test sum(f -> mesh.elements[f[1]][2] - mesh.elements[f[1]][1], facets) ≈ 2.0
    end

    #Facet integration on the boundary of the graded mesh
    fv = BezierFacetValues(FacetQuadratureRule{RefQuadrilateral}(4), IGAInterpolation{RefQuadrilateral,3}())
    for (set, len) in (("bottom", 2.0), ("right", 1.0), ("top", 2.0), ("left", 1.0))
        L = 0.0
        for (cellid, facetid) in getfacetset(grid, set)
            reinit!(fv, getcoordinates(grid, cellid), facetid)
            for qp in 1:getnquadpoints(fv)
                L += getdetJdV(fv, qp)
                @test sum(i -> shape_value(fv, qp, i), 1:getnbasefunctions(fv)) ≈ 1.0
            end
        end
        @test L ≈ len
    end

    #Even degree
    mesh2 = generate_tspline_patch(:graded_rectangle, [3, 6], [2, 4], 2)
    @test is_analysis_suitable(mesh2)
    @test is_linearly_independent(mesh2)
    @test all(ien -> length(ien) == 9, mesh2.IEN)
end

@testset "T-splines: analysis-suitability for odd and even degrees" begin
    configurations = [
        ([(4.5, 4.0, 5.0)], [(4.5, 4.0, 5.0)]),                                 #A face split into four
        ([(3.5, 3.0, 5.0), (4.5, 3.0, 5.0)], [(3.5, 3.0, 5.0), (4.5, 3.0, 5.0)]), #A refined 2×2 block
        ([(2.5, 2.0, 5.0)], [(5.5, 3.0, 6.0)]),                                 #Perpendicular T-junctions close to each other
        ([(2.5, 1.0, 4.0)], [(6.5, 4.0, 7.0)]),                                 #Perpendicular T-junctions far apart
        ([(2.5, 2.0, 5.0)], [(4.5, 3.0, 6.0)]),
        ([(2.5, 2.0, 4.0)], [(4.5, 5.0, 7.0)]),
        ([(4.5, 2.0, 6.0), (5.5, 3.0, 5.0)], NTuple{3,Float64}[]),              #Only vertical T-junctions
    ]
    expected = Dict(3 => [false, false, false, true, false, true, true], 2 => [false, false, false, true, false, true, true])
    for p in (2, 3)
        for (k, (vsegs, hsegs)) in enumerate(configurations)
            mesh = _tmesh_with_segments(vsegs, hsegs; p)
            AS = is_analysis_suitable(mesh)
            @test AS == expected[p][k]
            #Analysis-suitable T-splines are linearly independent and locally linearly independent, 
            #i.e. there are (p+1)^2 functions on each element. (For these meshes, the reverse also holds.)
            @test AS == all(ien -> length(ien) == (p+1)^2, mesh.IEN)
            AS && @test is_linearly_independent(mesh)
        end
    end
    #Analysis-suitable meshes for degree p are analysis-suitable for lower degrees (same T-mesh for p = 2 and 3)
    @test all(k -> !expected[3][k] || expected[2][k], eachindex(configurations))
end

@testset "T-splines: non analysis-suitable mesh" begin
    #A vertical and a horizontal T-junction close to each other
    mesh = _tmesh_with_segments([(2.5, 2.0, 5.0)], [(5.5, 3.0, 6.0)])
    @test !is_analysis_suitable(mesh)
    @test is_linearly_independent(mesh)
    @test sort(unique(length.(mesh.IEN))) == [16, 17]

    #The grid can be created, with one SubDofHandler per cell type
    grid = BezierGrid(mesh)
    @test !isconcretetype(eltype(grid.cells))
    dh = DofHandler(grid)
    add!(dh, :u, IGAInterpolation{RefQuadrilateral,3}()^2)
    close!(dh)
    @test length(dh.subdofhandlers) == 2
    @test sum(sdh -> count(_ -> true, CellIterator(sdh)), dh.subdofhandlers) == getncells(grid)
    for sdh in dh.subdofhandlers, cc in CellIterator(sdh)
        @test length(celldofs(cc)) == 2 * Ferrite.nnodes(grid.cells[cellid(cc)])
    end

    #BezierCellValues only support square extraction operators, i.e. cells with (p+1)^2 functions
    cv = BezierCellValues(QuadratureRule{RefQuadrilateral}(4), IGAInterpolation{RefQuadrilateral,3}())
    coords16 = getcoordinates(grid, findfirst(c -> Ferrite.nnodes(c) == 16, grid.cells))
    reinit!(cv, coords16)
    @test all(qp -> sum(i -> shape_value(cv, qp, i), 1:16) ≈ 1.0, 1:getnquadpoints(cv))
    coords17 = getcoordinates(grid, findfirst(c -> Ferrite.nnodes(c) == 17, grid.cells))
    @test_throws ArgumentError reinit!(cv, coords17)
end

@testset "T-splines: periodic graded mesh" begin
    for p in (2, 3)
        open_mesh = generate_tspline_patch(:graded_rectangle, [4, 8], [2, 4], p; size = (2.0, 1.0))
        mesh = generate_tspline_patch(:graded_rectangle, [4, 8], [2, 4], p; size = (2.0, 1.0), periodic = true)
        @test mesh.domain == ((0.0, 2.0), (0.0, 1.0))
        @test is_analysis_suitable(mesh)
        @test is_linearly_independent(mesh)
        @test all(ien -> length(ien) == (p+1)^2, mesh.IEN)
        #The knots extend outside of the domain, but the elements do not
        @test mesh.knots[1][1] < 0.0 && mesh.knots[1][end] > 2.0
        @test all(e -> 0.0 <= e[1] && e[2] <= 2.0, mesh.elements)
        @test sum(e -> (e[2]-e[1])*(e[4]-e[3]), mesh.elements) ≈ 2.0
        @test getncells(mesh) == getncells(open_mesh) && getnnodes(mesh) == getnnodes(open_mesh)
        for beo in compute_bezier_extraction_operators(mesh)
            @test sum(Vector(r) for r in beo.C) ≈ ones((p+1)^2)
        end
        #p columns of basis functions are non-zero on the left and right sides
        grid = BezierGrid(mesh)
        @test length(getnodeset(grid, "left")) == length(getnodeset(grid, "right")) == p * length(getnodeset(BezierGrid(open_mesh), "left"))
        @test length(getnodeset(grid, "bottom")) == 4 + p
    end
end

@testset "T-splines: invalid T-meshes" begin
    knots = ([0.0, 0.0, 1.0, 2.0, 2.0], [0.0, 0.0, 1.0, 2.0, 2.0])
    faces = [(i, i+1, j, j+1) for j in 1:4 for i in 1:4]
    @test TSplineMesh(knots, 3, faces) isa TSplineMesh
    @test_throws ArgumentError TSplineMesh(knots, 3, faces[1:end-1])                 #Gap
    @test_throws ArgumentError TSplineMesh(knots, 3, [faces; (1, 3, 1, 2)])          #Overlap
    @test_throws ArgumentError TSplineMesh(knots, 3, [faces[3:end]; (1, 3, 1, 2)])   #Zero interval boundary line not full
    @test_throws ArgumentError TSplineMesh(knots, 3, faces; control_points = zeros(Vec{2,Float64}, 3))
    @test_throws ArgumentError TSplineMesh(knots, 3, faces; domain = ((0.0, 1.5), (0.0, 2.0)))       #Domain cuts elements
end

@testset "T-splines: Poisson" begin
    #Patch test: linear solutions are reproduced exactly
    for p in (1, 2, 3)
        grid = BezierGrid(generate_tspline_patch(:graded_rectangle, [2, 4, 8], [1, 2, 4], p; size = (2.0, 1.0), cornerpos = (-1.0, 0.5)))
        err, _ = _solve_tspline_poisson(grid, p, x -> 1 + 2x[1] + 3x[2], x -> 0.0)
        @test err < 1e-10
    end

    #Convergence for a smooth solution on graded meshes (fine at the top, coarse at the bottom)
    uex(x) = sin(π*x[1]) * sin(π*x[2])
    f(x) = 2π^2 * uex(x)
    errors = map((2, 4, 8)) do k
        grid = BezierGrid(generate_tspline_patch(:graded_rectangle, k .* [2, 4, 8], k .* [1, 2, 4], 3))
        first(_solve_tspline_poisson(grid, 3, uex, f))
    end
    rates = log2.(errors[1:end-1] ./ errors[2:end])
    @test rates[end] > 3.7
end

@testset "T-splines: caches, iterators and VTK export" begin
    #Analysis-suitable mesh with both horizontal and vertical T-junctions
    mesh = _tmesh_with_segments([(2.5, 1.0, 4.0)], [(6.5, 4.0, 7.0)])
    @test is_analysis_suitable(mesh)
    grid = BezierGrid(mesh)
    ip = IGAInterpolation{RefQuadrilateral,3}()
    dh = DofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    cc = IGACellCache(dh)
    for cellid in 1:getncells(grid)
        reinit!(cc, cellid)
        N = Ferrite.nnodes(grid.cells[cellid])
        @test length(getcoordinates(cc).x) == N
        @test length(getcoordinates(cc).xb) == 16
        @test celldofs(cc) == celldofs(dh, cellid)
        @test length(celldofs(cc)) == 2N
    end
    #A single cell type, so Ferrite's CellIterator works as usual
    @test isconcretetype(eltype(grid.cells))
    @test count(_ -> true, CellIterator(dh)) == getncells(grid)

    #Vector valued field: interpolate the identity map, x = Σ Rᴬ xᴬ
    a = zeros(ndofs(dh))
    for cellid in 1:getncells(grid)
        dofs = celldofs(dh, cellid)
        for (i, node) in pairs(grid.cells[cellid].nodes)
            a[dofs[2i-1:2i]] .= grid.nodes[node].x
        end
    end
    cv = BezierCellValues(QuadratureRule{RefQuadrilateral}(4), ip^2)
    for cellid in 1:getncells(grid)
        coords = getcoordinates(grid, cellid)
        reinit!(cv, coords)
        ae = a[celldofs(dh, cellid)]
        for qp in 1:getnquadpoints(cv)
            @test function_value(cv, qp, ae) ≈ spatial_coordinate(cv, qp, coords)
            @test function_gradient(cv, qp, ae) ≈ one(Tensor{2,2})
        end
    end

    mktempdir() do dir
        @test_nowarn VTKIGAFile(joinpath(dir, "tspline"), grid) do vtk
            write_solution(vtk, dh, a)
            write_node_data(vtk, [x[1] for x in mesh.control_points], "x")
        end
    end
end
