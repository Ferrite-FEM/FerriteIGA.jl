#Solve -Δu + c u = f, with the constraints added by `add_constraints!(ch, grid)`. Returns the L2 error, the dof handler and the solution.
function _solve_with_constraints(grid, p, uex, f, add_constraints!; c = 0.0)
    ip = IGAInterpolation{RefQuadrilateral,p}()
    cv = BezierCellValues(QuadratureRule{RefQuadrilateral}(p+2), ip)
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)
    ch = ConstraintHandler(dh)
    add_constraints!(ch, grid)
    close!(ch)

    K = allocate_matrix(dh, ch)
    F = zeros(ndofs(dh))
    assembler = start_assemble(K, F)
    for cellid in 1:getncells(grid)
        coords = getcoordinates(grid, cellid)
        reinit!(cv, coords)
        n = getnbasefunctions(cv)
        Ke = zeros(n, n)
        fe = zeros(n)
        for qp in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, qp)
            x = spatial_coordinate(cv, qp, coords)
            for i in 1:n
                δu = shape_value(cv, qp, i)
                fe[i] += f(x) * δu * dΩ
                for j in 1:n
                    Ke[i,j] += (shape_gradient(cv, qp, i) ⋅ shape_gradient(cv, qp, j) + c * δu * shape_value(cv, qp, j)) * dΩ
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
    return sqrt(err), dh, ch, a
end

function _nurbs_square_grid(nel, p)
    grid = BezierGrid(generate_nurbs_patch(:rectangle, (nel, nel), p; size = (1.0, 1.0)))
    addfacetset!(grid, "left",   x -> x[1] ≈ 0.0)
    addfacetset!(grid, "right",  x -> x[1] ≈ 1.0)
    addfacetset!(grid, "bottom", x -> x[2] ≈ 0.0)
    addfacetset!(grid, "top",    x -> x[2] ≈ 1.0)
    return grid
end

const _boundary_names = ("left", "right", "bottom", "top")

@testset "ProjectedIGADirichlet" begin
    #A quadratic solution is in the spline space for p >= 2. With projected boundary values, the solution is exact,
    #while evaluating the boundary values at the control points is not.
    uex(x) = x[1]^2 + x[2]^2
    f(x) = -4.0
    grids = [
        (2, BezierGrid(generate_tspline_patch(:graded_rectangle, [2, 4], [1, 2], 2; size = (2.0, 1.0)))),
        (3, BezierGrid(generate_tspline_patch(:graded_rectangle, [2, 4, 8], [1, 2, 4], 3; size = (2.0, 1.0)))),
        (2, _nurbs_square_grid(4, 2)),
    ]
    for (p, grid) in grids
        err_projected, _ = _solve_with_constraints(grid, p, uex, f, (ch, grid) -> begin
            add!(ch, ProjectedIGADirichlet(:u, union((getfacetset(grid, name) for name in _boundary_names)...), uex))
        end)
        err_controlpoints, _ = _solve_with_constraints(grid, p, uex, f, (ch, grid) -> begin
            add!(ch, Dirichlet(:u, union((getfacetset(grid, name) for name in _boundary_names)...), uex))
        end)
        @test err_projected < 1e-10
        @test err_controlpoints > 1e-6
    end

    #Vector field: only the second component is constrained, and the values are kept by update!
    grid = BezierGrid(generate_tspline_patch(:graded_rectangle, [2, 4], [1, 2], 2))
    dh = DofHandler(grid)
    add!(dh, :u, IGAInterpolation{RefQuadrilateral,2}()^2)
    close!(dh)
    ch = ConstraintHandler(dh)
    add!(ch, ProjectedIGADirichlet(:u, getfacetset(grid, "bottom"), x -> 3.0; components = [2]))
    close!(ch)
    update!(ch, 1.0)
    bottom_nodes = getnodeset(grid, "bottom")
    node_dofs = FerriteIGA._node_component_dofs(dh, :u, bottom_nodes)
    @test sort(ch.prescribed_dofs) == sort([node_dofs[A][2] for A in bottom_nodes])
    @test all(≈(3.0), ch.inhomogeneities)

    @test_throws ArgumentError add!(ConstraintHandler(dh), ProjectedIGADirichlet(:u, getfacetset(grid, "bottom"), x -> (1.0, 2.0); components = [3]))
end

@testset "PeriodicIGADirichlet" begin
    #Periodic in x, zero at the bottom and top, on a graded T-spline mesh
    uex(x) = cos(2π*x[1]) * sin(π*x[2])
    f(x) = 5π^2 * uex(x)
    grid = BezierGrid(generate_tspline_patch(:graded_rectangle, [8, 16], [4, 8], 3))
    err, dh, ch, a = _solve_with_constraints(grid, 3, uex, f, (ch, grid) -> begin
        add!(ch, Dirichlet(:u, union(getfacetset(grid, "bottom"), getfacetset(grid, "top")), x -> 0.0))
        add!(ch, PeriodicIGADirichlet(:u, getfacetset(grid, "left"), getfacetset(grid, "right")))
    end)
    @test err < 1e-3
    left, right = collect(getnodeset(grid, "left")), collect(getnodeset(grid, "right"))
    node_dofs = FerriteIGA._node_component_dofs(dh, :u, Set([left; right]))
    for A in left
        B = right[findfirst(B -> grid.nodes[B].x[2] ≈ grid.nodes[A].x[2], right)]
        @test a[only(node_dofs[A])] ≈ a[only(node_dofs[B])]
    end

    #Doubly periodic NURBS square. The corners are tied in both directions without nested constraints.
    uex2(x) = cos(2π*x[1]) * cos(2π*x[2])
    f2(x) = (8π^2 + 1) * uex2(x)
    grid = _nurbs_square_grid(8, 3)
    err, dh, ch, a = _solve_with_constraints(grid, 3, uex2, f2, (ch, grid) -> begin
        add!(ch, PeriodicIGADirichlet(:u, getfacetset(grid, "left"), getfacetset(grid, "right")))
        add!(ch, PeriodicIGADirichlet(:u, getfacetset(grid, "bottom"), getfacetset(grid, "top")))
    end; c = 1.0)
    @test err < 1e-3
    corners = [findfirst(n -> n.x ≈ Vec(x), grid.nodes) for x in ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0))]
    corner_dofs = only.(values(FerriteIGA._node_component_dofs(dh, :u, Set(corners))))
    @test allequal(a[corner_dofs])

    #The sides must be identical
    grid = BezierGrid(generate_tspline_patch(:graded_rectangle, [2, 4], [1, 2], 3))
    dh = DofHandler(grid)
    add!(dh, :u, IGAInterpolation{RefQuadrilateral,3}())
    close!(dh)
    @test_throws ArgumentError add!(ConstraintHandler(dh), PeriodicIGADirichlet(:u, getfacetset(grid, "bottom"), getfacetset(grid, "top")))
    @test_throws ArgumentError add!(ConstraintHandler(dh), ProjectedDirichlet(:u, getfacetset(grid, "bottom"), x -> 0.0))
end
