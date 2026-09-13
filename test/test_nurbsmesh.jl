

@testset "test nurbsmesh" begin
    
    mesh = generate_nurbs_patch(:plate_with_hole, (4,4), (2,2))

    x = eval_parametric_coordinate(mesh, Vec((-1.0,-1.0)))
    @test x == Vec((-1.0, 0.0))

    x = eval_parametric_coordinate(mesh, Vec((1.0,1.0)))
    @test x == Vec((0.0, 4.0))

    x = eval_parametric_coordinate(mesh, Vec((-1.0,1.0)))
    @test x == Vec((-4.0, 0.0))

end

@testset "knotinsertion!" begin

    # creates a single element patch 
    coarse = generate_nurbs_patch(:hypercube, (1, 1), (2, 2); cornerpos=(0.0, 0.0), size=(2.0, 3.0))
    kv = (copy(coarse.knot_vectors[1]), copy(coarse.knot_vectors[2]))
    cp = copy(coarse.control_points)
    w = copy(coarse.weights)
    orders = (2, 2)
    x = eval_parametric_coordinate(coarse, Vec(0.5, 0.75))

    # splits one span in ξ, geometry should be same and cell count should double
    FerriteIGA.knotinsertion!(kv, orders, cp, w, 0.0, dir=1)
    refined = NURBSMesh(kv, orders, cp, w)
    @test eval_parametric_coordinate(refined, Vec(0.5, 0.75)) ≈ x
    @test getncells(refined) == 2

    # splits one span in η, result should match a 2x2 patch with 4 total elements
    FerriteIGA.knotinsertion!(kv, orders, cp, w, 0.0, dir=2)
    manual = NURBSMesh(kv, orders, cp, w)
    reference = generate_nurbs_patch(:hypercube, (2, 2), (2, 2); cornerpos=(0.0, 0.0), size=(2.0, 3.0))
    @test getncells(manual) == getncells(reference)
    @test manual.knot_vectors == reference.knot_vectors
end

@testset "orderelevation!" begin
    coarse = generate_nurbs_patch(:hypercube, (1, 1), (2, 2); cornerpos=(0.0, 0.0), size=(2.0, 3.0))
    x = eval_parametric_coordinate(coarse, Vec(0.5, 0.75))

    #single span test should elevate from p=2 to p=3
    kv = (copy(coarse.knot_vectors[1]), copy(coarse.knot_vectors[2]))
    cp = copy(coarse.control_points)
    w = copy(coarse.weights)
    orders = (2, 2)
    orders = FerriteIGA.orderelevation!(kv, orders, cp, w, dir=1)
    elevated = NURBSMesh(kv, orders, cp, w)
    @test eval_parametric_coordinate(elevated, Vec(0.5, 0.75)) ≈ x
    @test kv[1] == [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0]

    #double element patch, so knot vector should change from [-1,-1,-1,0,1,1,1] to [-1,-1,-1,-1,0,0,1,1,1,1]
    kv = (copy(coarse.knot_vectors[1]), copy(coarse.knot_vectors[2]))
    cp = copy(coarse.control_points)
    w = copy(coarse.weights)
    orders = (2, 2)
    FerriteIGA.knotinsertion!(kv, orders, cp, w, 0.0, dir=1)
    orders = FerriteIGA.orderelevation!(kv, orders, cp, w, dir=1)
    elevated = NURBSMesh(kv, orders, cp, w)
    @test eval_parametric_coordinate(elevated, Vec(0.5, 0.75)) ≈ x
    @test getncells(elevated) == 2
    @test kv[1] == [-1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    #hypercube test that area and boundary lengths are unchanged after p-refinement
    mesh = generate_nurbs_patch(:hypercube, (4, 4), (2, 2); cornerpos=(-1.0, -1.0), size=(2.0, 3.0))
    grid = BezierGrid(mesh)
    addcellset!(grid, "all", (x) -> true)
    addfacetset!(grid, "left", (x) -> x[1] ≈ -1.0)
    addfacetset!(grid, "right", (x) -> x[1] ≈ 1.0)
    addfacetset!(grid, "top", (x) -> x[2] ≈ 2.0)
    addfacetset!(grid, "bottom", (x) -> x[2] ≈ -1.0)
    bern_ip = IGAInterpolation{RefQuadrilateral, 2}()
    qr = QuadratureRule{RefQuadrilateral}(5)
    cv = BezierCellValues(qr, bern_ip)
    fqr = FacetQuadratureRule{RefQuadrilateral}(5)
    fv = BezierFacetValues(fqr, bern_ip)
    V0 = _calculate_volume(cv, grid, getcellset(grid, "all"))
    A0_left = _calculate_area(fv, grid, getfacetset(grid, "left"))
    A0_right = _calculate_area(fv, grid, getfacetset(grid, "right"))
    A0_top = _calculate_area(fv, grid, getfacetset(grid, "top"))
    A0_bottom = _calculate_area(fv, grid, getfacetset(grid, "bottom"))

    kv = (copy(mesh.knot_vectors[1]), copy(mesh.knot_vectors[2]))
    cp = copy(mesh.control_points)
    w = copy(mesh.weights)
    orders = (2, 2)
    orders = FerriteIGA.orderelevation!(kv, orders, cp, w, dir=1)
    orders = FerriteIGA.orderelevation!(kv, orders, cp, w, dir=2)
    elevated = NURBSMesh(kv, orders, cp, w)
    @test orders == (3, 3)

    grid = BezierGrid(elevated)
    addcellset!(grid, "all", (x) -> true)
    addfacetset!(grid, "left", (x) -> x[1] ≈ -1.0)
    addfacetset!(grid, "right", (x) -> x[1] ≈ 1.0)
    addfacetset!(grid, "top", (x) -> x[2] ≈ 2.0)
    addfacetset!(grid, "bottom", (x) -> x[2] ≈ -1.0)
    bern_ip = IGAInterpolation{RefQuadrilateral, 3}()
    cv = BezierCellValues(qr, bern_ip)
    fv = BezierFacetValues(fqr, bern_ip)
    @test _calculate_volume(cv, grid, getcellset(grid, "all")) ≈ V0
    @test _calculate_area(fv, grid, getfacetset(grid, "left")) ≈ A0_left
    @test _calculate_area(fv, grid, getfacetset(grid, "right")) ≈ A0_right
    @test _calculate_area(fv, grid, getfacetset(grid, "top")) ≈ A0_top
    @test _calculate_area(fv, grid, getfacetset(grid, "bottom")) ≈ A0_bottom

    #plate with hole test that volume and boundary lengths unchanged after p-refinement
    L = 4.0
    r = 1.0
    mesh = generate_nurbs_patch(:plate_with_hole, (4, 4), (2, 2))
    grid = BezierGrid(mesh)
    addcellset!(grid, "all", (x) -> true)
    addfacetset!(grid, "left", (x) -> x[1] ≈ -L)
    addfacetset!(grid, "top", (x) -> x[2] ≈ L)
    addfacetset!(grid, "right", (x) -> x[1] ≈ 0.0)
    addfacetset!(grid, "bottom", (x) -> x[2] ≈ 0.0)
    addfacetset!(grid, "circle", (x) -> r * 0.9 < norm(x) < r * 1.1)
    bern_ip = IGAInterpolation{RefQuadrilateral, 2}()
    cv = BezierCellValues(qr, bern_ip)
    fv = BezierFacetValues(fqr, bern_ip)
    V0 = _calculate_volume(cv, grid, getcellset(grid, "all"))
    A0_left = _calculate_area(fv, grid, getfacetset(grid, "left"))
    A0_right = _calculate_area(fv, grid, getfacetset(grid, "right"))
    A0_top = _calculate_area(fv, grid, getfacetset(grid, "top"))
    A0_bottom = _calculate_area(fv, grid, getfacetset(grid, "bottom"))
    A0_circle = _calculate_area(fv, grid, getfacetset(grid, "circle"))

    kv = (copy(mesh.knot_vectors[1]), copy(mesh.knot_vectors[2]))
    cp = copy(mesh.control_points)
    w = copy(mesh.weights)
    orders = (2, 2)
    orders = FerriteIGA.orderelevation!(kv, orders, cp, w, dir=1)
    orders = FerriteIGA.orderelevation!(kv, orders, cp, w, dir=2)
    elevated = NURBSMesh(kv, orders, cp, w)
    @test orders == (3, 3)

    grid = BezierGrid(elevated)
    addcellset!(grid, "all", (x) -> true)
    addfacetset!(grid, "left", (x) -> x[1] ≈ -L)
    addfacetset!(grid, "top", (x) -> x[2] ≈ L)
    addfacetset!(grid, "right", (x) -> x[1] ≈ 0.0)
    addfacetset!(grid, "bottom", (x) -> x[2] ≈ 0.0)
    addfacetset!(grid, "circle", (x) -> r * 0.9 < norm(x) < r * 1.1)
    bern_ip = IGAInterpolation{RefQuadrilateral, 3}()
    cv = BezierCellValues(qr, bern_ip)
    fv = BezierFacetValues(fqr, bern_ip)
    @test _calculate_volume(cv, grid, getcellset(grid, "all")) ≈ V0
    @test _calculate_area(fv, grid, getfacetset(grid, "left")) ≈ A0_left
    @test _calculate_area(fv, grid, getfacetset(grid, "right")) ≈ A0_right
    @test _calculate_area(fv, grid, getfacetset(grid, "top")) ≈ A0_top
    @test _calculate_area(fv, grid, getfacetset(grid, "bottom")) ≈ A0_bottom
    @test _calculate_area(fv, grid, getfacetset(grid, "circle")) ≈ A0_circle
end

@testset "smoothnesselevation!" begin
    coarse = generate_nurbs_patch(:hypercube, (1, 1), (2, 2); cornerpos=(0.0, 0.0), size=(2.0, 3.0))
    x = eval_parametric_coordinate(coarse, Vec(0.5, 0.75))

    #k-refinement, elevate like p-refinement, insert like h-refinement
    kv = (copy(coarse.knot_vectors[1]), copy(coarse.knot_vectors[2]))
    cp = copy(coarse.control_points)
    w = copy(coarse.weights)
    orders = (2, 2)
    orders = FerriteIGA.smoothnesselevation!(kv, orders, cp, w, [0.0]; dir=1)
    refined = NURBSMesh(kv, orders, cp, w)
    @test eval_parametric_coordinate(refined, Vec(0.5, 0.75)) ≈ x
    @test orders == (3, 2)
    @test getncells(refined) == 2
    @test kv[1] == [-1.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    #should match calling orderelevation! then knotinsertion! by hand
    kv_ref = (copy(coarse.knot_vectors[1]), copy(coarse.knot_vectors[2]))
    cp_ref = copy(coarse.control_points)
    w_ref = copy(coarse.weights)
    orders_ref = (2, 2)
    orders_ref = FerriteIGA.orderelevation!(kv_ref, orders_ref, cp_ref, w_ref; dir=1)
    FerriteIGA.knotinsertion!(kv_ref, orders_ref, cp_ref, w_ref, 0.0; dir=1)
    @test orders == orders_ref
    @test kv == kv_ref
    @test cp ≈ cp_ref
    @test w ≈ w_ref

    #smoothness check, the new knot from k-refinement should be one more derivative continuous, while h-then-p raises multiplicity and continuity remains the same
    @test count(==(0.0), kv[1]) == 1
    kv_hp = (copy(coarse.knot_vectors[1]), copy(coarse.knot_vectors[2]))
    cp_hp = copy(coarse.control_points)
    w_hp = copy(coarse.weights)
    orders_hp = (2, 2)
    FerriteIGA.knotinsertion!(kv_hp, orders_hp, cp_hp, w_hp, 0.0; dir=1)
    orders_hp = FerriteIGA.orderelevation!(kv_hp, orders_hp, cp_hp, w_hp; dir=1)
    @test orders_hp == (3, 2)
    @test count(==(0.0), kv_hp[1]) == 2
    @test kv_hp[1] == [-1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
end

@testset "Test grid to BezierGrid convertion" begin
    
    grid = Ferrite.generate_grid(QuadraticQuadrilateral, (4,4))
    bgrid = BezierGrid(grid)

    cellid = 1
    nnodes = length(bgrid.cells[cellid].nodes)
    w = zeros(Float64, nnodes)
    getweights!(w, bgrid, cellid)
    @test all(w .== 1.0)

    C = get_extraction_operator(bgrid, cellid)
    @test FerriteIGA.beo2matrix(C) == diagm(ones(Float64, nnodes))

    @test getncells(bgrid) == getncells(grid)
    @test getnnodes(bgrid) == getnnodes(grid)
end

@testset "NURBSMesh refinement functions" begin
    coarse = generate_nurbs_patch(:hypercube, (1, 1), (2, 2); cornerpos=(0.0, 0.0), size=(2.0, 3.0))

    refined = FerriteIGA.knotinsertion(coarse, 0.0; dir=1)
    kv = (copy(coarse.knot_vectors[1]), copy(coarse.knot_vectors[2]))
    cp = copy(coarse.control_points)
    w = copy(coarse.weights)
    FerriteIGA.knotinsertion!(kv, coarse.orders, cp, w, 0.0; dir=1)
    raw = NURBSMesh(kv, coarse.orders, cp, w)
    @test refined.knot_vectors == raw.knot_vectors
    @test refined.orders == raw.orders
    @test refined.control_points ≈ raw.control_points
    @test refined.weights ≈ raw.weights
    @test getncells(refined) == getncells(raw)
    @test getncells(coarse) == 1

    elevated = FerriteIGA.orderelevation(coarse; dir=1)
    kv = (copy(coarse.knot_vectors[1]), copy(coarse.knot_vectors[2]))
    cp = copy(coarse.control_points)
    w = copy(coarse.weights)
    orders = FerriteIGA.orderelevation!(kv, coarse.orders, cp, w; dir=1)
    raw = NURBSMesh(kv, orders, cp, w)
    @test elevated.knot_vectors == raw.knot_vectors
    @test elevated.orders == raw.orders
    @test elevated.control_points ≈ raw.control_points
    @test elevated.weights ≈ raw.weights
    @test getncells(elevated) == getncells(raw)

    smooth = FerriteIGA.smoothnesselevation(coarse, [0.0]; dir=1)
    kv = (copy(coarse.knot_vectors[1]), copy(coarse.knot_vectors[2]))
    cp = copy(coarse.control_points)
    w = copy(coarse.weights)
    orders = FerriteIGA.smoothnesselevation!(kv, coarse.orders, cp, w, [0.0]; dir=1)
    raw = NURBSMesh(kv, orders, cp, w)
    @test smooth.knot_vectors == raw.knot_vectors
    @test smooth.orders == raw.orders
    @test smooth.control_points ≈ raw.control_points
    @test smooth.weights ≈ raw.weights
    @test getncells(smooth) == getncells(raw)
end