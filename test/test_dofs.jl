@testset "Multiple subdofhandlers" begin
    grid = generate_grid(BezierCell{RefQuadrilateral,2}, (2,2))

    ip = IGAInterpolation{RefQuadrilateral,2}()
    dh = DofHandler(grid)
    sdh1 = SubDofHandler(dh, 1:2)
    add!(sdh1, :u, ip^2)
    sdh2 = SubDofHandler(dh, 3:4)
    add!(sdh2, :p, ip)
    close!(dh)

    urange = dof_range(sdh1, :u)
    prange = dof_range(sdh2, :p)
    a = zeros(Float64, ndofs(dh))

    for cellid in 1:2
        dofs = celldofs(dh, cellid)
        a[dofs[urange]] .= 2.0
    end

    for cellid in 3:4
        dofs = celldofs(dh, cellid)
        a[dofs[prange]] .= 1.0
    end

    mktempdir() do dir
        @test_nowarn VTKIGAFile(joinpath(dir, "test_export"), grid) do vtk
            write_solution(vtk, dh, a)
        end
    end
end

@testset "Multi-fields for iga" begin
    grid = generate_grid(BezierCell{RefQuadrilateral,2}, (2,2))

    ip = IGAInterpolation{RefQuadrilateral,2}()
    dh = DofHandler(grid)
    add!(dh, :u, ip^2)
    add!(dh, :p, ip)
    close!(dh)

    urange = dof_range(dh, :u)
    prange = dof_range(dh, :p)
    a = zeros(Float64, ndofs(dh))

    dofs = zeros(Int, ndofs_per_cell(dh))
    for cellid in 1:getncells(grid)
        celldofs!(dofs, dh, cellid)
        a[dofs[urange]] .= 1.0
        a[dofs[prange]] .= 2.0
    end

    mktempdir() do dir
        @test_nowarn VTKIGAFile(joinpath(dir, "test_export"), grid) do vtk
            write_solution(vtk, dh, a)
        end
    end
end