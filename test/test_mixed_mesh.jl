function _offset_cell(cell, offset)
    nodes = ntuple(i -> cell.nodes[i] + offset, length(cell.nodes))
    return typeof(cell)(nodes)
end

function gridmerge(grids...)
    CT = Union{map(g -> eltype(g.cells), grids)...}

    cells = CT[]
    nodes = Node{2,Float64}[]
    beos = Union{BezierExtractionOperator{Float64}, Nothing}[]
    weights = Float64[]

    node_offset = 0

    for grid in grids
        # Append cells with updated node indices
        append!(cells, [_offset_cell(cell, node_offset) for cell in grid.cells])

        # Append nodes
        append!(nodes, grid.nodes)

        # Append Bezier extraction operators
        if hasproperty(grid, :beo)
            append!(beos, grid.beo)
        else
            append!(beos, fill(nothing, getncells(grid)))
        end

        # Append weights
        if hasproperty(grid, :weights)
            append!(weights, grid.weights)
        else
            append!(weights, fill(1.0, getnnodes(grid)))
        end

        node_offset += getnnodes(grid)
    end

    return BezierGrid(cells, nodes, weights, beos)
end

@testset "Mixed meshes" begin
    grid1 = generate_grid(BezierCell{RefQuadrilateral, 2}, (2, 2), Vec((0.0, 0.0)), Vec((1.0, 1.0)))
    grid2 = generate_grid(Quadrilateral, (2, 2), Vec((2.0, 2.0)), Vec((3.0, 3.0)))
    grid3 = generate_grid(Triangle, (1, 1), Vec((4.0, 4.0)), Vec((5.0, 5.0)))

    grid = gridmerge(grid1, grid2, grid3)
    #
    ip_iga = IGAInterpolation{RefQuadrilateral,2}()
    ip_fem = Lagrange{RefQuadrilateral,2}()
    ip_fem2 = Lagrange{RefTriangle,2}()
    dh = DofHandler(grid)
    sdh1 = SubDofHandler(dh, 1:4)
    add!(sdh1, :u, ip_iga^2)
    sdh2 = SubDofHandler(dh, 5:8)
    add!(sdh2, :u, ip_fem^2)
    add!(sdh2, :p, ip_fem)
    sdh3 = SubDofHandler(dh, 9:10)
    add!(sdh3, :u, ip_fem2^2)
    add!(sdh3, :p, ip_fem2)
    close!(dh)

    urange = dof_range(sdh1, :u)
    prange = dof_range(sdh2, :p)
    a = zeros(Float64, ndofs(dh))

    for cellid in sdh1.cellset
        dofs = celldofs(dh, cellid)
        a[dofs[urange]] .= 2.0
    end

    for cellid in sdh2.cellset
        dofs = celldofs(dh, cellid)
        a[dofs[prange]] .= 1.0
    end

    #For mixed meshes, we need to use the VTKIGAFile to export both the FEM and IGA elements
    mktempdir() do dir
        @test_nowarn VTKIGAFile("test_export2", grid) do vtk
            write_solution(vtk, dh, a)
        end
    end
end