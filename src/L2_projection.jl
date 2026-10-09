#assembles the L2 matrix for BezierCellValues
function Ferrite._assemble_L2_matrix!(assembler, cellvalues::BezierCellValues, sdh::Ferrite.SubDofHandler)
    grid = Ferrite.get_grid(sdh.dh)
    grid isa BezierGrid || error("BezierCellValues L2 projection requires a BezierGrid")
    #reuse the dof buffer across cells
    bcoords = getcoordinates(grid, first(sdh.cellset))
    dofs = zeros(Int, ndofs_per_cell(sdh.dh, first(sdh.cellset)))
    Me = zeros(0, 0)

    for cellid in sdh.cellset
        #updates the current cell's extraction operator and transforms the Bernstein polynomials to splines
        getcoordinates!(bcoords, grid, cellid)
        reinit!(cellvalues, bcoords)
        n = getnbasefunctions(cellvalues)
        size(Me, 1) == n || (Me = zeros(n, n))
        fill!(Me, 0)
        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            for j in 1:n
                v = shape_value(cellvalues, q_point, j)
                for i in 1:j
                    u = shape_value(cellvalues, q_point, i)
                    Me[i, j] += v ⋅ u * dΩ
                end
            end
        end

        _symmetrize_to_lower!(Me)
        celldofs!(dofs, sdh.dh, cellid)
        assemble!(assembler, dofs, Me)
    end
    return assembler
end

#rhs of the l2 projection
function Ferrite.assemble_proj_rhs!(
    f::Matrix, cellvalues::BezierCellValues, sdh::Ferrite.SubDofHandler,
    vars::Union{AbstractVector, AbstractDict})

    grid = Ferrite.get_grid(sdh.dh)
    grid isa BezierGrid || error("BezierCellValues L2 projection requires a BezierGrid")

    M = size(f, 2)
    nqp = getnquadpoints(cellvalues)
    bcoords = getcoordinates(grid, first(sdh.cellset))
    dofs = zeros(Int, ndofs_per_cell(sdh.dh, first(sdh.cellset)))
    fe = zeros(0, M)

    for cellid in sdh.cellset
        cell_vars = vars[cellid]
        length(cell_vars) == nqp || error("The number of variables per cell doesn't match the number of quadrature points")

        #rhs uses the given spline basis and geometry from the BezierCoords
        getcoordinates!(bcoords, grid, cellid)
        reinit!(cellvalues, bcoords)
        n = getnbasefunctions(cellvalues)
        size(fe, 1) == n || (fe = zeros(n, M))
        fill!(fe, 0)

        for q_point in 1:nqp
            dΩ = getdetJdV(cellvalues, q_point)
            qp_vars = cell_vars[q_point]
            for i in 1:n
                v = shape_value(cellvalues, q_point, i)
                for j in 1:M
                    fe[i, j] += v * _projection_component(qp_vars, j) * dΩ
                end
            end
        end

        #scatters local to global 
        celldofs!(dofs, sdh.dh, cellid)
        for (i, dof) in pairs(dofs)
            f[dof, :] += fe[i, :]
        end
    end
    return f
end

_projection_component(x::Number, _) = x
_projection_component(x::AbstractTensor, i) = x.data[i]

#mirrors the upper to a lower trianglur matrix
function _symmetrize_to_lower!(K::Matrix)
    for i in 1:size(K, 1)
        for j in (i + 1):size(K, 1)
            K[j, i] = K[i, j]
        end
    end
    return K
end
