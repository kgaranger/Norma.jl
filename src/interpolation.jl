function barycentricD2N3(ξ::Vector{Float64})
    N = [1.0 - ξ[1] - ξ[2], ξ[1], ξ[2]]
    dN = [-1 1 0
        -1 0 1] / 1.0
    return N, dN
end

function barycentricD2N3G1()
    w = 0.5 * ones(1, 1)
    N = zeros(3, 1)
    dN = zeros(2, 3, 1)
    ξ = ones(2) / 3.0
    N, dN[:, :, 1] = barycentricD2N3(ξ)
    return N, dN, w
end

function barycentricD2N3G3()
    w = ones(1, 3) / 6.0
    N = zeros(3, 3)
    dN = zeros(2, 3, 3)
    ξ = [4 1 1 1
        1 4 1 1
        1 1 4 1] / 6.0
    for i ∈ 1:3
        N[:, i], dN[:, :, i] = barycentricD2N3(ξ[:, i])
    end
    return N, dN, w
end

function barycentricD3N4(ξ::Vector{Float64})
    N = [1.0 - ξ[1] - ξ[2] - ξ[3],
        ξ[1],
        ξ[2],
        ξ[3]]
    dN = [-1 1 0 0
        -1 0 1 0
        -1 0 0 1] / 1.0
    return N, dN
end

function barycentricD3N4G1()
    w = ones(1, 1) / 6.0
    N = zeros(4, 1)
    dN = zeros(3, 4, 1)
    ξ = 0.25 * ones(3, 1)
    N, dN[:, :, 1] = barycentricD3N4(ξ)
    return N, dN, w
end

function barycentricD3N4G4()
    w = ones(1, 4) / 24.0
    N = zeros(4, 4)
    dN = zeros(3, 4, 4)
    s = sqrt(5.0)
    a = 5.0 + 3.0 * s
    b = 5.0 - s
    ξ = [a b b b
        b a b b
        b b a b] / 20.0
    for i ∈ 1:4
        N[:, i], dN[:, :, i] = barycentricD3N4(ξ[:, i])
    end
    return N, dN, w
end

function lagrangianD1N2(ξ::Float64)
    N = [0.5 * (1.0 - ξ), 0.5 * (1.0 + ξ)]
    dN = [-0.5, 0.5]
    return N, dN
end

function lagrangianD1N2G1()
    w = 2.0 * ones(1, 1)
    N = zeros(2, 1)
    dN = zeros(1, 2, 1)
    N, dN[:, :, 1] = lagrangianD1N2(0.0)
    return N, dN, w
end

function lagrangianD1N2G2()
    w = ones(1, 2)
    N = zeros(2, 2)
    dN = zeros(1, 2, 2)
    g = sqrt(3.0) / 3.0
    ξ = [-g, g]
    for i ∈ 1:2
        N[:, i], dN[:, :, i] = lagrangianD1N2(ξ[i])
    end
    return N, dN, w
end

function lagrangianD2N4(ξ::Vector{Float64})
    r = ξ[1]
    s = ξ[2]
    ra = [-1 1 1 -1] / 1.0
    sa = [-1 -1 1 1] / 1.0
    N = zeros(1, 4)
    dN = zeros(2, 4)
    for i ∈ 1:4
        N[i] = 0.25 * (1.0 + ra[i] * r) * (1.0 + sa[i] * s)
        dN[1, i] = 0.25 * ra[i] * (1 + sa[i] * s)
        dN[2, i] = 0.25 * (1 + ra[i] * r) * sa[i]
    end
    return N, dN
end

function lagrangianD2N4G4()
    w = ones(1, 4)
    N = zeros(4, 4)
    dN = zeros(2, 4, 4)
    g = sqrt(3.0) / 3.0
    ξ = g * [-1 1 1 -1
        -1 -1 1 1]
    for i ∈ 1:4
        N[:, i], dN[:, :, i] = lagrangianD2N4(ξ[:, i])
    end
    return N, dN, w
end

function lagrangianD3N8(ξ::Vector{Float64})
    r = ξ[1]
    s = ξ[2]
    t = ξ[3]
    ra = [-1 1 1 -1 -1 1 1 -1] / 1.0
    sa = [-1 -1 1 1 -1 -1 1 1] / 1.0
    ta = [-1 -1 -1 -1 1 1 1 1] / 1.0
    N = zeros(1, 8)
    dN = zeros(3, 8)
    for i ∈ 1:8
        N[i] = 0.125 * (1.0 + ra[i] * r) * (1.0 + sa[i] * s) * (1.0 + ta[i] * t)
        dN[1, i] = 0.125 * ra[i] * (1.0 + sa[i] * s) * (1.0 + ta[i] * t)
        dN[2, i] = 0.125 * (1.0 + ra[i] * r) * sa[i] * (1.0 + ta[i] * t)
        dN[3, i] = 0.125 * (1.0 + ra[i] * r) * (1.0 + sa[i] * s) * ta[i]
    end
    return N, dN
end

function lagrangianD3N8G8()
    w = ones(1, 8)
    N = zeros(8, 8)
    dN = zeros(3, 8, 8)
    g = sqrt(3.0) / 3.0
    ξ = g * [-1 1 1 -1 -1 1 1 -1
        -1 -1 1 1 -1 -1 1 1
        -1 -1 -1 -1 1 1 1 1]
    for i ∈ 1:8
        N[:, i], dN[:, :, i] = lagrangianD3N8(ξ[:, i])
    end
    return N, dN, w
end

function default_num_int_pts(element_type)
    if element_type == "BAR2"
        return 1
    elseif element_type == "TRI3"
        return 1
    elseif element_type == "QUAD4"
        return 4
    elseif element_type == "TETRA4"
        return 1
    elseif element_type == "HEX8"
        return 8
    else
        error("Invalid element type: ", element_type)
    end
end

#
# Compute isoparametric interpolation functions, their parametric
# derivatives and integration weights.
#
function isoparametric(element_type, num_int)
    msg1 = "Invalid number of integration points: "
    msg2 = " for element type: "
    str_int = string(num_int)
    if element_type == "BAR2"
        if num_int == 1
            N, dN, w = lagrangianD1N2G1()
        elseif num_int == 2
            N, dN, w = lagrangianD1N2G2()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "TRI3"
        if num_int == 1
            N, dN, w = barycentricD2N3G1()
        elseif num_int == 3
            N, dN, w = barycentricD2N3G3()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "QUAD4"
        if num_int == 4
            N, dN, w = lagrangianD2N4G4()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "TETRA4"
        if num_int == 1
            N, dN, w = barycentricD3N4G1()
        elseif num_int == 4
            N, dN, w = barycentricD4N4G4()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "HEX8"
        if num_int == 8
            N, dN, w = lagrangianD3N8G8()
        else
            error(msg1, num_int, msg2, element_type)
        end
    else
        error("Invalid element type: ", element_type)
    end
    return N, dN, w
end

function gradient_operator(dNdX)
    dim, nen = size(dNdX)
    B = zeros(dim * dim, nen * dim)
    for i ∈ 1:dim
        for j ∈ 1:dim
            p = dim * (i - 1) + j
            for a ∈ 1:nen
                for k ∈ 1:dim
                    q = dim * (a - 1) + k
                    B[p, q] = I[i, k] * dNdX[j, a]
                end
            end
        end
    end
    return B
end

function two_dims_from_three_dims(coordinates3D::Matrix{Float64})
    O = coordinates3D[:, 1]
    A = coordinates3D[:, 2]
    B = coordinates3D[:, 3]
    OA = A - O
    Ax = sqrt(OA' * OA)
    OB = B - O
    lengthOB2 = OB' * OB
    n = OA / Ax
    Bx = OB' * n
    By = sqrt(lengthOB2 - Bx * Bx)
    coordinates2D = zeros(2, 3)
    coordinates2D[1, 2] = Ax
    coordinates2D[1, 3] = Bx
    coordinates2D[2, 3] = By
    return coordinates2D
end

function get_side_set_nodal_forces(coordinates::Matrix{Float64}, expr::Any, time::Float64)
    _, num_nodes = size(coordinates)
    if num_nodes == 3
        return get_triangle_nodal_forces(coordinates, expr, time)
    elseif num_nodes == 4
        indices1 = [1,2,3]
        indices2 = [1,3,4]
        triangle1 = coordinates[:, indices1]
        forces1 = get_triangle_nodal_forces(triangle1, expr, time)
        triangle2 = coordinates[:, indices2]
        forces2 = get_triangle_nodal_forces(triangle2, expr, time)
        forces = zeros(4)
        forces[indices1] = forces1
        forces[indices2] += forces2
        return forces
    else
        error("Unknown side topology with number of nodes: ", num_nodes)
    end
end

function get_triangle_nodal_forces(coordinates::Matrix{Float64}, expr::Any, time::Float64)
    centroid = (coordinates[:, 1] + coordinates[:, 2] + coordinates[:, 3]) / 3.0
    coordinates2D = two_dims_from_three_dims(coordinates)
    Nₚ, dNdξ, elem_weights = barycentricD2N3G1()
    point = 1
    dNdξₚ = dNdξ[:, :, point]
    dXdξ = dNdξₚ * coordinates2D'
    j = det(dXdξ)
    w = elem_weights[point]
    global t = time
    global x = centroid[1]
    global y = centroid[2]
    global z = centroid[3]
    traction = eval(expr)
    nodal_force_component = traction * Nₚ * j * w
    return nodal_force_component
end