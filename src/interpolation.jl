function barycentricD2N3(ξ::Vector{Float64})
    N = [1.0 - ξ[1] - ξ[2], ξ[1], ξ[2]]
    dN = [
        -1 1 0
        -1 0 1
    ] / 1.0
    ddN = zeros(2, 2, 3)
    return N, dN, ddN
end

function barycentricD2N3G1()
    w = 0.5 * ones(1)
    N = zeros(3, 1)
    dN = zeros(2, 3, 1)
    ξ = ones(2) / 3
    N, dN[:, :, 1], _ = barycentricD2N3(ξ)
    return N, dN, w, ξ
end

function barycentricD2N3G3()
    w = ones(3) / 6.0
    N = zeros(3, 3)
    dN = zeros(2, 3, 3)
    ξ = [
        1 4 1
        1 1 4
    ] / 6
    for p ∈ 1:3
        N[:, p], dN[:, :, p], _ = barycentricD2N3(ξ[:, p])
    end
    return N, dN, w, ξ
end

function barycentricD3N4(ξ::Vector{Float64})
    N = [1 - ξ[1] - ξ[2] - ξ[3], ξ[1], ξ[2], ξ[3]]
    dN = [
        -1 1 0 0
        -1 0 1 0
        -1 0 0 1
    ] / 1.0
    ddN = zeros(3, 3, 4)
    return N, dN, ddN
end

function barycentricD3N10(ξ::Vector{Float64})
    t0 = 1 - ξ[1] - ξ[2] - ξ[3]
    t1 = ξ[1]
    t2 = ξ[2]
    t3 = ξ[3]
    N = [
        t0 * (2 * t0 - 1),
        t1 * (2 * t1 - 1),
        t2 * (2 * t2 - 1),
        t3 * (2 * t3 - 1),
        4 * t0 * t1,
        4 * t1 * t2,
        4 * t2 * t0,
        4 * t0 * t3,
        4 * t1 * t3,
        4 * t2 * t3,
    ]
    dN = [
        1-4*t0 4*t1-1 0 0 4*(t0-t1) 4*t2 -4*t2 -4*t3 4*t3 0
        1-4*t0 0 4*t2-1 0 -4*t1 4*t1 4*(t0-t2) -4*t3 0 4*t3
        1-4*t0 0 0 4*t3-1 -4*t1 0 -4*t2 4*(t0-t3) 4*t1 4*t2
    ]
    ddN = zeros(3, 3, 10)
    ddN[1, 1, :] = [4 4 0 0 -8 0 0 0 0 0] * 1.0
    ddN[1, 2, :] = [4 0 0 0 -4 4 -4 0 0 0] * 1.0
    ddN[1, 3, :] = [4 0 0 0 -4 0 0 -4 4 0] * 1.0
    ddN[2, 1, :] = [4 0 0 0 -4 4 -4 0 0 0] * 1.0
    ddN[2, 2, :] = [4 0 4 0 0 0 -8 0 0 0] * 1.0
    ddN[2, 3, :] = [4 0 0 0 0 0 -4 -4 0 4] * 1.0
    ddN[3, 1, :] = [4 0 0 0 -4 0 0 -4 4 0] * 1.0
    ddN[3, 2, :] = [4 0 0 0 0 0 -4 -4 0 4] * 1.0
    ddN[3, 3, :] = [4 0 0 4 0 0 0 -8 0 0] * 1.0
    return N, dN, ddN
end

function barycentricD3N4G1()
    w = ones(1) / 6.0
    N = zeros(4, 1)
    dN = zeros(3, 4, 1)
    ξ = 0.25 * ones(3)
    N, dN[:, :, 1], _ = barycentricD3N4(ξ)
    return N, dN, w, ξ
end

function barycentricD3N4G4()
    w = ones(4) / 24.0
    N = zeros(4, 4)
    dN = zeros(3, 4, 4)
    s = sqrt(5.0)
    a = 5.0 + 3.0 * s
    b = 5.0 - s
    ξ = [
        b a b b
        b b a b
        b b b a
    ] / 20.0
    for p ∈ 1:4
        N[:, p], dN[:, :, p], _ = barycentricD3N4(ξ[:, p])
    end
    return N, dN, w, ξ
end

function barycentricD3N10G4()
    w = ones(4) / 24.0
    N = zeros(10, 4)
    dN = zeros(3, 10, 4)
    s = sqrt(5.0)
    a = 5.0 + 3.0 * s
    b = 5.0 - s
    ξ = [
        b a b b
        b b a b
        b b b a
    ] / 20.0
    for p ∈ 1:4
        N[:, p], dN[:, :, p], _ = barycentricD3N10(ξ[:, p])
    end
    return N, dN, w, ξ
end

function barycentricD3N10G5()
    a = -2 / 15
    b = 3 / 40
    w = [a b b b b]
    N = zeros(10, 5)
    dN = zeros(3, 10, 5)
    ξ = [
        1/4 1/6 1/6 1/6 1/2
        1/4 1/6 1/6 1/2 1/6
        1/4 1/6 1/2 1/6 1/6
    ]
    for p ∈ 1:5
        N[:, p], dN[:, :, p], _ = barycentricD3N4(ξ[:, p])
    end
    return N, dN, w, ξ
end

# Computes the nodes x and weights w
# for n-point Gauss-Legendre quadrature.
# Reference:
# G. H. Golub and J. H. Welsch, Calculation of Gauss quadrature
# rules, Math. Comp., 23(106):221-230, 1969.
function gauss_legendre(n::Integer)
    if n == 1
        return zeros(1), 2.0 * ones(1)
    elseif n == 2
        g = sqrt(3.0 / 3.0)
        return [-g, g], ones(2)
    elseif n == 3
        w = [5.0, 8.0, 5.0] / 9.0
        g = sqrt(3.0 / 5.0)
        ξ = [-g, 0, g]
        return ξ, w
    elseif n == 4
        a = sqrt(3.0 / 7.0 + 2.0 * sqrt(6.0 / 5.0) / 7.0)
        b = sqrt(3.0 / 7.0 - 2.0 * sqrt(6.0 / 5.0) / 7.0)
        c = (18.0 - sqrt(30.0)) / 36.0
        d = (18.0 + sqrt(30.0)) / 36.0
        w = [c, d, d, c]
        ξ = [-a, -b, b, a]
        return ξ, w
    end
    i = 1:n-1
    v = i ./ sqrt.(4.0 .* i .* i .- 1.0)
    vv = eigen(diagm(1 => v, -1 => v))
    ξ = vv.values
    w = 2.0 * vv.vectors[1, :] .^ 2
    return ξ, w
end

function gauss_legendreD1(n::Integer)
    return gauss_legendre(n)
end

function gauss_legendreD2(n::Integer)
    if n ∉ [1, 4, 9]
        error("Order must be in [1,4,9] : ", n)
    end
    if n == 1
        return zeros(2, 1), 4.0 * ones(1)
    elseif n == 4
        w = ones(4)
        g = sqrt(3.0) / 3.0
        ξ = g * [
            -1 1 1 -1
            -1 -1 1 1
        ]
        return ξ, w
    elseif n == 9
        x, ω = gauss_legendreD1(3)
        ξ = [
            x[1] x[3] x[3] x[1] x[2] x[3] x[2] x[1] x[2]
            x[1] x[1] x[3] x[3] x[1] x[2] x[3] x[2] x[2]
        ]
        w = [
            ω[1] * ω[1],
            ω[3] * ω[1],
            ω[3] * ω[3],
            ω[1] * ω[3],
            ω[2] * ω[1],
            ω[3] * ω[2],
            ω[2] * ω[3],
            ω[1] * ω[2],
            ω[2] * ω[2],
        ]
    end
end

function lagrangianD1N2(ξ::Float64)
    N = [0.5 * (1.0 - ξ), 0.5 * (1.0 + ξ)]
    dN = [-0.5, 0.5]
    ddN = zeros(1, 1, 2)
    return N, dN, ddN
end

function lagrangianD1N2G1()
    N = zeros(2, 1)
    dN = zeros(1, 2, 1)
    ξ, w = gauss_legendreD1(1)
    N, dN[:, :, 1], _ = lagrangianD1N2(ξ[1])
    return N, dN, w, ξ
end

function lagrangianD1N2G2()
    N = zeros(2, 2)
    dN = zeros(1, 2, 2)
    ξ, w = gauss_legendreD1(2)
    for p ∈ 1:2
        N[:, p], dN[:, :, p], _ = lagrangianD1N2(ξ[p])
    end
    return N, dN, w, ξ
end

function lagrangianD1N2G3()
    N = zeros(2, 3)
    dN = zeros(1, 2, 3)
    ξ, w = gauss_legendreD1(3)
    for p ∈ 1:3
        N[:, p], dN[:, :, p], _ = lagrangianD1N2(ξ[p])
    end
    return N, dN, w, ξ
end

function lagrangianD1N2G4()
    N = zeros(2, 4)
    dN = zeros(1, 2, 4)
    ξ, w = gauss_legendreD1(4)
    for p ∈ 1:4
        N[:, p], dN[:, :, p], _ = lagrangianD1N2(ξ[p])
    end
    return N, dN, w, ξ
end

function lagrangianD2N4(ξ::Vector{Float64})
    r = ξ[1]
    s = ξ[2]
    ra = [-1 1 1 -1] / 1.0
    sa = [-1 -1 1 1] / 1.0
    N = zeros(4)
    dN = zeros(2, 4)
    ddN = zeros(2, 2, 4)
    for p ∈ 1:4
        N[p] = 0.25 * (1.0 + ra[p] * r) * (1.0 + sa[p] * s)
        dN[1, p] = 0.25 * ra[p] * (1 + sa[p] * s)
        dN[2, p] = 0.25 * (1 + ra[p] * r) * sa[p]
        ddN[1, 1, p] = 0.0
        ddN[1, 2, p] = 0.25 * ra[p] * sa[p]
        ddN[2, 1, p] = 0.25 * ra[p] * sa[p]
        ddN[2, 2, p] = 0.0
    end
    return N, dN, ddN
end

function lagrangianD2N4G4()
    N = zeros(4, 4)
    dN = zeros(2, 4, 4)
    ξ, w = gauss_legendreD2(4)
    for p ∈ 1:4
        N[:, p], dN[:, :, p], _ = lagrangianD2N4(ξ[:, p])
    end
    return N, dN, w, ξ
end

function lagrangianD2N4G9()
    N = zeros(4, 9)
    dN = zeros(2, 4, 9)
    ξ, w = gauss_legendreD2(9)
    for p ∈ 1:4
        N[:, p], dN[:, :, p], _ = lagrangianD2N4(ξ[:, p])
    end
    return N, dN, w, ξ
end

function lagrangianD3N8(ξ::Vector{Float64})
    r = ξ[1]
    s = ξ[2]
    t = ξ[3]
    ra = [-1 1 1 -1 -1 1 1 -1] / 1.0
    sa = [-1 -1 1 1 -1 -1 1 1] / 1.0
    ta = [-1 -1 -1 -1 1 1 1 1] / 1.0
    N = zeros(8)
    dN = zeros(3, 8)
    ddN = zeros(3, 3, 8)
    for p ∈ 1:8
        N[p] = 0.125 * (1.0 + ra[p] * r) * (1.0 + sa[p] * s) * (1.0 + ta[p] * t)
        dN[1, p] = 0.125 * ra[p] * (1.0 + sa[p] * s) * (1.0 + ta[p] * t)
        dN[2, p] = 0.125 * (1.0 + ra[p] * r) * sa[p] * (1.0 + ta[p] * t)
        dN[3, p] = 0.125 * (1.0 + ra[p] * r) * (1.0 + sa[p] * s) * ta[p]
        ddN[1, 1, p] = 0.0
        ddN[1, 2, p] = 0.125 * ra[p] * sa[p] * (1.0 + ta[p] * t)
        ddN[1, 3, p] = 0.125 * ra[p] * ta[p] * (1.0 + sa[p] * s)
        ddN[2, 1, p] = 0.125 * ra[p] * sa[p] * (1.0 + ta[p] * t)
        ddN[2, 2, p] = 0.0
        ddN[2, 3, p] = 0.125 * sa[p] * ta[p] * (1.0 + ra[p] * r)
        ddN[3, 1, p] = 0.125 * ra[p] * ta[p] * (1.0 + sa[p] * s)
        ddN[3, 2, p] = 0.125 * sa[p] * ta[p] * (1.0 + ra[p] * r)
        ddN[3, 3, p] = 0.0
    end
    return N, dN, ddN
end

function lagrangianD3N8G8()
    w = ones(8)
    N = zeros(8, 8)
    dN = zeros(3, 8, 8)
    g = sqrt(3.0) / 3.0
    ξ = g * [
        -1 1 1 -1 -1 1 1 -1
        -1 -1 1 1 -1 -1 1 1
        -1 -1 -1 -1 1 1 1 1
    ]
    for p ∈ 1:8
        N[:, p], dN[:, :, p], _ = lagrangianD3N8(ξ[:, p])
    end
    return N, dN, w, ξ
end

function default_num_int_pts(element_type::String)
    if uppercase(element_type) == "BAR2"
        return 1
    elseif uppercase(element_type) == "TRI3"
        return 3
    elseif uppercase(element_type) == "QUAD4"
        return 4
    elseif uppercase(element_type) == "TETRA4" || uppercase(element_type) == "TETRA"
        return 4
    elseif uppercase(element_type) == "TETRA10"
        return 4
    elseif uppercase(element_type) == "HEX8" || uppercase(element_type) == "HEX"
        return 8
    else
        error("Invalid element type: ", element_type)
    end
end

function get_element_type(dim::Integer, num_nodes::Integer)
    if dim == 1 && num_nodes == 2
        return "BAR2"
    elseif dim == 2 && num_nodes == 3
        return "TRI3"
    elseif dim == 2 && num_nodes == 4
        return "QUAD4"
    elseif dim == 3 && num_nodes == 4
        return "TETRA4"
    elseif dim == 3 && num_nodes == 10
        return "TETRA10"
    elseif dim == 3 && num_nodes == 8
        return "HEX8"
    else
        error("Invalid dimension : ", dim, " and number of nodes : ", num_nodes)
    end
end

#
# Compute isoparametric interpolation functions, their parametric
# derivatives, integration weights, and integration point locations
#
function isoparametric(element_type::String, num_int::Integer)
    msg1 = "Invalid number of integration points: "
    msg2 = " for element type: "
    if uppercase(element_type) == "BAR2"
        if num_int == 1
            return lagrangianD1N2G1()
        elseif num_int == 2
            return lagrangianD1N2G2()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif uppercase(element_type) == "TRI3"
        if num_int == 1
            return barycentricD2N3G1()
        elseif num_int == 3
            return barycentricD2N3G3()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif uppercase(element_type) == "QUAD4"
        if num_int == 4
            return lagrangianD2N4G4()
        elseif num_int == 9
            return lagrangianD2N4G9()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif uppercase(element_type) == "TETRA4" || uppercase(element_type) == "TETRA"
        if num_int == 1
            return barycentricD3N4G1()
        elseif num_int == 4
            return barycentricD3N4G4()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif uppercase(element_type) == "TETRA10"
        if num_int == 4
            return barycentricD3N10G4()
        elseif num_int == 5
            return barycentricD3N10G5()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif uppercase(element_type) == "HEX8" || uppercase(element_type) == "HEX"
        if num_int == 8
            return lagrangianD3N8G8()
        else
            error(msg1, num_int, msg2, element_type)
        end
    else
        error("Invalid element type: ", element_type)
    end
end

function gradient_operator(dNdX::Matrix{Float64})
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

using Symbolics
@variables t, x, y, z

function get_side_set_nodal_forces(
    nodal_coord::Matrix{Float64},
    traction_num::Num,
    time::Float64
)
    _, num_side_nodes = size(nodal_coord)
    element_type = get_element_type(2, num_side_nodes)
    num_int_points = default_num_int_pts(element_type)
    N, dNdξ, w, _ = isoparametric(element_type, num_int_points)
    nodal_force_component = zeros(num_side_nodes)
    for point ∈ 1:num_int_points
        Nₚ = N[:, point]
        dNdξₚ = dNdξ[:, :, point]
        dXdξ = dNdξₚ * nodal_coord'
        j = norm(cross(dXdξ[1, :], dXdξ[2, :]))
        wₚ = w[point]
        point_coord = nodal_coord * Nₚ
        values =
            Dict(t => time, x => point_coord[1], y => point_coord[2], z => point_coord[3])
        traction_sym = substitute(traction_num, values)
        traction_val = extract_value(traction_sym)
        nodal_force_component += traction_val * Nₚ * j * wₚ
    end
    return nodal_force_component
end

# Helper function to calculate the signed volume of a tetrahedron
function signed_volume(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64}, p4::Vector{Float64})
    return det([p2 - p1 p3 - p1 p4 - p1])
end

# Check if point is inside: all volumes should have the same sign as V_total within tolerance
function same_sign_or_within_tolerance(v1::Float64, v2::Float64, tol::Float64)
    return abs(v1 - v2) ≤ tol || sign(v1) == sign(v2)
end

function is_point_in_tetrahedron(P::Vector{Float64}, A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64}, D::Vector{Float64}, tol::Float64)

    # Calculate signed volumes
    V_total = signed_volume(A, B, C, D)
    V1 = signed_volume(P, B, C, D)
    V2 = signed_volume(A, P, C, D)
    V3 = signed_volume(A, B, P, D)
    V4 = signed_volume(A, B, C, P)

    return (same_sign_or_within_tolerance(V1, V_total, tol) &&
            same_sign_or_within_tolerance(V2, V_total, tol) &&
            same_sign_or_within_tolerance(V3, V_total, tol) &&
            same_sign_or_within_tolerance(V4, V_total, tol))
end

function is_point_in_hexahedron(P::Vector{Float64}, A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64}, D::Vector{Float64}, E::Vector{Float64}, F::Vector{Float64}, G::Vector{Float64}, H::Vector{Float64}, tol::Float64)

    # Check the point against each tetrahedron
    return (
        is_point_in_tetrahedron(P, A, B, C, E, tol) ||
        is_point_in_tetrahedron(P, C, B, F, E, tol) ||
        is_point_in_tetrahedron(P, C, F, G, E, tol) ||
        is_point_in_tetrahedron(P, C, G, H, E, tol) ||
        is_point_in_tetrahedron(P, C, H, D, E, tol)
    )
end

function map_to_parametric(
    element_type::String,
    nodes::Matrix{Float64},
    point::Vector{Float64}
)
    tol = 1.0e-08
    max_iters = 32
    dim = length(point)
    ξ = zeros(dim)
    hessian = zeros(dim, dim)
    for _ ∈ 1:max_iters
        N, dN, _ = interpolate(element_type, ξ)
        trial_point = nodes * N
        residual = trial_point - point
        hessian = nodes * dN'
        δ = -hessian \ residual
        ξ = ξ + δ
        error = norm(δ)
        if error <= tol
            break
        end
    end
    return ξ
end

function interpolate(element_type::String, ξ::Vector{Float64})
    if uppercase(element_type) == "BAR2"
        return lagrangianD1N2(ξ)
    elseif uppercase(element_type) == "TRI3"
        return barycentricD2N3(ξ)
    elseif uppercase(element_type) == "QUAD4"
        return lagrangianD2N4(ξ)
    elseif uppercase(element_type) == "TETRA4" || uppercase(element_type) == "TETRA"
        return barycentricD3N4(ξ)
    elseif uppercase(element_type) == "TETRA10"
        return barycentricD3N10(ξ)
    elseif uppercase(element_type) == "HEX8" || uppercase(element_type) == "HEX"
        return lagrangianD3N8(ξ)
    else
        error("Invalid element type: ", element_type)
    end
end

function is_inside_parametric(element_type::String, ξ::Vector{Float64}, tol::Float64 = 1.0e-06)
    factor = 1.0 + tol
    if uppercase(element_type) == "BAR2"
        return -factor ≤ ξ ≤ factor
    elseif uppercase(element_type) == "TRI3"
        return reduce(*, -tol * ones(2) .≤ ξ .≤ factor * ones(2))
    elseif uppercase(element_type) == "QUAD4"
        return reduce(*, -factor * ones(2) .≤ ξ .≤ factor * ones(2))
    elseif uppercase(element_type) == "TETRA4" || uppercase(element_type) == "TETRA10" || uppercase(element_type) == "TETRA"
        return reduce(*, -tol * ones(3) .≤ ξ .≤ factor * ones(3))
    elseif uppercase(element_type) == "HEX8" || uppercase(element_type) == "HEX"
        return reduce(*, -factor * ones(3) .≤ ξ .≤ factor * ones(3))
    else
        error("Invalid element type: ", element_type)
    end
end

function is_inside(element_type::String, nodes::Matrix{Float64}, point::Vector{Float64}, tol::Float64 = 1.0e-06)
    ξ = map_to_parametric(element_type, nodes, point)
    return is_inside_parametric(element_type, ξ, tol)
end

#function is_inside(element_type::String, nodes::Matrix{Float64}, point::Vector{Float64}, tol::Float64 = 1.0e-06)
#    if uppercase(element_type) == "TETRA4" || uppercase(element_type) == "TETRA10" || uppercase(element_type) == "TETRA"
#        return is_point_in_tetrahedron(point, nodes[:,1], nodes[:,2], nodes[:,3], nodes[:,4], tol)
#    elseif uppercase(element_type) == "HEX8"
#        return is_point_in_hexahedron(point, nodes[:,1], nodes[:,2], nodes[:,3], nodes[:,4], nodes[:,5], nodes[:,6], nodes[:,7], nodes[:,8], tol)
#    else
#        error("Invalid element type: ", element_type)
#    end
#end

function find_and_project(
    point::Vector{Float64},
    mesh::ExodusDatabase,
    side_set_id::Integer,
    model::SolidMechanics,
)
    #we assume that we know the contact surfaces in advance 
    num_nodes_per_sides, side_set_node_indices =
        Exodus.read_side_set_node_list(mesh, side_set_id)
    ss_node_index = 1
    point_new = point
    closest_face_nodes = Array{Float64}(undef, 0)
    closest_face_node_indices = Array{Int64}(undef, 0)
    space_dim = length(point)
    parametric_dim = space_dim - 1
    closest_ξ = zeros(parametric_dim)
    closest_normal = zeros(space_dim)
    minimum_absolute_distance = Inf
    closest_distance = 0.0
    for num_nodes_side ∈ num_nodes_per_sides
        face_node_indices =
            side_set_node_indices[ss_node_index:ss_node_index+num_nodes_side-1]
        face_nodes = model.current[:, face_node_indices]
        trial_point, ξ, distance, normal =
            closest_point_projection(parametric_dim, face_nodes, point)
        distance_centr = get_distance_to_centroid(face_nodes, point)
        if abs(distance_centr) < minimum_absolute_distance
            minimum_absolute_distance = abs(distance_centr)
            point_new = trial_point
            closest_face_nodes = face_nodes
            closest_face_node_indices = face_node_indices
            closest_normal = normal
            closest_ξ = ξ
            closest_distance = distance
        end
        ss_node_index += num_nodes_side
    end
    return point_new,
    closest_ξ,
    closest_face_nodes,
    closest_face_node_indices,
    closest_normal,
    closest_distance
end

function get_distance_to_centroid(nodes::Matrix{Float64}, x::Vector{Float64})
    _, num_nodes = size(nodes)
    x_coords = 0.0
    y_coords = 0.0
    z_coords = 0.0
    for i ∈ 1:num_nodes
        x_coords = x_coords + nodes[1, i]
        y_coords = y_coords + nodes[2, i]
        z_coords = z_coords + nodes[3, i]
    end
    x_centroid = x_coords / num_nodes
    y_centroid = y_coords / num_nodes
    z_centroid = z_coords / num_nodes
    centroid = [x_centroid, y_centroid, z_centroid]
    distance = norm(centroid - x)
    return distance
end

function get_side_set_global_to_local_map(mesh::ExodusDatabase, side_set_id::Integer)
    num_nodes_per_sides, side_set_node_indices =
        Exodus.read_side_set_node_list(mesh, side_set_id)
    unique_node_indices = unique(side_set_node_indices)
    num_nodes = length(unique_node_indices)
    global_to_local_map = Dict{Int64,Int64}()
    for i ∈ 1:num_nodes
        global_to_local_map[Int64(unique_node_indices[i])] = i
    end
    return global_to_local_map, num_nodes_per_sides, side_set_node_indices
end

function get_side_set_local_to_global_map(mesh::ExodusDatabase, side_set_id::Integer)
    side_set_node_indices = Exodus.read_side_set_node_list(mesh, side_set_id)[2]
    unique_node_indices = unique(side_set_node_indices)
    num_nodes = length(unique_node_indices)
    local_to_global_map = zeros(Int64, num_nodes)
    for i ∈ 1:num_nodes
        local_to_global_map[i] = Int64(unique_node_indices[i])
    end
    return local_to_global_map
end

function get_square_projection_matrix(
    mesh::ExodusDatabase,
    model::SolidMechanics,
    side_set_id::Integer,
)
    global_to_local_map, num_nodes_sides, side_set_node_indices =
        get_side_set_global_to_local_map(mesh, side_set_id)
    num_nodes = length(global_to_local_map)
    coords = model.current
    square_projection_matrix = zeros(num_nodes, num_nodes)
    side_set_node_index = 1
    for num_nodes_side ∈ num_nodes_sides
        side_nodes =
            side_set_node_indices[side_set_node_index:side_set_node_index+num_nodes_side-1]
        side_coordinates = coords[:, side_nodes]
        element_type = get_element_type(2, Int64(num_nodes_side))
        num_int_points = default_num_int_pts(element_type)
        N, dNdξ, w, _ = isoparametric(element_type, num_int_points)
        side_matrix = zeros(num_nodes_side, num_nodes_side)
        for point ∈ 1:num_int_points
            Nₚ = N[:, point]
            dNdξₚ = dNdξ[:, :, point]
            dXdξ = dNdξₚ * side_coordinates'
            j = norm(cross(dXdξ[1, :], dXdξ[2, :]))
            wₚ = w[point]
            side_matrix += Nₚ * Nₚ' * j * wₚ
        end
        local_indices = get.(Ref(global_to_local_map), side_nodes, 0)
        square_projection_matrix[local_indices, local_indices] += side_matrix
        side_set_node_index += num_nodes_side
    end
    return square_projection_matrix
end

function get_rectangular_projection_matrix(
    dst_mesh::ExodusDatabase,
    dst_model::SolidMechanics,
    dst_side_set_id::Integer,
    src_mesh::ExodusDatabase,
    src_model::SolidMechanics,
    src_side_set_id::Integer,
)
    dst_global_to_local_map, dst_num_nodes_sides, dst_side_set_node_indices =
        get_side_set_global_to_local_map(dst_mesh, dst_side_set_id)
    dst_num_nodes = length(dst_global_to_local_map)
    dst_coords = dst_model.current
    src_global_to_local_map, _, _ =
        get_side_set_global_to_local_map(src_mesh, src_side_set_id)
    src_num_nodes = length(src_global_to_local_map)
    rectangular_projection_matrix = zeros(dst_num_nodes, src_num_nodes)
    dst_local_indices = Array{Int64}(undef, 0)
    src_local_indices = Array{Int64}(undef, 0)
    dst_side_set_node_index = 1
    for dst_num_nodes_side ∈ dst_num_nodes_sides
        dst_side_nodes =
            dst_side_set_node_indices[dst_side_set_node_index:dst_side_set_node_index+dst_num_nodes_side-1]
        dst_local_indices = get.(Ref(dst_global_to_local_map), dst_side_nodes, 0)
        dst_side_coordinates = dst_coords[:, dst_side_nodes]
        dst_element_type = get_element_type(2, Int64(dst_num_nodes_side))
        dst_num_int_points = default_num_int_pts(dst_element_type)
        dst_N, dst_dNdξ, dst_w, _ = isoparametric(dst_element_type, dst_num_int_points)
        for dst_point ∈ 1:dst_num_int_points
            dst_Nₚ = dst_N[:, dst_point]
            dst_dNdξₚ = dst_dNdξ[:, :, dst_point]
            dst_dXdξ = dst_dNdξₚ * dst_side_coordinates'
            dst_j = norm(cross(dst_dXdξ[1, :], dst_dXdξ[2, :]))
            dst_wₚ = dst_w[dst_point]
            dst_int_point_coord = dst_side_coordinates * dst_Nₚ
            is_inside = false
            _, ξ, src_side_coordinates, src_side_nodes, _, _ =
                find_and_project(dst_int_point_coord, src_mesh, src_side_set_id, src_model)
            src_side_element_type = get_element_type(2, size(src_side_coordinates)[2])
            src_Nₚ, _, _ = interpolate(src_side_element_type, ξ)
            src_local_indices = get.(Ref(src_global_to_local_map), src_side_nodes, 0)
            rectangular_projection_matrix[dst_local_indices, src_local_indices] +=
                dst_Nₚ * src_Nₚ' * dst_j * dst_wₚ
        end
        dst_side_set_node_index += dst_num_nodes_side
    end
    return rectangular_projection_matrix
end

function compute_normal(mesh::ExodusDatabase, side_set_id::Int64, model::SolidMechanics)
    global_to_local_map, num_nodes_sides, side_set_node_indices =
        get_side_set_global_to_local_map(mesh, side_set_id)
    coords = model.current
    num_nodes = length(global_to_local_map)
    space_dim, _ = size(coords)
    normals = zeros(space_dim, num_nodes)
    local_indices = Array{Int64}(undef, 0)
    side_set_node_index = 1
    for num_nodes_side ∈ num_nodes_sides
        side_nodes =
            side_set_node_indices[side_set_node_index:side_set_node_index+num_nodes_side-1]
        local_indices = get.(Ref(global_to_local_map), side_nodes, 0)
        coordinates = coords[:, side_nodes]
        point_A = coordinates[:, 1]
        point_B = coordinates[:, 2]
        point_C = coordinates[:, end]
        BA = point_B - point_A
        CA = point_C - point_A
        N = cross(BA, CA)
        normals[:, local_indices] .= N / norm(N)
        side_set_node_index += num_nodes_side
    end
    return normals
end

function interpolate(
    tᵃ::Float64,
    tᵇ::Float64,
    xᵃ::Vector{Float64},
    xᵇ::Vector{Float64},
    t::Float64,
)
    Δt = tᵇ - tᵃ
    if Δt == 0.0
        return 0.5 * (xᵃ + xᵇ)
    end
    p = (tᵇ - t) / Δt
    q = 1.0 - p
    return p * xᵃ + q * xᵇ
end

function interpolate(
    param_hist::Vector{Float64},
    value_hist::Vector{Vector{Float64}},
    param::Float64,
)
    if param < param_hist[1]
        param = param_hist[1]
    end
    if param > param_hist[end]
        param = param_hist[end]
    end
    index = 1
    size = length(param_hist)
    while param_hist[index] < param
        if index == size
            break
        end
        index += 1
    end
    if index == 1
        return value_hist[1]
    elseif index == size
        return value_hist[size]
    else
        return interpolate(
            param_hist[index],
            param_hist[index+1],
            value_hist[index],
            value_hist[index+1],
            param,
        )
    end
end

using Einsum

function closest_point_projection(
    parametric_dim::Integer,
    nodes::Matrix{Float64},
    x::Vector{Float64},
)
    space_dim, num_nodes = size(nodes)
    element_type = get_element_type(parametric_dim, num_nodes)
    ξ = zeros(parametric_dim)
    residual = zeros(parametric_dim)
    hessian = zeros(parametric_dim, parametric_dim)
    y = x
    yx = zeros(space_dim)
    tol = 1.0e-12
    normal = zeros(space_dim)
    iteration = 1
    max_iterations = 64
    while true
        N, dN, ddN = interpolate(element_type, ξ)
        y = nodes * N
        dydξ = dN * nodes'
        yx = y - x
        residual = dydξ * yx
        @einsum ddyddξ[i, j, k] := ddN[i, j, l] * nodes[k, l]
        @einsum ddyddξyx[i, j] := ddyddξ[i, j, k] * yx[k]
        hessian = ddyddξyx + dydξ * dydξ'
        δ = -hessian \ residual
        ξ = ξ + δ
        error = norm(δ)
        if error <= tol
            break
        end
        iteration += 1
        if iteration > max_iterations
            error("Closest point projection failed to converge")
        end
    end
    _, dN, _ = interpolate(element_type, ξ)
    dxdξ = dN * nodes'
    perp_vec = cross(dxdξ[1, :], dxdξ[2, :])
    normal = perp_vec / norm(perp_vec)
    distance = -copysign(norm(yx), dot(yx, normal))
    return y, ξ, distance, normal
end
