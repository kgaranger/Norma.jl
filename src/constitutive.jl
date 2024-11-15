function elastic_constants(params::Dict{Any,Any})
    E = 0.0
    ν = 0.0
    κ = 0.0
    λ = 0.0
    μ = 0.0
    if haskey(params, "elastic modulus") == true
        E = params["elastic modulus"]
        if haskey(params, "Poisson's ratio") == true
            ν = params["Poisson's ratio"]
            κ = E / 3(1 - 2ν)
            λ = E * ν / (1 + ν) / (1 - 2ν)
            μ = E / 2(1 + ν)
        elseif haskey(params, "bulk modulus") == true
            κ = params["bulk modulus"]
            ν = (3κ - E) / 6κ
            λ = (3κ * (3κ - E)) / (9κ - E)
            μ = 3κ * E / (9κ - E)
        elseif haskey(params, "Lamé's first constant") == true
            λ = params["Lamé's first constant"]
            R = sqrt(E^2 + 9λ^2 + 2E * λ)
            ν = 2λ / (E + λ + R)
            κ = (E + 3λ + R) / 6
            μ = (E - 3λ + R) / 4
        elseif haskey(params, "shear modulus") == true
            μ = params["shear modulus"]
            ν = E / 2μ - 1
            κ = E * μ / 3(3μ - E)
            λ = μ * (E - 2μ) / (3μ - E)
        else
            error("Two elastic constants are required but only elastic modulus found")
        end
    elseif haskey(params, "Poisson's ratio") == true
        ν = params["Poisson's ratio"]
        if haskey(params, "bulk modulus") == true
            κ = params["bulk modulus"]
            E = 3κ * (1 - 2ν)
            λ = 3κ * ν / (1 + ν)
            μ = 3κ * (1 - 2ν) / 2(1 + ν)
        elseif haskey(params, "Lamé's first constant") == true
            λ = params["Lamé's first constant"]
            E = λ * (1 + ν) * (1 - 2ν) / ν
            κ = λ * (1 + ν) / 3ν
            μ = λ * (1 - 2ν) / 2ν
        elseif haskey(params, "shear modulus") == true
            μ = params["shear modulus"]
            E = 2μ * (1 + ν)
            κ = 2μ * (1 + ν) / 3(1 - 2ν)
            λ = 2μ * ν / (1 - 2ν)
        else
            error("Two elastic constants are required but only Poisson's ratio found")
        end
    elseif haskey(params, "bulk modulus") == true
        κ = params["bulk modulus"]
        if haskey(params, "Lamé's first constant") == true
            λ = params["Lamé's first constant"]
            E = 9κ * (κ - λ) / (3κ - λ)
            ν = λ / (3κ - λ)
            μ = 3(κ - λ) / 2
        elseif haskey(params, "shear modulus") == true
            μ = params["shear modulus"]
            E = 9κ * μ / (3κ + μ)
            ν = (3κ - 2μ) / 2(3κ + μ)
            λ = κ - 2μ / 3
        else
            error("Two elastic constants are required but only bulk modulus found")
        end
    elseif haskey(params, "Lamé's first constant") == true
        λ = params["Lamé's first constant"]
        if haskey(params, "shear modulus") == true
            μ = params["shear modulus"]
            E = μ * (3λ + 2μ) / (λ + μ)
            ν = λ / 2(λ + μ)
            κ = λ + 2μ / 3
        else
            error("Two elastic constants are required but only Lamé's first constant found")
        end
    elseif haskey(params, "shear modulus") == true
        error("Two elastic constants are required but only shear modulus found")
    else
        error("Two elastic constants are required but none found")
    end
    return E, ν, κ, λ, μ
end

mutable struct SaintVenant_Kirchhoff <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function SaintVenant_Kirchhoff(params::Dict{Any,Any})
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = params["density"]
        new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct Linear_Elastic <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function Linear_Elastic(params::Dict{Any,Any})
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = params["density"]
        new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct Neohookean <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function Neohookean(params::Dict{Any,Any})
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = params["density"]
        new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct NeohookeanAD <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function NeohookeanAD(params::Dict{Any,Any})
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = params["density"]
        new(E, ν, κ, λ, μ, ρ)
    end
end

mutable struct J2 <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    Y₀::Float64
    n::Float64
    ε₀::Float64
    Sᵥᵢₛ₀::Float64
    m::Float64
    ∂ε∂t₀::Float64
    Cₚ::Float64
    β::Float64
    T₀::Float64
    Tₘ::Float64
    M::Float64
    function J2(params::Dict{Any,Any})
        E, ν, κ, λ, μ = elastic_constants(params)
        ρ = params["density"]
        Y₀ = params["yield stress"]
        n = get(params, "hardening exponent", 0.0)
        ε₀ = get(params, "reference plastic strain", 0.0)
        Sᵥᵢₛ₀ = get(params, "reference viscoplastic stress", 0.0)
        m = get(params, "rate dependence exponent", 0.0)
        ∂ε∂t₀ = get(params, "reference plastic strain rate", 0.0)
        Cₚ = get(params, "specific heat capacity", 0.0)
        β = get(params, "Taylor-Quinney coefficient", 0.0)
        T₀ = get(params, "reference temperature", 0.0)
        Tₘ = get(params, "melting temperature", 0.0)
        M = get(params, "thermal softening exponent", 0.0)
        κ = E / (1.0 - 2.0 * ν) / 3.0
        μ = E / (1.0 + ν) / 2.0
        new(E, ν, κ, λ, μ, ρ, Y₀, n, ε₀, Sᵥᵢₛ₀, m, ∂ε∂t₀, Cₚ, β, T₀, Tₘ, M)
    end
end

function temperature_multiplier(material::J2, T::Float64)
    T₀ = material.T₀
    Tₘ = material.Tₘ
    M = material.M
    M > 0.0 ? 1.0 - ((T - T₀) / (Tₘ - T₀))^M : 1.0
end

function hardening_potential(material::J2, ε::Float64)
    Y₀ = material.Y₀
    n = material.n
    ε₀ = material.ε₀
    exponent = (1.0 + n) / n
    n > 0.0 ? Y₀ * ε₀ / exponent * ((1.0 + ε / ε₀)^exponent - 1.0) : Y₀ * ε
end

function hardening_rate(material::J2, ε::Float64)
    Y₀ = material.Y₀
    n = material.n
    ε₀ = material.ε₀
    exponent = (1.0 - n) / n
    n > 0.0 ? Y₀ / ε₀ / n * (1.0 + ε / ε₀)^exponent : 0.0
end

function flow_strength(material::J2, ε::Float64)
    Y₀ = material.Y₀
    n = material.n
    ε₀ = material.ε₀
    n > 0.0 ? Y₀ * (1.0 + ε / ε₀)^(1.0 / n) : Y₀
end

function viscoplastic_dual_kinetic_potential(material::J2, Δε::Float64, Δt::Float64)
    Sᵥᵢₛ₀ = material.Sᵥᵢₛ₀
    m = material.m
    ∂ε∂t₀ = material.∂ε∂t₀
    exponent = (1.0 + m) / m
    Sᵥᵢₛ₀ > 0.0 && Δt > 0.0 && Δε > 0.0 ?
    Δt * Sᵥᵢₛ₀ * ∂ε∂t₀ / exponent * (Δε / Δt / ∂ε∂t₀)^exponent : 0.0
end

function viscoplastic_stress(material::J2, Δε::Float64, Δt::Float64)
    Sᵥᵢₛ₀ = material.Sᵥᵢₛ₀
    m = material.m
    ∂ε∂t₀ = material.∂ε∂t₀
    Sᵥᵢₛ₀ > 0.0 && Δt > 0.0 && Δε > 0.0 ?
    Sᵥᵢₛ₀ / ∂ε∂t₀ / Δt / m * (Δε / Δt / ∂ε∂t₀)^((1.0 - m) / m) : 0.0
end

function viscoplastic_hardening_rate(material::J2, Δε::Float64, Δt::Float64)
    Sᵥᵢₛ₀ = material.Sᵥᵢₛ₀
    m = material.m
    ∂ε∂t₀ = material.∂ε∂t₀
    Sᵥᵢₛ₀ > 0.0 && Δt > 0.0 && Δε > 0.0 ? Sᵥᵢₛ₀ * (Δε / Δt / ∂ε∂t₀)^(1.0 / m) : 0.0
end

function vol(A::Matrix{Float64})
    return tr(A) * I(3) / 3.0
end

function dev(A::Matrix{Float64})
    return A - vol(A)
end

function stress_update(
    material::J2,
    F::Matrix{Float64},
    Fᵖ::Matrix{Float64},
    εᵖ::Float64,
    Δt::Float64,
)
    max_rma_iter = 64
    max_ls_iter = 64

    κ = material.κ
    μ = material.μ
    λ = material.λ
    J = det(F)

    Fᵉ = F * inv(Fᵖ)
    Cᵉ = Fᵉ' * Fᵉ
    Eᵉ = 0.5 * log(Cᵉ)
    M = λ * tr(Eᵉ) * I(3) + 2.0 * μ * Eᵉ
    Mᵈᵉᵛ = dev(M)
    σᵛᵐ = sqrt(1.5) * norm(Mᵈᵉᵛ)
    σᵛᵒˡ = κ * vol(Eᵉ)

    Y = flow_strength(material, εᵖ)
    r = σᵛᵐ - Y
    r0 = r

    Δεᵖ = 0.0
    r_tol = 1e-10
    Δεᵖ_tol = 1e-10

    rma_iter = 0
    rma_converged = r ≤ r_tol
    while rma_converged == false
        if rma_iter == max_rma_iter
            break
        end
        Δεᵖ₀ = Δεᵖ
        merit_old = r * r
        H =
            hardening_rate(material, εᵖ + Δεᵖ) +
            viscoplastic_hardening_rate(material, Δεᵖ, Δt)
        ∂r = -3.0 * μ - H
        δεᵖ = -r / ∂r

        # line search
        ls_iter = 0
        α = 1.0
        backtrack_factor = 0.1
        decrease_factor = 1.0e-05
        ls_converged = false
        while ls_converged == false
            if ls_iter == max_ls_iter
                # line search has failed to satisfactorily improve newton step
                # just take the full newton step and hope for the best
                α = 1
                break
            end
            ls_iter += 1
            Δεᵖ = max(Δεᵖ₀ + α * δεᵖ, 0.0)
            Y = flow_strength(material, εᵖ + Δεᵖ) + viscoplastic_stress(material, Δεᵖ, Δt)
            r = σᵛᵐ - 3.0 * μ * Δεᵖ - Y

            merit_new = r * r
            decrease_tol = 1.0 - 2.0 * α * decrease_factor
            if merit_new <= decrease_tol * merit_old
                merit_old = merit_new
                ls_converged = true
            else
                α₀ = α
                α = α₀ * α₀ * merit_old / (merit_new - merit_old + 2.0 * α₀ * merit_old)
                if backtrack_factor * α₀ > α
                    α = backtrack_factor * α₀
                end
            end
        end
        rma_converged = abs(r / r0) < r_tol || Δεᵖ < Δεᵖ_tol
        rma_iter += 1
    end
    if rma_converged == false
        println("J2 stress update did not converge to specified tolerance")
    end

    Nᵖ = σᵛᵐ > 0.0 ? 1.5 * Mᵈᵉᵛ / σᵛᵐ : zeros(3, 3)
    ΔFᵖ = exp(Δεᵖ * Nᵖ)
    Fᵖ = ΔFᵖ * Fᵖ
    εᵖ += Δεᵖ

    ΔEᵉ = Δεᵖ * Nᵖ
    Mᵈᵉᵛ -= 2.0 * μ * ΔEᵉ
    σᵛᵐ = sqrt(1.5) * norm(Mᵈᵉᵛ)
    σᵈᵉᵛ = inv(Fᵉ)' * Mᵈᵉᵛ * Fᵉ' / J
    σ = σᵈᵉᵛ + σᵛᵒˡ

    eʸ = (σᵛᵐ - Y) / Y
    Fᵉ = F * inv(Fᵖ)
    Cᵉ = Fᵉ' * Fᵉ
    Eᵉ = 0.5 * log(Cᵉ)
    M = λ * tr(Eᵉ) * I(3) + 2.0 * μ * Eᵉ
    eᴹ = norm(Mᵈᵉᵛ - dev(M)) / norm(Mᵈᵉᵛ)
    return Fᵉ, Fᵖ, εᵖ, σ
end

struct Linear_Isotropic <: Thermal
    κ::Float64
    function Linear_Isotropic(params::Dict{Any,Any})
        κ = params["thermal conductivity"]
        new(κ)
    end
end

function odot(A::Matrix{Float64}, B::Matrix{Float64})
    n, _ = size(A)
    C = zeros(n, n, n, n)
    for a ∈ 1:n
        for b ∈ 1:n
            for c ∈ 1:n
                for d ∈ 1:n
                    C[a, b, c, d] = A[a, c] * B[b, d] + A[a, d] * B[b, c]
                end
            end
        end
    end
    C = 0.5 * C
    return C
end

function ox(A::Matrix{Float64}, B::Matrix{Float64})
    n, _ = size(A)
    C = zeros(n, n, n, n)
    for a ∈ 1:n
        for b ∈ 1:n
            for c ∈ 1:n
                for d ∈ 1:n
                    C[a, b, c, d] = A[a, b] * B[c, d]
                end
            end
        end
    end
    return C
end

function oxI(A::Matrix{Float64})
    n, _ = size(A)
    C = zeros(n, n, n, n)
    for a ∈ 1:n
        for b ∈ 1:n
            for c ∈ 1:n
                for d ∈ 1:n
                    C[a, b, c, d] = A[a, b] * I[c, d]
                end
            end
        end
    end
    return C
end

function Iox(B::Matrix{Float64})
    n, _ = size(B)
    C = zeros(n, n, n, n)
    for a ∈ 1:n
        for b ∈ 1:n
            for c ∈ 1:n
                for d ∈ 1:n
                    C[a, b, c, d] = I[a, b] * B[c, d]
                end
            end
        end
    end
    return C
end

function convect_tangent(CC::Array{Float64}, S::Matrix{Float64}, F::Matrix{Float64})
    n, _ = size(F)
    AA = zeros(n, n, n, n)
    for i ∈ 1:n
        for j ∈ 1:n
            for k ∈ 1:n
                for l ∈ 1:n
                    s = 0.0
                    for p ∈ 1:n
                        for q ∈ 1:n
                            s = s + F[i, p] * CC[p, j, l, q] * F[k, q]
                        end
                    end
                    AA[i, j, k, l] = S[l, j] * I[i, k] + s
                end
            end
        end
    end
    return AA
end

function second_from_fourth(AA::Array{Float64})
    n, _, _, _ = size(AA)
    A = zeros(n * n, n * n)
    for i ∈ 1:n
        for j ∈ 1:n
            p = n * (i - 1) + j
            for k ∈ 1:n
                for l ∈ 1:n
                    q = n * (k - 1) + l
                    A[p, q] = AA[i, j, k, l]
                end
            end
        end
    end
    return A
end

function constitutive(material::SaintVenant_Kirchhoff, F::Matrix{Float64})
    C = F' * F
    E = 0.5 * (C - I)
    λ = material.λ
    μ = material.μ
    trE = tr(E)
    W = 0.5 * λ * trE * trE + μ * tr(E * E)
    S = λ * trE * I + 2.0 * μ * E
    CC = zeros(3, 3, 3, 3)
    for i = 1:3
        for j = 1:3
            δᵢⱼ = I[i, j]
            for k = 1:3
                δᵢₖ = I[i, k]
                δⱼₖ = I[j, k]
                for l = 1:3
                    δᵢₗ = I[i, l]
                    δⱼₗ = I[j, l]
                    δₖₗ = I[k, l]
                    CC[i, j, k, l] = λ * δᵢⱼ * δₖₗ + μ * (δᵢₖ * δⱼₗ + δᵢₗ * δⱼₖ)
                end
            end
        end
    end
    P = F * S
    AA = convect_tangent(CC, S, F)
    return W, P, AA
end

function constitutive(material::Linear_Elastic, F::Matrix{Float64})
    ∇u = F - I
    ϵ = MiniTensor.symm(∇u)
    λ = material.λ
    μ = material.μ
    trϵ = tr(ϵ)
    W = 0.5 * λ * trϵ * trϵ + μ * tr(ϵ * ϵ)
    σ = λ * trϵ * I + 2.0 * μ * ϵ
    CC = zeros(3, 3, 3, 3)
    for i = 1:3
        for j = 1:3
            δᵢⱼ = I[i, j]
            for k = 1:3
                δᵢₖ = I[i, k]
                δⱼₖ = I[j, k]
                for l = 1:3
                    δᵢₗ = I[i, l]
                    δⱼₗ = I[j, l]
                    δₖₗ = I[k, l]
                    CC[i, j, k, l] = λ * δᵢⱼ * δₖₗ + μ * (δᵢₖ * δⱼₗ + δᵢₗ * δⱼₖ)
                end
            end
        end
    end
    return W, σ, CC
end

function constitutive(material::Neohookean, F::Matrix{Float64})
    C = F' * F
    J2 = det(C)
    Jm23 = 1.0 / cbrt(J2)
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1)
    Wdev = 0.5 * μ * (Jm23 * trC - 3)
    W = Wvol + Wdev
    IC = inv(C)
    Svol = 0.5 * κ * (J2 - 1) * IC
    Sdev = μ * Jm23 * (I - IC * trC / 3)
    S = Svol + Sdev
    ICxIC = ox(IC, IC)
    ICoIC = odot(IC, IC)
    μJ2n = 2.0 * μ * Jm23 / 3
    CCvol = κ * (J2 * ICxIC - (J2 - 1) * ICoIC)
    CCdev = μJ2n * (trC * (ICxIC / 3 + ICoIC) - oxI(IC) - Iox(IC))
    CC = CCvol + CCdev
    P = F * S
    AA = convect_tangent(CC, S, F)
    return W, P, AA
end

function energy_neohookean(F::Matrix{Float64}, material::NeohookeanAD)
    C = F' * F
    J2 = det(C)
    Jm23 = 1.0 / cbrt(J2)
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1)
    Wdev = 0.5 * μ * (Jm23 * trC - 3)
    W = Wvol + Wdev
    return W
end

function constitutive(material::NeohookeanAD, F::Matrix{Float64})
    C = F' * F
    J2 = det(C)
    Jm23 = 1.0 / cbrt(J2)
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1)
    Wdev = 0.5 * μ * (Jm23 * trC - 3)
    energy(F) = Wvol + Wdev
    W, P, AA = constitutive(material, energy, F)
    return W, P, AA
end

function constitutive(material::Solid, energy::Function, F::Matrix{Float64})
    C = MiniTensor.dot(MiniTensor.transpose(F), F)
    J2 = MiniTensor.determinant(C)
    Jm23 = 1.0 / cbrt(J2)
    trC = MiniTensor.trace(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1)
    Wdev = 0.5 * μ * (Jm23 * trC - 3)
    f(F) = Wvol + Wdev
    #f(F) = MiniTensor.determinant(F)
    W = energy(F)
    #P = reshape(collect(gradient(Forward, f, F)), 3, 3)
    println("*** DEBUG DEFGRAD : ", F)
    println("*** DEBUG ENERGY : ", typeof(f), f(F))
    println("*** DEBUG STRESS : ", gradient(Forward, f, F))
    stop
    AA = zeros(3,3,3,3)
    #AA = reshape(collect(gradient(Forward, (FF) -> reshape(collect(gradient(Forward, energy, FF)), 3, 3), F)), 3, 3, 3, 3)
    return W, P, AA
end

function create_material(params::Dict{Any,Any})
    model_name = params["model"]
    if model_name == "linear elastic"
        return Linear_Elastic(params)
    elseif model_name == "Saint-Venant Kirchhoff"
        return SaintVenant_Kirchhoff(params)
    elseif model_name == "neohookean"
        return Neohookean(params)
    elseif model_name == "neohookeanAD"
        return NeohookeanAD(params)
    else
        error("Unknown material model : ", model_name)
    end
end

function get_p_wave_modulus(material::Solid)
    return material.λ + 2.0 * material.μ
end
