using Distributions
using Polynomials
using SpecialPolynomials
using LegendrePolynomials
using QuadGK: quadgk
using LinearAlgebra
using Tullio
using FastGaussQuadrature

using Turing

include("fnl-utils.jl")

function G(k, μ, σ_fog)
    return 1 ./ (1 .+ (k .* μ .* σ_fog).^2 ./ 2)
end

function b_tot(b₁, f_nl, p, α_k)
    f_term = f_nl .* (b₁ .- p) .* α_k
    return b₁ .+ f_term
end

function b_tot(b₁, f_nl, α_k)
    f_term = f_nl .* α_k
    return b₁ .+ f_term
end

function P_kmu(θ, k, μ, Pₘ, p, α_k, f)
    #Inference parameters
    f_nl = θ[1]
    b₁ = θ[2]
    σ_fog = θ[3]
    N = θ[4]

    #Theory
    G_k = G(k, μ, σ_fog)
    b_totk = p === nothing ? b_tot(b₁, f_nl, α_k) : b_tot(b₁, f_nl, p, α_k)

    P = G_k.^2 .* (b_totk .+ f .* μ^2).^2 .* Pₘ .+ N
    return P
end

function P_kmu_nofnl2(θ, k, μ, Pₘ, p, α_k, f)
    #Inference parameters
    f_nl = θ[1]
    b₁ = θ[2]
    σ_fog = θ[3]
    N = θ[4]

    #Theory
    G_k = G(k, μ, σ_fog)
    f_term = f_nl .* (b₁ .- p) .* α_k

    P = G_k.^2 .* (b₁.^2 .+ (f .* μ.^2).^2 .+ 2 .* b₁ .* (f .* μ.^2) .+ 2 .* f_term .* (b₁ .+ (f .* μ.^2))) .* Pₘ .+ N
end

function P₀(θ::Vector{T}, k, Pₘ, p, α_k, f) where T
    L₀ = Legendre([1])

    Pk = zeros(T, size(k))
    
    for i in eachindex(k)
        P_integrand = μ -> P_kmu(θ, k[i], μ, Pₘ[i], p, α_k[i], f)
        
        integrand(μ) = P_integrand(μ) * L₀(μ)
        P, err = quadgk(integrand, -1, 1, maxevals=1e7)

        Pk[i] = 1. / 2. * P
    
    end

    return Pk
end

function Pₗ_check(θ::Vector{T}, k, Pₘ, p, α_k, f, l) where T
    coeff = zeros(l+1)
    coeff[end] = 1
    Lₗ = Legendre(coeff)

    Pk = zeros(T, size(k))

    for i in eachindex(k)
        P_integrand = μ -> P_kmu(θ, k[i], μ, Pₘ[i], p, α_k[i], f)

        integrand(μ) = P_integrand(μ) * Lₗ(μ)
        P, err = quadgk(integrand, -1, 1, maxevals=1e7)

        Pk[i] = (2 * l + 1) / 2 * P
    
    end

    return Pk
end

function Pₗ(θ::Vector{T}, k, Pₘ, p, α_k, f, l) where T

    Pk = zeros(T, size(k))

    n_glb = 5
    nodes, weights = gausslobatto(n_glb*2)

    μ_nodes = nodes[1:n_glb]
    μ_weights = weights[1:n_glb]

    Pl_l = Pl.(μ_nodes, l)

    for i in eachindex(k)
        temp = zeros(T, size(μ_nodes))

        for j in eachindex(μ_nodes)
            temp[j] = P_kmu(θ, k[i], μ_nodes[j], Pₘ[i], p, α_k[i], f)
        end

        Pk[i] = (2 * l + 1) * _mygemmavx(μ_weights, temp, Pl_l)    
    end

    return Pk
end

@model function P_qso_unconvolved(data, k, Pₘ, p, α_k, f, Σ)
    #prior
    f_nl ~ Uniform(-500, 500)
    b₁ ~ Uniform(0.1, 6)
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)

    θ = [f_nl, b₁, σ_fog, N]
    prediction = P₀(θ, k, Pₘ, p, α_k, f)

    data ~ MvNormal(prediction, Σ)

    return nothing
end

@model function P_qso_convolved(data, kₚ, Pₘ, p, α_k, f, Σ, Qₗ)
    #prior
    f_nl ~ Uniform(-500, 500)
    b₁ ~ Uniform(0.1, 6)
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)

    θ = [f_nl, b₁, σ_fog, N]
    #likelihood
    #Pₘ is computed over the kₚs
    P0 = Pₗ(θ, kₚ, Pₘ, p, α_k, f, 0)
    P2 = Pₗ(θ, kₚ, Pₘ, p, α_k, f, 2)
    P4 = Pₗ(θ, kₚ, Pₘ, p, α_k, f, 4)

    #QN is a matrix k_eff x kₚ

    Q0 = Qₗ[:,:,1]
    Q2 = Qₗ[:,:,2]
    Q4 = Qₗ[:,:,3]

    @tullio prediction[i] := Q0[i,k] * P0[i,k] + Q2[i,k] * P2[i,k] + Q4[i,k] * P4[i,k]

    data ~ MvNormal(prediction, Σ)

    return nothing

end

#Global IC models

@model function P_qso_convolved_IC(data, kₚ, Pₘ, p, α_k, f, Σ, Qₗ, W₀, start_k)
    #prior
    f_nl ~ Uniform(-500, 500)
    b₁ ~ Uniform(0.1, 6)
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)

    θ = [f_nl, b₁, σ_fog, N]
    #likelihood
    #Pₘ is computed over the kₚs
    P0 = Pₗ(θ, kₚ, Pₘ, p, α_k, f, 0)
    P2 = Pₗ(θ, kₚ, Pₘ, p, α_k, f, 2)
    P4 = Pₗ(θ, kₚ, Pₘ, p, α_k, f, 4)

    #QN is a matrix k_eff x kₚ

    Q0 = Qₗ[:,:,1]
    Q2 = Qₗ[:,:,2]
    Q4 = Qₗ[:,:,3]

    @tullio convolvedPk[i] := Q0[i,k] * P0[i,k] + Q2[i,k] * P2[i,k] + Q4[i,k] * P4[i,k]

    #IC
    Pof0 = convolvedPk[1]
    convolvedPk = front_cut_ks(start_k, convolvedPk)
    prediction = convolvedPk .- Pof0 .* W₀

    data ~ MvNormal(prediction, Σ)

    return nothing

end

@model function P_qso_convolved_IC_joint(dataN, dataS, kₚN, kₚS, PₘN, PₘS, p, α_kN, α_kS, fN, fS, ΣN, ΣS, QₗN, QₗS, W₀N, W₀S, start_kN, start_kS)
    #prior
    f_nl ~ Uniform(-500, 500)
    b₁N ~ Uniform(0.1, 6)
    σ_fogN ~ Uniform(0, 20.)
    NN ~ Uniform(-5e3, 5e3)
    b₁S ~ Uniform(0.1, 6)
    σ_fogS ~ Uniform(0, 20.)
    NS ~ Uniform(-5e3, 5e3)

    θN = [f_nl, b₁N, σ_fogN, NN]
    θS = [f_nl, b₁S, σ_fogS, NS]
    #likelihood
    #Pₘ is computed over the kₚs
    P0N = Pₗ(θN, kₚN, PₘN, p, α_kN, fN, 0)
    P2N = Pₗ(θN, kₚN, PₘN, p, α_kN, fN, 2)
    P4N = Pₗ(θN, kₚN, PₘN, p, α_kN, fN, 4)

    P0S = Pₗ(θS, kₚS, PₘS, p, α_kS, fS, 0)
    P2S = Pₗ(θS, kₚS, PₘS, p, α_kS, fS, 2)
    P4S = Pₗ(θS, kₚS, PₘS, p, α_kS, fS, 4)

    #QN is a matrix k_eff x kₚ

    Q0N = QₗN[:,:,1]
    Q2N = QₗN[:,:,2]
    Q4N = QₗN[:,:,3]

    Q0S = QₗS[:,:,1]
    Q2S = QₗS[:,:,2]
    Q4S = QₗS[:,:,3]

    @tullio convolvedPkN[i] := Q0N[i,k] * P0N[i,k] + Q2N[i,k] * P2N[i,k] + Q4N[i,k] * P4N[i,k]
    @tullio convolvedPkS[i] := Q0S[i,k] * P0S[i,k] + Q2S[i,k] * P2S[i,k] + Q4S[i,k] * P4S[i,k]

    #IC
    Pof0N = convolvedPkN[1]
    convolvedPkN = front_cut_ks(start_kN, convolvedPkN)
    predictionN = convolvedPkN .- Pof0N .* W₀N

    Pof0S = convolvedPkS[1]
    convolvedPkS = front_cut_ks(start_kS, convolvedPkS)
    predictionS = convolvedPkS .- Pof0S .* W₀S

    dataN ~ MvNormal(predictionN, ΣN)
    dataS ~ MvNormal(predictionS, ΣS)

    return nothing

end

#Radial IC models
@model function P_qso_convolved_IC(data, kₚ, Pₘ, p, α_k, f, Σ, Qₗ, W₀, start_k, Wric)
    #prior
    f_nl ~ Uniform(-500, 500)
    b₁ ~ Uniform(0.1, 6)
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)

    θ = [f_nl, b₁, σ_fog, N]
    #likelihood
    #Pₘ is computed over the kₚs
    P0 = Pₗ(θ, kₚ, Pₘ, p, α_k, f, 0)
    P2 = Pₗ(θ, kₚ, Pₘ, p, α_k, f, 2)
    P4 = Pₗ(θ, kₚ, Pₘ, p, α_k, f, 4)

    #QN is a matrix k_eff x kₚ

    Q0 = Qₗ[:,:,1]
    Q2 = Qₗ[:,:,2]
    Q4 = Qₗ[:,:,3]

    @tullio convolvedPk[i] := Q0[i,k] * P0[i,k] + Q2[i,k] * P2[i,k] + Q4[i,k] * P4[i,k]

    #IC
    Pof0 = convolvedPk[1]
    convolvedPk = front_cut_ks(start_k, convolvedPk)
    prediction = convolvedPk .- Pof0 .* W₀ .- convolvedPk .* Wric

    data ~ MvNormal(prediction, Σ)

    return nothing

end

@model function P_qso_convolved_IC_joint(dataN, dataS, kₚN, kₚS, PₘN, PₘS, p, α_kN, α_kS, fN, fS, ΣN, ΣS, QₗN, QₗS, W₀N, W₀S, start_kN, start_kS, WricN, WricS)
    #prior
    f_nl ~ Uniform(-500, 500)
    b₁N ~ Uniform(0.1, 6)
    σ_fogN ~ Uniform(0, 20.)
    NN ~ Uniform(-5e3, 5e3)
    b₁S ~ Uniform(0.1, 6)
    σ_fogS ~ Uniform(0, 20.)
    NS ~ Uniform(-5e3, 5e3)

    θN = [f_nl, b₁N, σ_fogN, NN]
    θS = [f_nl, b₁S, σ_fogS, NS]
    #likelihood
    #Pₘ is computed over the kₚs
    P0N = Pₗ(θN, kₚN, PₘN, p, α_kN, fN, 0)
    P2N = Pₗ(θN, kₚN, PₘN, p, α_kN, fN, 2)
    P4N = Pₗ(θN, kₚN, PₘN, p, α_kN, fN, 4)

    P0S = Pₗ(θS, kₚS, PₘS, p, α_kS, fS, 0)
    P2S = Pₗ(θS, kₚS, PₘS, p, α_kS, fS, 2)
    P4S = Pₗ(θS, kₚS, PₘS, p, α_kS, fS, 4)

    #QN is a matrix k_eff x kₚ

    Q0N = QₗN[:,:,1]
    Q2N = QₗN[:,:,2]
    Q4N = QₗN[:,:,3]

    Q0S = QₗS[:,:,1]
    Q2S = QₗS[:,:,2]
    Q4S = QₗS[:,:,3]

    @tullio convolvedPkN[i] := Q0N[i,k] * P0N[i,k] + Q2N[i,k] * P2N[i,k] + Q4N[i,k] * P4N[i,k]
    @tullio convolvedPkS[i] := Q0S[i,k] * P0S[i,k] + Q2S[i,k] * P2S[i,k] + Q4S[i,k] * P4S[i,k]

    #IC
    Pof0N = convolvedPkN[1]
    convolvedPkN = front_cut_ks(start_kN, convolvedPkN)
    predictionN = convolvedPkN .- Pof0N .* W₀N .- convolvedPkN .* WricN

    Pof0S = convolvedPkS[1]
    convolvedPkS = front_cut_ks(start_kS, convolvedPkS)
    predictionS = convolvedPkS .- Pof0S .* W₀S .- convolvedPkS .* WricS

    dataN ~ MvNormal(predictionN, ΣN)
    dataS ~ MvNormal(predictionS, ΣS)

    return nothing

end