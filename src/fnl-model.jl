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
include("I_a.jl")

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

#not integral computation implementaion
#window+GIC+RIC
@model function P_qso_convolved_IC(data, kₚ, Pₘ, p, α_k, f, Σ, Qₗ, W₀, start_k, Wric, fast::Bool)
    #prior
    f_nl ~ Uniform(-500, 500)
    b₁ ~ Uniform(0.1, 6)
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)

    #likelihood

    b_totk = p === nothing ? b_tot(b₁, f_nl, α_k) : b_tot(b₁, f_nl, p, α_k)
    I0 = compute_I_a(0, σ_fog, kₚ)
    I2 = compute_I_a(2, σ_fog, kₚ)
    I4 = compute_I_a(4, σ_fog, kₚ)
    I6 = compute_I_a(6, σ_fog, kₚ)
    I8 = compute_I_a(8, σ_fog, kₚ)

    #Pₘ is computed over the kₚs
    P0 = 0.5 .* Pₘ .* (b_totk.^2 .* I0 .+ 
                       2 .* b_totk .* f .* I2 .+ 
                       f.^2 .* I4) .+ N
    P2 = 5 .* 0.5 .* Pₘ .* (0.5 .* b_totk.^2 .* (3 .* I2 .- I0) .+ 
                            b_totk .* f .* (3 .* I4 .- I2) .+ 
                            0.5 .* f.^2 .* (3 .* I6 .- I4))
    P4 = 9 .* 0.5 .* Pₘ .* (0.125 .* b_totk.^2 .* (35 .* I4 .- 30 .* I2 .+ 3 .* I0) .+
                            0.25 .* b_totk .* f .* (35 .* I6 .- 30 .* I4 .+ 3 .* I2) .+
                            0.125 .* f.^2 .* (35 .* I8 .- 30 .* I6 .+ 3 .* I4))

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

@model function P_qso_convolved_IC_joint(dataN, dataS, kₚN, kₚS, PₘN, PₘS, p, α_kN, α_kS, fN, fS, ΣN, ΣS, QₗN, QₗS, W₀N, W₀S, start_kN, start_kS, WricN, WricS, fast::Bool)
    #prior
    f_nl ~ Uniform(-500, 500)
    b₁N ~ Uniform(0.1, 6)
    σ_fogN ~ Uniform(0, 20.)
    NN ~ Uniform(-5e3, 5e3)
    b₁S ~ Uniform(0.1, 6)
    σ_fogS ~ Uniform(0, 20.)
    NS ~ Uniform(-5e3, 5e3)

    #likelihood

    b_totkN = p === nothing ? b_tot(b₁N, f_nl, α_kN) : b_tot(b₁N, f_nl, p, α_kN)
    I0N = compute_I_a(0, σ_fogN, kₚN)
    I2N = compute_I_a(2, σ_fogN, kₚN)
    I4N = compute_I_a(4, σ_fogN, kₚN)
    I6N = compute_I_a(6, σ_fogN, kₚN)
    I8N = compute_I_a(8, σ_fogN, kₚN)

    b_totkS = p === nothing ? b_tot(b₁S, f_nl, α_kS) : b_tot(b₁S, f_nl, p, α_kS)
    I0S = compute_I_a(0, σ_fogS, kₚS)
    I2S = compute_I_a(2, σ_fogS, kₚS)
    I4S = compute_I_a(4, σ_fogS, kₚS)
    I6S = compute_I_a(6, σ_fogS, kₚS)
    I8S = compute_I_a(8, σ_fogS, kₚS)

    #Pₘ is computed over the kₚs
    P0N = 0.5 .* PₘN .* (b_totkN.^2 .* I0N .+ 
                         2 .* b_totkN .* fN .* I2N .+ 
                         fN.^2 .* I4N) .+ NN
    P2N = 5 .* 0.5 .* PₘN .* (0.5 .* b_totkN.^2 .* (3 .* I2N .- I0N) .+ 
                              b_totkN .* fN .* (3 .* I4N .- I2N) .+ 
                              0.5 .* fN.^2 .* (3 .* I6N .- I4N))
    P4N = 9 .* 0.5 .* PₘN .* (0.125 .* b_totkN.^2 .* (35 .* I4N .- 30 .* I2N .+ 3 .* I0N) .+
                              0.25 .* b_totkN .* fN .* (35 .* I6N .- 30 .* I4N .+ 3 .* I2N) .+
                              0.125 .* fN.^2 .* (35 .* I8N .- 30 .* I6N .+ 3 .* I4N))

    P0S = 0.5 .* PₘS .* (b_totkS.^2 .* I0S .+ 
                         2 .* b_totkS .* fS .* I2S .+ 
                         fS.^2 .* I4S) .+ NS
    P2S = 5 .* 0.5 .* PₘS .* (0.5 .* b_totkS.^2 .* (3 .* I2S .- I0S) .+ 
                              b_totkS .* fS .* (3 .* I4S .- I2S) .+ 
                              0.5 .* fS.^2 .* (3 .* I6S .- I4S))
    P4S = 9 .* 0.5 .* PₘS .* (0.125 .* b_totkS.^2 .* (35 .* I4S .- 30 .* I2S .+ 3 .* I0S) .+
                              0.25 .* b_totkS .* fS .* (35 .* I6S .- 30 .* I4S .+ 3 .* I2S) .+
                              0.125 .* fS.^2 .* (35 .* I8S .- 30 .* I6S .+ 3 .* I4S))

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