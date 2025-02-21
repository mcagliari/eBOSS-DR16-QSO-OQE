include("../src/fnl-model.jl")
include("fnl-bispectrum_model.jl")

function get_b₁_corr_covariance(σs, ϵ)
    Σ = I(length(σs)) .* σs
    corr = zeros((length(σs), length(σs)))

    for i in 1:length(σs)
        for j in 1:length(σs)
            if i != j
                corr[i,j] = 1 - 0.5 * ϵ^2
                corr[i,j] *= sqrt(Σ[i,i]) * sqrt(Σ[j,j])
            end
        end
    end

    Σ_corr = Σ .+ corr
    return Σ_corr
end

#single field
@model function PpB_qso_convolved_IC(dataPB, p, iΓPB, kP, kB, PPₘ, PBₘ, α_kP, MinvB, fP, fB, QₗP, W₀P, start_kP, WricP, WricB, IabcB, PsnB)
    #prior: common to P and B
    fNL ~ Uniform(-500, 500)
    #b₁P ~ Uniform(0.2, 6) 
    #b₁B ~ Uniform(0.1, b₁P)
    ϵ = 0.2
    μbs = [2.30, 2.28] # [b₁P, b₁B]
    σbs = [2., 2.]
    Σb₁ = get_b₁_corr_covariance(σbs, ϵ)
    b₁s ~ MvNormal(μbs, Σb₁)
    #prior: P only
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)
    #prior: B only
    b₂ ~ Uniform(-4, 4) 
    bₛ ~ Uniform(-4, 4)
    c₁ ~ Uniform(-1, 1)
    α₁ ~ Uniform(-1.5, 1.5)
    α₂ ~ Uniform(-1.5, 1.5)
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕ = 2 * δ_c * (b₁s[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁s[2] - 1)) - b₁s[2] + 1) #usiamo Universal relation

    #parameter vectors
    θ = [fNL, b₁s[1], σ_fog, N] #for P
    bias = [b₁s[2], b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂] #for B

    #P likelihood
    #PPₘ is computed over the kPs
    P0 = Pₗ(θ, kP, PPₘ, p, α_kP, fP, 0)
    P2 = Pₗ(θ, kP, PPₘ, p, α_kP, fP, 2)
    P4 = Pₗ(θ, kP, PPₘ, p, α_kP, fP, 4)

    #QN is a matrix k_eff x kₚ

    Q0 = QₗP[:,:,1]
    Q2 = QₗP[:,:,2]
    Q4 = QₗP[:,:,3]

    @tullio convolvedPk[i] := Q0[i,k] * P0[i,k] + Q2[i,k] * P2[i,k] + Q4[i,k] * P4[i,k]

    #IC
    Pof0 = convolvedPk[1]
    convolvedPk = front_cut_ks(start_kP, convolvedPk)
    predictionP = convolvedPk .- Pof0 .* W₀P .- convolvedPk .* WricP

    #B likelihood
    B0 = compute_B0(kB, PBₘ, MinvB, IabcB, bias, fB, fNL, PsnB)

    #RIC contribution (GIC contribution in P0_model)
    predictionB = B0 .- B0 .* WricB

    predictionPB = vcat(vec(predictionP), predictionB)
    predictionPB_recast = iΓPB * predictionPB

    dataPB ~ MvNormal(predictionPB_recast, I)

    return nothing

end

#does not recompute B0_terms
@model function PpB_qso_convolved_IC(dataPB, p, iΓPB, kP, PPₘ, α_kP, B0_terms, BGIC_terms, fP, fB, QₗP, W₀P, start_kP, WricP, WricB, PsnB)
    #prior: common to P and b
    fNL ~ Uniform(-500, 500)
    #b₁P ~ Uniform(0.2, 6) 
    #b₁B ~ Uniform(0.1, b₁P)
    ϵ = 0.2
    μbs = [2.30, 2.28] # [b₁P, b₁B]
    σbs = [2., 2.]
    Σb₁ = get_b₁_corr_covariance(σbs, ϵ)
    b₁s ~ MvNormal(μbs, Σb₁)
    #prior: P only
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)
    #prior: B only
    b₂ ~ Uniform(-4, 4)
    bₛ ~ Uniform(-4, 4)
    c₁ ~ Uniform(-1, 1) 
    α₁ ~ Uniform(-1.5, 1.5) 
    α₂ ~ Uniform(-1.5, 1.5)
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕ = 2 * δ_c * (b₁s[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁s[2] - 1)) - b₁s[2] + 1) #usiamo Universal relation

    #parameter vectors
    θ = [fNL, b₁s[1], σ_fog, N] #for P
    bias = [b₁s[2], b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂] #for B

    #P likelihood
    #PPₘ is computed over the kPs
    P0 = Pₗ(θ, kP, PPₘ, p, α_kP, fP, 0)
    P2 = Pₗ(θ, kP, PPₘ, p, α_kP, fP, 2)
    P4 = Pₗ(θ, kP, PPₘ, p, α_kP, fP, 4)

    #QN is a matrix k_eff x kₚ

    Q0 = QₗP[:,:,1]
    Q2 = QₗP[:,:,2]
    Q4 = QₗP[:,:,3]

    @tullio convolvedPk[i] := Q0[i,k] * P0[i,k] + Q2[i,k] * P2[i,k] + Q4[i,k] * P4[i,k]

    #IC
    Pof0 = convolvedPk[1]
    convolvedPk = front_cut_ks(start_kP, convolvedPk)
    predictionP = convolvedPk .- Pof0 .* W₀P .- convolvedPk .* WricP

    #B likelihood
    B0 = compute_B0(B0_terms, bias, fB, fNL, PsnB)
    BGIC = compute_B0(BGIC_terms, bias, fB, fNL, PsnB)

    #RIC contribution (GIC contribution in P0_model)
    predictionB = BGIC .- B0 .* WricB

    predictionPB = vcat(vec(predictionP), predictionB)
    predictionPB_recast = iΓPB * predictionPB

    dataPB ~ MvNormal(predictionPB_recast, I)

    return nothing

end

#No IC in B model
@model function PpB_qso_convolved(dataPB, p, iΓPB, kP, PPₘ, α_kP, B0_terms, fP, fB, QₗP, W₀P, start_kP, WricP, PsnB)
    #prior: common to P and b
    fNL ~ Uniform(-500, 500)
    #b₁P ~ Uniform(0.2, 6) 
    #b₁B ~ Uniform(0.1, b₁P)
    ϵ = 0.2
    μbs = [2.30, 2.28] # [b₁P, b₁B]
    σbs = [2., 2.]
    Σb₁ = get_b₁_corr_covariance(σbs, ϵ)
    b₁s ~ MvNormal(μbs, Σb₁)
    #prior: P only
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)
    #prior: B only
    b₂ ~ Uniform(-4, 4)
    bₛ ~ Uniform(-4, 4)
    c₁ ~ Uniform(-1, 1) 
    α₁ ~ Uniform(-1.5, 1.5) 
    α₂ ~ Uniform(-1.5, 1.5)
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕ = 2 * δ_c * (b₁s[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁s[2] - 1)) - b₁s[2] + 1) #usiamo Universal relation

    #parameter vectors
    θ = [fNL, b₁s[1], σ_fog, N] #for P
    bias = [b₁s[2], b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂] #for B

    #P likelihood
    #PPₘ is computed over the kPs
    P0 = Pₗ(θ, kP, PPₘ, p, α_kP, fP, 0)
    P2 = Pₗ(θ, kP, PPₘ, p, α_kP, fP, 2)
    P4 = Pₗ(θ, kP, PPₘ, p, α_kP, fP, 4)

    #QN is a matrix k_eff x kₚ

    Q0 = QₗP[:,:,1]
    Q2 = QₗP[:,:,2]
    Q4 = QₗP[:,:,3]

    @tullio convolvedPk[i] := Q0[i,k] * P0[i,k] + Q2[i,k] * P2[i,k] + Q4[i,k] * P4[i,k]

    #IC
    Pof0 = convolvedPk[1]
    convolvedPk = front_cut_ks(start_kP, convolvedPk)
    predictionP = convolvedPk .- Pof0 .* W₀P .- convolvedPk .* WricP

    #B likelihood
    predictionB = compute_B0(B0_terms, bias, fB, fNL, PsnB)

    predictionPB = vcat(vec(predictionP), predictionB)
    predictionPB_recast = iΓPB * predictionPB

    dataPB ~ MvNormal(predictionPB_recast, I)

    return nothing

end

#joint analysis NGC and SGC
@model function PpB_qso_convolved_IC(dataPBN, dataPBS, p, iΓPBN, iΓPBS, kPN, kPS, kBN, kBS, PPₘN, PPₘS, PBₘN, PBₘS, α_kPN, α_kPS, MinvBN, MinvBS, fPN, fPS, fBN, fBS, QₗPN, QₗPS, W₀PN, W₀PS, start_kPN, start_kPS, WricPN, WricPS, WricBN, WricBS, IabcBN, IabcBS, PsnBN, PsnBS)
    #prior: common to P and b
    fNL ~ Uniform(-500, 500)
    #b₁PN ~ Uniform(0.2, 6)
    #b₁PS ~ Uniform(0.2, 6) 
    #b₁BN ~ Uniform(0.1, b₁PN)
    #b₁BS ~ Uniform(0.1, b₁PS)
    ϵ = 0.2
    μbsN = [2.30, 2.28] # [b₁P, b₁B]
    σbsN = [2., 2.]
    Σb₁N = get_b₁_corr_covariance(σbsN, ϵ)
    b₁sN ~ MvNormal(μbsN, Σb₁N)
    μbsS = [2.30, 2.28] # [b₁P, b₁B]
    σbsS = [2., 2.]
    Σb₁S = get_b₁_corr_covariance(σbsS, ϵ)
    b₁sS ~ MvNormal(μbsS, Σb₁S)
    #prior: P only
    σ_fogN ~ Uniform(0, 20.)
    NN ~ Uniform(-5e3, 5e3)
    σ_fogS ~ Uniform(0, 20.)
    NS ~ Uniform(-5e3, 5e3)
    #prior: B only
    b₂N ~ Uniform(-4, 4)
    bₛN ~ Uniform(-4, 4) 
    c₁N ~ Uniform(-1, 1) 
    α₁N ~ Uniform(-1.5, 1.5)
    α₂N ~ Uniform(-1.5, 1.5)
    b₂S ~ Uniform(-4, 4)
    bₛS ~ Uniform(-4, 4)
    c₁S ~ Uniform(-1, 1)
    α₁S ~ Uniform(-1.5, 1.5)
    α₂S ~ Uniform(-1.5, 1.5)
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕN = 2 * δ_c * (b₁sN[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδN = bϕN + 2 * (δ_c * (b₂N - 8 / 21 * (b₁sN[2] - 1)) - b₁sN[2] + 1) #usiamo Universal relation
    bϕS = 2 * δ_c * (b₁sS[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδS = bϕS + 2 * (δ_c * (b₂S - 8 / 21 * (b₁sS[2] - 1)) - b₁sS[2] + 1) #usiamo Universal relation

    #parameter vectors
    θN = [fNL, b₁sN[1], σ_fogN, NN] #for P
    θS = [fNL, b₁sS[1], σ_fogS, NS] #for P
    biasN = [b₁sN[2], b₂N, bₛN, bϕN, bϕδN, c₁N, α₁N, α₂N] #for B
    biasS = [b₁sS[2], b₂S, bₛS, bϕS, bϕδS, c₁S, α₁S, α₂S] #for B

    #P likelihood
    #Pₘ is computed over the kₚs
    P0N = Pₗ(θN, kPN, PPₘN, p, α_kPN, fPN, 0)
    P2N = Pₗ(θN, kPN, PPₘN, p, α_kPN, fPN, 2)
    P4N = Pₗ(θN, kPN, PPₘN, p, α_kPN, fPN, 4)

    P0S = Pₗ(θS, kPS, PPₘS, p, α_kPS, fPS, 0)
    P2S = Pₗ(θS, kPS, PPₘS, p, α_kPS, fPS, 2)
    P4S = Pₗ(θS, kPS, PPₘS, p, α_kPS, fPS, 4)

    #QN is a matrix k_eff x kₚ

    Q0N = QₗPN[:,:,1]
    Q2N = QₗPN[:,:,2]
    Q4N = QₗPN[:,:,3]

    Q0S = QₗPS[:,:,1]
    Q2S = QₗPS[:,:,2]
    Q4S = QₗPS[:,:,3]

    @tullio convolvedPkN[i] := Q0N[i,k] * P0N[i,k] + Q2N[i,k] * P2N[i,k] + Q4N[i,k] * P4N[i,k]
    @tullio convolvedPkS[i] := Q0S[i,k] * P0S[i,k] + Q2S[i,k] * P2S[i,k] + Q4S[i,k] * P4S[i,k]

    #IC
    Pof0N = convolvedPkN[1]
    convolvedPkN = front_cut_ks(start_kPN, convolvedPkN)
    predictionPN = convolvedPkN .- Pof0N .* W₀PN .- convolvedPkN .* WricPN

    Pof0S = convolvedPkS[1]
    convolvedPkS = front_cut_ks(start_kPS, convolvedPkS)
    predictionPS = convolvedPkS .- Pof0S .* W₀PS .- convolvedPkS .* WricPS

    #B likelihood
    B0N = compute_B0(kBN, PBₘN, MinvBN, IabcBN, biasN, fBN, fNL, PsnBN)
    B0S = compute_B0(kBS, PBₘS, MinvBS, IabcBS, biasS, fBS, fNL, PsnBS)

    #RIC contribution (GIC contribution in P0_model)
    predictionBN = B0N .- B0N .* WricBN
    predictionBS = B0S .- B0S .* WricBS

    predictionPBN = vcat(vec(predictionPN), predictionBN)
    predictionPBS = vcat(vec(predictionPS), predictionBS)
    predictionPBN_recast = iΓPBN * predictionPBN
    predictionPBS_recast = iΓPBS * predictionPBS

    dataPBN ~ MvNormal(predictionPBN_recast, I)
    dataPBS ~ MvNormal(predictionPBS_recast, I)

    return nothing
    
end

#B terms not recomputed
@model function PpB_qso_convolved_IC(dataPBN, dataPBS, p, iΓPBN, iΓPBS, kPN, kPS, PPₘN, PPₘS, α_kPN, α_kPS, B0_termsN, B0_termsS, BGIC_termsN, BGIC_termsS, fPN, fPS, fBN, fBS, QₗPN, QₗPS, W₀PN, W₀PS, start_kPN, start_kPS, WricPN, WricPS, WricBN, WricBS, PsnBN, PsnBS)
    #prior: common to P and b
    fNL ~ Uniform(-500, 500)
    #b₁PN ~ Uniform(0.2, 6)
    #b₁PS ~ Uniform(0.2, 6)
    #b₁BN ~ Uniform(0.1, b₁PN)
    #b₁BS ~ Uniform(0.1, b₁PS)
    ϵ = 0.2
    μbsN = [2.30, 2.28] # [b₁P, b₁B]
    σbsN = [2., 2.]
    Σb₁N = get_b₁_corr_covariance(σbsN, ϵ)
    b₁sN ~ MvNormal(μbsN, Σb₁N)
    μbsS = [2.30, 2.28] # [b₁P, b₁B]
    σbsS = [2., 2.]
    Σb₁S = get_b₁_corr_covariance(σbsS, ϵ)
    b₁sS ~ MvNormal(μbsS, Σb₁S)
    #prior: P only
    σ_fogN ~ Uniform(0, 20.)
    NN ~ Uniform(-5e3, 5e3)
    σ_fogS ~ Uniform(0, 20.)
    NS ~ Uniform(-5e3, 5e3)
    #prior: B only
    b₂N ~ Uniform(-4, 4) 
    bₛN ~ Uniform(-4, 4)
    c₁N ~ Uniform(-1, 1) 
    α₁N ~ Uniform(-1.5, 1.5)
    α₂N ~ Uniform(-1.5, 1.5)
    b₂S ~ Uniform(-4, 4) 
    bₛS ~ Uniform(-4, 4)
    c₁S ~ Uniform(-1, 1) 
    α₁S ~ Uniform(-1.5, 1.5)
    α₂S ~ Uniform(-1.5, 1.5)
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕN = 2 * δ_c * (b₁sN[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδN = bϕN + 2 * (δ_c * (b₂N - 8 / 21 * (b₁sN[2] - 1)) - b₁sN[2] + 1) #usiamo Universal relation
    bϕS = 2 * δ_c * (b₁sS[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδS = bϕS + 2 * (δ_c * (b₂S - 8 / 21 * (b₁sS[2] - 1)) - b₁sS[2] + 1) #usiamo Universal relation

    #parameter vectors
    θN = [fNL, b₁sN[1], σ_fogN, NN] #for P
    θS = [fNL, b₁sS[1], σ_fogS, NS] #for P
    biasN = [b₁sN[2], b₂N, bₛN, bϕN, bϕδN, c₁N, α₁N, α₂N] #for B
    biasS = [b₁sS[2], b₂S, bₛS, bϕS, bϕδS, c₁S, α₁S, α₂S] #for B

    #P likelihood
    #Pₘ is computed over the kₚs
    P0N = Pₗ(θN, kPN, PPₘN, p, α_kPN, fPN, 0)
    P2N = Pₗ(θN, kPN, PPₘN, p, α_kPN, fPN, 2)
    P4N = Pₗ(θN, kPN, PPₘN, p, α_kPN, fPN, 4)

    P0S = Pₗ(θS, kPS, PPₘS, p, α_kPS, fPS, 0)
    P2S = Pₗ(θS, kPS, PPₘS, p, α_kPS, fPS, 2)
    P4S = Pₗ(θS, kPS, PPₘS, p, α_kPS, fPS, 4)

    #QN is a matrix k_eff x kₚ

    Q0N = QₗPN[:,:,1]
    Q2N = QₗPN[:,:,2]
    Q4N = QₗPN[:,:,3]

    Q0S = QₗPS[:,:,1]
    Q2S = QₗPS[:,:,2]
    Q4S = QₗPS[:,:,3]

    @tullio convolvedPkN[i] := Q0N[i,k] * P0N[i,k] + Q2N[i,k] * P2N[i,k] + Q4N[i,k] * P4N[i,k]
    @tullio convolvedPkS[i] := Q0S[i,k] * P0S[i,k] + Q2S[i,k] * P2S[i,k] + Q4S[i,k] * P4S[i,k]

    #IC
    Pof0N = convolvedPkN[1]
    convolvedPkN = front_cut_ks(start_kPN, convolvedPkN)
    predictionPN = convolvedPkN .- Pof0N .* W₀PN .- convolvedPkN .* WricPN

    Pof0S = convolvedPkS[1]
    convolvedPkS = front_cut_ks(start_kPS, convolvedPkS)
    predictionPS = convolvedPkS .- Pof0S .* W₀PS .- convolvedPkS .* WricPS

    #B likelihood
    B0N = compute_B0(B0_termsN, biasN, fBN, fNL, PsnBN)
    B0S = compute_B0(B0_termsS, biasS, fBS, fNL, PsnBS)
    BGICN = compute_B0(BGIC_termsN, biasN, fBN, fNL, PsnBN)
    BGICS = compute_B0(BGIC_termsS, biasS, fBS, fNL, PsnBS)

    #RIC contribution (GIC contribution in P0_model)
    predictionBN = BGICN .- B0N .* WricBN
    predictionBS = BGICS .- B0S .* WricBS

    predictionPBN = vcat(vec(predictionPN), predictionBN)
    predictionPBS = vcat(vec(predictionPS), predictionBS)
    predictionPBN_recast = iΓPBN * predictionPBN
    predictionPBS_recast = iΓPBS * predictionPBS

    dataPBN ~ MvNormal(predictionPBN_recast, I)
    dataPBS ~ MvNormal(predictionPBS_recast, I)

    return nothing
    
end

#Fast P implementation
#single field
@model function PpB_qso_convolved_IC(dataPB, p, iΓPB, kP, PPₘ, α_kP, B0_terms, BGIC_terms, fP, fB, QₗP, W₀P, start_kP, WricP, WricB, PsnB, fast::Bool)
    #prior: common to P and b
    fNL ~ Uniform(-500, 500)
    #b₁P ~ Uniform(0.2, 6) 
    #b₁B ~ Uniform(0.1, b₁P)
    ϵ = 0.2
    μbs = [2.30, 2.28] # [b₁P, b₁B]
    σbs = [2., 2.]
    Σb₁ = get_b₁_corr_covariance(σbs, ϵ)
    b₁s ~ MvNormal(μbs, Σb₁)
    #prior: P only
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)
    #prior: B only
    b₂ ~ Uniform(-4, 4) #Uniform(-20, 20) #large prior test
    bₛ ~ Uniform(-4, 4)
    c₁ ~ Uniform(-1, 1) #Uniform(-100, 100) #large prior test 
    α₁ ~ Uniform(-1.5, 1.5) 
    α₂ ~ Uniform(-1.5, 1.5)
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕ = 2 * δ_c * (b₁s[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁s[2] - 1)) - b₁s[2] + 1) #usiamo Universal relation

    #parameter vectors
    θ = [fNL, b₁s[1], σ_fog, N] #for P
    bias = [b₁s[2], b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂] #for B

    #P likelihood
   
    b_totk = p === nothing ? b_tot(b₁s[1], fNL, α_kP) : b_tot(b₁s[1], fNL, p, α_kP)
    I0 = compute_I_a(0, σ_fog, kP)
    I2 = compute_I_a(2, σ_fog, kP)
    I4 = compute_I_a(4, σ_fog, kP)
    I6 = compute_I_a(6, σ_fog, kP)
    I8 = compute_I_a(8, σ_fog, kP)

    #PPₘ is computed over the kₚs
    P0 = 0.5 .* PPₘ .* (b_totk.^2 .* I0 .+ 
                       2 .* b_totk .* fP .* I2 .+ 
                       fP.^2 .* I4) .+ N
    P2 = 5 .* 0.5 .* PPₘ .* (0.5 .* b_totk.^2 .* (3 .* I2 .- I0) .+ 
                            b_totk .* fP .* (3 .* I4 .- I2) .+ 
                            0.5 .* fP.^2 .* (3 .* I6 .- I4))
    P4 = 9 .* 0.5 .* PPₘ .* (0.125 .* b_totk.^2 .* (35 .* I4 .- 30 .* I2 .+ 3 .* I0) .+
                            0.25 .* b_totk .* fP .* (35 .* I6 .- 30 .* I4 .+ 3 .* I2) .+
                            0.125 .* fP.^2 .* (35 .* I8 .- 30 .* I6 .+ 3 .* I4))

    #QN is a matrix k_eff x kₚ

    Q0 = QₗP[:,:,1]
    Q2 = QₗP[:,:,2]
    Q4 = QₗP[:,:,3]

    @tullio convolvedPk[i] := Q0[i,k] * P0[i,k] + Q2[i,k] * P2[i,k] + Q4[i,k] * P4[i,k]

    #IC
    Pof0 = convolvedPk[1]
    convolvedPk = front_cut_ks(start_kP, convolvedPk)
    predictionP = convolvedPk .- Pof0 .* W₀P .- convolvedPk .* WricP

    #B likelihood
    B0 = compute_B0(B0_terms, bias, fB, fNL, PsnB)
    BGIC = compute_B0(BGIC_terms, bias, fB, fNL, PsnB)

    #RIC contribution (GIC contribution in P0_model)
    predictionB = BGIC .- B0 .* WricB

    predictionPB = vcat(vec(predictionP), predictionB)
    predictionPB_recast = iΓPB * predictionPB

    dataPB ~ MvNormal(predictionPB_recast, I)

    return nothing

end

#No IC in B model
@model function PpB_qso_convolved(dataPB, p, iΓPB, kP, PPₘ, α_kP, B0_terms, fP, fB, QₗP, W₀P, start_kP, WricP, PsnB, fast::Bool)
    #prior: common to P and b
    fNL ~ Uniform(-500, 500)
    #b₁P ~ Uniform(0.2, 6) 
    #b₁B ~ Uniform(0.1, b₁P)
    ϵ = 0.2
    μbs = [2.30, 2.28] # [b₁P, b₁B]
    σbs = [2., 2.]
    Σb₁ = get_b₁_corr_covariance(σbs, ϵ)
    b₁s ~ MvNormal(μbs, Σb₁)
    #prior: P only
    σ_fog ~ Uniform(0, 20.)
    N ~ Uniform(-5e3, 5e3)
    #prior: B only
    b₂ ~ Uniform(-4, 4)
    bₛ ~ Uniform(-4, 4)
    c₁ ~ Uniform(-1, 1) 
    α₁ ~ Uniform(-1.5, 1.5) 
    α₂ ~ Uniform(-1.5, 1.5)
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕ = 2 * δ_c * (b₁s[2] - p) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁s[2] - 1)) - b₁s[2] + 1) #usiamo Universal relation

    #parameter vectors
    θ = [fNL, b₁s[1], σ_fog, N] #for P
    bias = [b₁s[2], b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂] #for B

    #P likelihood
    b_totk = p === nothing ? b_tot(b₁s[1], fNL, α_kP) : b_tot(b₁s[1], fNL, p, α_kP)
    I0 = compute_I_a(0, σ_fog, kP)
    I2 = compute_I_a(2, σ_fog, kP)
    I4 = compute_I_a(4, σ_fog, kP)
    I6 = compute_I_a(6, σ_fog, kP)
    I8 = compute_I_a(8, σ_fog, kP)

    #PPₘ is computed over the kₚs
    P0 = 0.5 .* PPₘ .* (b_totk.^2 .* I0 .+ 
                       2 .* b_totk .* fP .* I2 .+ 
                       fP.^2 .* I4) .+ N
    P2 = 5 .* 0.5 .* PPₘ .* (0.5 .* b_totk.^2 .* (3 .* I2 .- I0) .+ 
                            b_totk .* fP .* (3 .* I4 .- I2) .+ 
                            0.5 .* fP.^2 .* (3 .* I6 .- I4))
    P4 = 9 .* 0.5 .* PPₘ .* (0.125 .* b_totk.^2 .* (35 .* I4 .- 30 .* I2 .+ 3 .* I0) .+
                            0.25 .* b_totk .* fP .* (35 .* I6 .- 30 .* I4 .+ 3 .* I2) .+
                            0.125 .* fP.^2 .* (35 .* I8 .- 30 .* I6 .+ 3 .* I4))

    #QN is a matrix k_eff x kₚ

    Q0 = QₗP[:,:,1]
    Q2 = QₗP[:,:,2]
    Q4 = QₗP[:,:,3]

    @tullio convolvedPk[i] := Q0[i,k] * P0[i,k] + Q2[i,k] * P2[i,k] + Q4[i,k] * P4[i,k]

    #IC
    Pof0 = convolvedPk[1]
    convolvedPk = front_cut_ks(start_kP, convolvedPk)
    predictionP = convolvedPk .- Pof0 .* W₀P .- convolvedPk .* WricP

    #B likelihood
    predictionB = compute_B0(B0_terms, bias, fB, fNL, PsnB)

    predictionPB = vcat(vec(predictionP), predictionB)
    predictionPB_recast = iΓPB * predictionPB

    dataPB ~ MvNormal(predictionPB_recast, I)

    return nothing

end

#joint
@model function PpB_qso_convolved_IC(dataPBN, dataPBS, p, iΓPBN, iΓPBS, kPN, kPS, PPₘN, PPₘS, α_kPN, α_kPS, B0_termsN, B0_termsS, BGIC_termsN, BGIC_termsS, fPN, fPS, fBN, fBS, QₗPN, QₗPS, W₀PN, W₀PS, start_kPN, start_kPS, WricPN, WricPS, WricBN, WricBS, PsnBN, PsnBS, fast::Bool)
    #prior: common to P and b
    fNL ~ Uniform(-500, 500)
    #b₁PN ~ Uniform(0.2, 6)
    #b₁PS ~ Uniform(0.2, 6)
    #b₁BN ~ Uniform(0.1, b₁PN)
    #b₁BS ~ Uniform(0.1, b₁PS)
    ϵ = 0.2
    μbsN = [2.30, 2.28] # [b₁P, b₁B]
    σbsN = [2., 2.]
    Σb₁N = get_b₁_corr_covariance(σbsN, ϵ)
    b₁sN ~ MvNormal(μbsN, Σb₁N)
    μbsS = [2.30, 2.28] # [b₁P, b₁B]
    σbsS = [2., 2.]
    Σb₁S = get_b₁_corr_covariance(σbsS, ϵ)
    b₁sS ~ MvNormal(μbsS, Σb₁S)
    #prior: P only
    σ_fogN ~ Uniform(0, 20.)
    NN ~ Uniform(-5e3, 5e3)
    σ_fogS ~ Uniform(0, 20.)
    NS ~ Uniform(-5e3, 5e3)
    #prior: B only
    b₂N ~ Uniform(-4, 4) 
    bₛN ~ Uniform(-4, 4)
    c₁N ~ Uniform(-1, 1) 
    α₁N ~ Uniform(-1.5, 1.5)
    α₂N ~ Uniform(-1.5, 1.5)
    b₂S ~ Uniform(-4, 4) 
    bₛS ~ Uniform(-4, 4)
    c₁S ~ Uniform(-1, 1) 
    α₁S ~ Uniform(-1.5, 1.5)
    α₂S ~ Uniform(-1.5, 1.5)
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕN = 2 * δ_c * (b₁sN[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδN = bϕN + 2 * (δ_c * (b₂N - 8 / 21 * (b₁sN[2] - 1)) - b₁sN[2] + 1) #usiamo Universal relation
    bϕS = 2 * δ_c * (b₁sS[2] - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδS = bϕS + 2 * (δ_c * (b₂S - 8 / 21 * (b₁sS[2] - 1)) - b₁sS[2] + 1) #usiamo Universal relation

    #parameter vectors
    θN = [fNL, b₁sN[1], σ_fogN, NN] #for P
    θS = [fNL, b₁sS[1], σ_fogS, NS] #for P
    biasN = [b₁sN[2], b₂N, bₛN, bϕN, bϕδN, c₁N, α₁N, α₂N] #for B
    biasS = [b₁sS[2], b₂S, bₛS, bϕS, bϕδS, c₁S, α₁S, α₂S] #for B

    #P likelihood
    
    b_totkN = p === nothing ? b_tot(b₁sN[1], fNL, α_kPN) : b_tot(b₁sN[1], fNL, p, α_kPN)
    I0N = compute_I_a(0, σ_fogN, kPN)
    I2N = compute_I_a(2, σ_fogN, kPN)
    I4N = compute_I_a(4, σ_fogN, kPN)
    I6N = compute_I_a(6, σ_fogN, kPN)
    I8N = compute_I_a(8, σ_fogN, kPN)

    b_totkS = p === nothing ? b_tot(b₁sS[1], fNL, α_kPS) : b_tot(b₁sS[1], fNL, p, α_kPS)
    I0S = compute_I_a(0, σ_fogS, kPS)
    I2S = compute_I_a(2, σ_fogS, kPS)
    I4S = compute_I_a(4, σ_fogS, kPS)
    I6S = compute_I_a(6, σ_fogS, kPS)
    I8S = compute_I_a(8, σ_fogS, kPS)

    #Pₘ is computed over the kₚs
    P0N = 0.5 .* PPₘN .* (b_totkN.^2 .* I0N .+ 
                         2 .* b_totkN .* fPN .* I2N .+ 
                         fPN.^2 .* I4N) .+ NN
    P2N = 5 .* 0.5 .* PPₘN .* (0.5 .* b_totkN.^2 .* (3 .* I2N .- I0N) .+ 
                              b_totkN .* fPN .* (3 .* I4N .- I2N) .+ 
                              0.5 .* fPN.^2 .* (3 .* I6N .- I4N))
    P4N = 9 .* 0.5 .* PPₘN .* (0.125 .* b_totkN.^2 .* (35 .* I4N .- 30 .* I2N .+ 3 .* I0N) .+
                              0.25 .* b_totkN .* fPN .* (35 .* I6N .- 30 .* I4N .+ 3 .* I2N) .+
                              0.125 .* fPN.^2 .* (35 .* I8N .- 30 .* I6N .+ 3 .* I4N))

    P0S = 0.5 .* PPₘS .* (b_totkS.^2 .* I0S .+ 
                         2 .* b_totkS .* fPS .* I2S .+ 
                         fPS.^2 .* I4S) .+ NS
    P2S = 5 .* 0.5 .* PPₘS .* (0.5 .* b_totkS.^2 .* (3 .* I2S .- I0S) .+ 
                              b_totkS .* fPS .* (3 .* I4S .- I2S) .+ 
                              0.5 .* fPS.^2 .* (3 .* I6S .- I4S))
    P4S = 9 .* 0.5 .* PPₘS .* (0.125 .* b_totkS.^2 .* (35 .* I4S .- 30 .* I2S .+ 3 .* I0S) .+
                              0.25 .* b_totkS .* fPS .* (35 .* I6S .- 30 .* I4S .+ 3 .* I2S) .+
                              0.125 .* fPS.^2 .* (35 .* I8S .- 30 .* I6S .+ 3 .* I4S))

    #QN is a matrix k_eff x kₚ

    Q0N = QₗPN[:,:,1]
    Q2N = QₗPN[:,:,2]
    Q4N = QₗPN[:,:,3]

    Q0S = QₗPS[:,:,1]
    Q2S = QₗPS[:,:,2]
    Q4S = QₗPS[:,:,3]

    @tullio convolvedPkN[i] := Q0N[i,k] * P0N[i,k] + Q2N[i,k] * P2N[i,k] + Q4N[i,k] * P4N[i,k]
    @tullio convolvedPkS[i] := Q0S[i,k] * P0S[i,k] + Q2S[i,k] * P2S[i,k] + Q4S[i,k] * P4S[i,k]

    #IC
    Pof0N = convolvedPkN[1]
    convolvedPkN = front_cut_ks(start_kPN, convolvedPkN)
    predictionPN = convolvedPkN .- Pof0N .* W₀PN .- convolvedPkN .* WricPN

    Pof0S = convolvedPkS[1]
    convolvedPkS = front_cut_ks(start_kPS, convolvedPkS)
    predictionPS = convolvedPkS .- Pof0S .* W₀PS .- convolvedPkS .* WricPS

    #B likelihood
    B0N = compute_B0(B0_termsN, biasN, fBN, fNL, PsnBN)
    B0S = compute_B0(B0_termsS, biasS, fBS, fNL, PsnBS)
    BGICN = compute_B0(BGIC_termsN, biasN, fBN, fNL, PsnBN)
    BGICS = compute_B0(BGIC_termsS, biasS, fBS, fNL, PsnBS)

    #RIC contribution (GIC contribution in P0_model)
    predictionBN = BGICN .- B0N .* WricBN
    predictionBS = BGICS .- B0S .* WricBS

    predictionPBN = vcat(vec(predictionPN), predictionBN)
    predictionPBS = vcat(vec(predictionPS), predictionBS)
    predictionPBN_recast = iΓPBN * predictionPBN
    predictionPBS_recast = iΓPBS * predictionPBS

    dataPBN ~ MvNormal(predictionPBN_recast, I)
    dataPBS ~ MvNormal(predictionPBS_recast, I)

    return nothing
    
end