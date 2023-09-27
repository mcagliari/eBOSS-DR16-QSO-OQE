include("fnl-model.jl")

function P_qso_joint(θ, kₚN, kₚS, PₘN, PₘS, p, α_kN, α_kS, fN, fS, QₗN, QₗS, W₀N, W₀S, start_kN, start_kS)
    f_nl = θ[1]
    b₁N = θ[2]
    σ_fogN = θ[3]
    NN = θ[4]
    b₁S = θ[5]
    σ_fogS = θ[6]
    NS = θ[7]

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


    #or do we need a realisation?
    return predictionN, predictionS

end

function P_qso_sample(θ, kₚ, Pₘ, p, α_k, f, Qₗ, W₀, start_k)
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
    prediction = convolvedPk #.- Pof0 .* W₀

    return prediction
end

function P_qso_sample(θ, kₚ, Pₘ, p, α_k, f, Qₗ, W₀, start_k, Wric)
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

    return prediction

end