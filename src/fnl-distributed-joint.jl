using Turing
using Distributed
using Plots
using StatsPlots
using PairPlots
using MCMCChains
using DataFrames
using NPZ
using Optim

include("fnl-utils.jl")
include("fnl-model.jl")

addprocs(12)

#include
@everywhere begin
    using Turing
    using Distributed
    using Plots
    using StatsPlots
    using PairPlots
    using MCMCChains
    using DataFrames
    using NPZ
    using Optim
    
    include("fnl-utils.jl")
    include("fnl-model.jl")
end

#model instance
@everywhere begin
    #joint analysis
    p_weight = 3.0 #weights used for P(k) (0.0 = FKP)
    p = 3.0 #p to be used in the inference, if nothing it makes inference on f_NL*b_phi
    NN_weights = false #either to use the data from the linear weight catalog or the NN weight catalogue
    RIC = true #either to appy or not the RIC, GIC is always applied

    kₚ = get_kₚ(NN_weights)

    #NCG
    Pk_modelN = get_Pk_model("N", p_weight, NN_weights)
    αk_modelN = get_alphak_model("N", p_weight, NN_weights)
    fz_modelN = get_fz_model("N", p_weight, NN_weights)
    

    QₗN = get_Ql("N", p_weight, NN_weights)

    kN, Pk_dataN = get_Pk_data("N", p_weight, NN_weights)
    ΣN = get_Σ("N", p_weight, NN_weights)

    #if p_weight == 3.0
    #    Pk_dataN .*= -1.
    #end

    k_startN = findfirst(x -> x ≈ round(kN[1], digits=8), kₚ[:,1] .* 10. ./ 2)

    k_WN, W₀N = get_W₀k("N", p_weight, NN_weights)
    W₀kN = interpolate_fk(k_WN, W₀N./W₀N[1], kN)

    if RIC
        WricN = get_Wric("N", p_weight, NN_weights)
        #if p_weight == 3.0
        #    WricN .*= -1.
        #end
    end

    if check_k_dimension(kN, ΣN)
        len_k = length(kN)

        kₚN = cut_ks(len_k + k_startN - 1, kₚ)
        Pk_modelN = cut_ks(len_k + k_startN - 1, Pk_modelN)
        αk_modelN = cut_ks(len_k + k_startN - 1, αk_modelN)
        QₗN = cut_ks(len_k + k_startN - 1, QₗN)
        if RIC
            WricN = cut_ks(len_k, WricN)
        end
    else
        len_k = length(ΣN[:,1])

        Pk_dataN = cut_ks(len_k, Pk_dataN)
        W₀kN = cut_ks(len_k, W₀kN)
        kₚN = cut_ks(len_k + k_startN - 1, kₚ)
        Pk_modelN = cut_ks(len_k + k_startN - 1, Pk_modelN)
        αk_modelN = cut_ks(len_k + k_startN - 1, αk_modelN)
        QₗN = cut_ks(len_k + k_startN - 1, QₗN)
        if RIC
            WricN = cut_ks(len_k, WricN)
        end
    end

    len_k = length(ΣN[:,1])
    WHN = (1000 - len_k - 2) / (1000 - 1) #wishhart factor

    #SGC
    Pk_modelS = get_Pk_model("S", p_weight, NN_weights)
    αk_modelS = get_alphak_model("S", p_weight, NN_weights)
    fz_modelS = get_fz_model("S", p_weight, NN_weights)
    

    QₗS = get_Ql("S", p_weight, NN_weights)

    kS, Pk_dataS = get_Pk_data("S", p_weight, NN_weights)
    ΣS = get_Σ("S", p_weight, NN_weights)

    #if p_weight == 3.0
    #    Pk_dataS .*= -1.
    #end

    k_startS = findfirst(x -> x ≈ round(kS[1], digits=8), kₚ[:,1] .* 10. ./ 2)

    k_WS, W₀S = get_W₀k("S", p_weight, NN_weights)
    W₀kS = interpolate_fk(k_WS, W₀S./W₀S[1], kS)

    if RIC
        WricS = get_Wric("S", p_weight, NN_weights)
        #if p_weight == 3.0
        #    WricS .*= -1.
        #end
    end

    if check_k_dimension(kS, ΣS)
        len_k = length(kS)

        kₚS = cut_ks(len_k + k_startS - 1, kₚ)
        Pk_modelS = cut_ks(len_k + k_startS - 1, Pk_modelS)
        αk_modelS = cut_ks(len_k + k_startS - 1, αk_modelS)
        QₗS = cut_ks(len_k + k_startS - 1, QₗS)
        if RIC
            WricS = cut_ks(len_k, WricS)
        end
    else
        len_k = length(ΣS[:,1])

        Pk_dataS = cut_ks(len_k, Pk_dataS)
        W₀kS = cut_ks(len_k, W₀kS)
        kₚS = cut_ks(len_k + k_startS - 1, kₚ)
        Pk_modelS = cut_ks(len_k + k_startS - 1, Pk_modelS)
        αk_modelS = cut_ks(len_k + k_startS - 1, αk_modelS)
        QₗS = cut_ks(len_k + k_startS - 1, QₗS)
        if RIC
            WricS = cut_ks(len_k, WricS)
        end
    end

    len_k = length(ΣS[:,1])
    WHS = (1000 - len_k - 2) / (1000 - 1) #wishhart factor

    #Model
    if RIC
        data_model = P_qso_convolved_IC_joint(Pk_dataN, Pk_dataS, kₚN, kₚS, Pk_modelN, Pk_modelS, p, αk_modelN, αk_modelS, fz_modelN, fz_modelS, ΣN ./ WHN, ΣS ./ WHS, QₗN, QₗS, W₀kN, W₀kS, k_startN, k_startS, vec(WricN), vec(WricS))
    else
        data_model = P_qso_convolved_IC_joint(Pk_dataN, Pk_dataS, kₚN, kₚS, Pk_modelN, Pk_modelS, p, αk_modelN, αk_modelS, fz_modelN, fz_modelS, ΣN ./ WHN, ΣS ./ WHS, QₗN, QₗS, W₀kN, W₀kS, k_startN, k_startS)
    end
end

map = optimize(data_model, MAP())

sampler = NUTS(1000, 0.65)

chain = sample(data_model, sampler, MCMCDistributed(), 3000, 12, init_theta = map.values.array)

describe(chain)

if p == p_weight
    output_folder = get_output_folder("joint", p, NN_weights)
elseif p === nothing
    output_folder = get_output_folder("joint", p_weight, NN_weights)
else
    output_folder = get_output_folder("joint", p, p_weight, NN_weights)
end
    
pl = StatsPlots.plot(chain)
savefig(pl, joinpath(output_folder, "traces.png"))

df = DataFrame(chain)
par = df[:,[:f_nl, :b₁N, :σ_fogN, :NN, :b₁S, :σ_fogS, :NS]]

c_plot = PairPlots.corner(par)
savefig(c_plot, joinpath(output_folder, "corner.png"))

npzwrite(joinpath(output_folder, "chain.npy"), chain.value.data)
