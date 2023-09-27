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
    #parsed_args = parse_commandline()
    GC = "N" #parsed_args["sample"] #Galaxy cap to analyse
    p_weight = 1.0 #weights used for P(k) (0.0 = FKP)
    p = 1.0 #parsed_args["p"] #p to be used in the inference, if nothing it makes inference on f_NL*b_phi
    NN_weights = true #either to use the data from the linear weight catalog or the NN weight catalogue
    RIC = true #either to appy or not the RIC, GIC is always applied

    Pk_model = get_Pk_model(GC, p_weight, NN_weights)
    αk_model = get_alphak_model(GC, p_weight, NN_weights)
    fz_model = get_fz_model(GC, p_weight, NN_weights)
    kₚ = get_kₚ(NN_weights)

    Qₗ = get_Ql(GC, p_weight, NN_weights)

    k, Pk_data = get_Pk_data(GC, p_weight, NN_weights)
    Σ = get_Σ(GC, p_weight, NN_weights)

    #if p_weight == 3.0
    #    Pk_data .*= -1.
    #end

    k_start = findfirst(x -> x ≈ round(k[1], digits=8), kₚ[:,1] .* 10. ./ 2)

    k_W, W₀ = get_W₀k(GC, p_weight, NN_weights)
    W₀k = interpolate_fk(k_W, W₀./W₀[1], k)

    if RIC
        Wric = get_Wric(GC, p_weight, NN_weights)
        #if p_weight == 3.0
        #    Wric .*= -1.
        #end
    end
    

    if check_k_dimension(k, Σ)
        len_k = length(k)

        kₚ = cut_ks(len_k + k_start - 1, kₚ)
        Pk_model = cut_ks(len_k + k_start - 1, Pk_model)
        αk_model = cut_ks(len_k + k_start - 1, αk_model)
        Qₗ = cut_ks(len_k + k_start - 1, Qₗ)
        if RIC
            Wric = cut_ks(len_k, Wric)
        end
    else
        len_k = length(Σ[:,1])

        Pk_data = cut_ks(len_k, Pk_data)
        W₀k = cut_ks(len_k, W₀k)
        kₚ = cut_ks(len_k + k_start - 1, kₚ)
        Pk_model = cut_ks(len_k + k_start - 1, Pk_model)
        αk_model = cut_ks(len_k + k_start - 1, αk_model)
        Qₗ = cut_ks(len_k + k_start - 1, Qₗ)
        if RIC
            Wric = cut_ks(len_k, Wric)
        end
    end

    len_k = length(Σ[:,1])
    WH = (1000 - len_k - 2) / (1000 - 1) #wishhart factor

    if RIC
        data_model = P_qso_convolved_IC(Pk_data, kₚ, Pk_model, p, αk_model, fz_model, Σ ./ WH, Qₗ, W₀k, k_start, vec(Wric))
    else
        data_model = P_qso_convolved_IC(Pk_data, kₚ, Pk_model, p, αk_model, fz_model, Σ ./ WH, Qₗ, W₀k, k_start)
    end
end

map = optimize(data_model, MAP())

sampler = NUTS(1000, 0.65)

chain = sample(data_model, sampler, MCMCDistributed(), 3000, 12, init_theta = map.values.array)

describe(chain)

if p == p_weight
    output_folder = get_output_folder(GC, p, NN_weights)
elseif p === nothing
    output_folder = get_output_folder(GC, p_weight, NN_weights)
else
    output_folder = get_output_folder(GC, p, p_weight, NN_weights)
end
    
pl = StatsPlots.plot(chain)
savefig(pl, joinpath(output_folder, "traces.png"))

df = DataFrame(chain)
par = df[:,[:f_nl, :b₁, :σ_fog, :N]]

c_plot = PairPlots.corner(par)
savefig(c_plot, joinpath(output_folder, "corner.png"))

npzwrite(joinpath(output_folder, "chain.npy"), chain.value.data)
