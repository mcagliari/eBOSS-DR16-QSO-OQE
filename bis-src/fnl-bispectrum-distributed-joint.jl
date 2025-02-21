#B only analysis joint
using Turing
using Distributed
using Plots
using StatsPlots
using PairPlots
using CairoMakie
using MCMCChains
using DataFrames
using Optim
using NPZ
using LinearAlgebra


include("../src/fnl-utils.jl")
include("fnl-bispectrum-utils.jl")
include("fnl-bispectrum_model.jl")

nprocs = Sys.CPU_THREADS
addprocs(nprocs)

#includes
@everywhere begin
    using Turing
    using Distributed
    using Plots
    using StatsPlots
    using PairPlots
    using CairoMakie
    using MCMCChains
    using DataFrames
    using Optim
    using NPZ
    using LinearAlgebra

    #Turing.setadbackend(:zygote)
    
    include("../src/fnl-utils.jl")
    include("fnl-bispectrum-utils.jl")
    include("fnl-bispectrum_model.jl")
end

#model instance
@everywhere begin
    #joint analysis
    p = 1.0 #or 1.6
    IC = true #or false
    window = true #or false not yet implemented in the reading

    #NGC
    GC = "N"
    #uploading models
    Pks_modelN = get_Pk_model(GC, window)
    Minv_modelN = get_Minv_model(GC)
    fz_modelN = get_fz_model(GC)

    #uploading data
    k1N, k2N, k3N, Bk_dataN = get_Bk_data(GC)
    PsnN = get_S0(GC)
    ΣN = get_Σ(GC, false)

    ksN = [k1N k2N k3N]

    if IC
        k_WN, W₀N = get_W₀k(GC, -1, false)
        W₀k1N = interpolate_fk(k_WN, W₀N./W₀N[1], k1N)
        W₀k2N = interpolate_fk(k_WN, W₀N./W₀N[1], k2N)
        W₀k3N = interpolate_fk(k_WN, W₀N./W₀N[1], k3N)

        #GIC applyed to Pks_model
        P0_modelN = [Pks_modelN[1,1] Pks_modelN[1,2] Pks_modelN[1,3]]
        W₀ksN = [W₀k1N W₀k2N W₀k3N]
        Pks_model_GICN = Pks_modelN .- P0_modelN .* W₀ksN
        
        WricN = get_Wric(GC)
    end

    #cuts if necessary
    #it should be necessary to cut stuff only if the data are longer than the covariance
    if !check_k_dimension(k1N, ΣN)
        len_k = length(ΣN[:,1])

        Pks_modelN = cut_ks(len_k, Pks_modelN)
        Minv_modelN = cut_ks(len_k, Minv_modelN)
        Bk_dataN = cut_ks(len_k, Bk_dataN)
        ksN = cut_ks(len_k, ksN)
        if IC
            Pks_model_GICN = cut_ks(len_k, Pks_model_GICN)
            WricN = cut_ks(len_k, WricN) #dovrebbe essere già della dimensine giusta
        end
    end

    IabcN = get_Is(ksN)

    len_kN = length(ΣN[:,1])
    WHN = (1000 - len_kN - 2) / (1000 - 1) #wishhart factor

    Σ_correctedN = ΣN ./ WHN
    ΓN = sqrt(Σ_correctedN)
    iΓN = inv(ΓN)
    Bk_recastN = iΓN * Bk_dataN

    #B0_terms one for all computation
    termsN = Bℓ(ksN, Pks_modelN, Minv_modelN, IabcN, 0)
    nkN = length(ksN[:,1])
    Bα₂N = ones(nkN)
    termsN = vcat(termsN, [Bα₂N])
    B0_termsN = reduce(vcat, termsN')

    #compute B0 terms with GIC
    if IC
        termsN = Bℓ(ksN, Pks_model_GICN, Minv_modelN, IabcN, 0)
        nkN = length(ksN[:,1])
        Bα₂N = ones(nkN)
        termsN = vcat(termsN, [Bα₂N])
        BGIC_termsN = reduce(vcat, termsN')
    end

    #SGC
    GC = "S"
    #uploading models
    Pks_modelS = get_Pk_model(GC, window)
    Minv_modelS = get_Minv_model(GC)
    fz_modelS = get_fz_model(GC)

    #uploading data
    k1S, k2S, k3S, Bk_dataS = get_Bk_data(GC)
    PsnS = get_S0(GC)
    ΣS = get_Σ(GC, false)

    ksS = [k1S k2S k3S]

    if IC
        k_WS, W₀S = get_W₀k(GC, -1, false)
        W₀k1S = interpolate_fk(k_WS, W₀S./W₀S[1], k1S)
        W₀k2S = interpolate_fk(k_WS, W₀S./W₀S[1], k2S)
        W₀k3S = interpolate_fk(k_WS, W₀S./W₀S[1], k3S)

        #GIC applyed to Pks_model
        P0_modelS = [Pks_modelS[1,1] Pks_modelS[1,2] Pks_modelS[1,3]]
        W₀ksS = [W₀k1S W₀k2S W₀k3S]
        Pks_model_GICS = Pks_modelS .- P0_modelS .* W₀ksS
        
        WricS = get_Wric(GC)
    end

    #cuts if necessary
    #it should be necessary to cut stuff only if the data are longer than the covariance
    if !check_k_dimension(k1S, ΣS)
        len_k = length(ΣS[:,1])

        Pks_modelS = cut_ks(len_k, Pks_modelS)
        Minv_modelS = cut_ks(len_k, Minv_modelS)
        Bk_dataS = cut_ks(len_k, Bk_dataS)
        ksS = cut_ks(len_k, ksS)
        if IC
            Pks_model_GICS = cut_ks(len_k, Pks_model_GICS)
            WricS = cut_ks(len_k, WricS) #dovrebbe essere già della dimensine giusta
        end
    end

    IabcS = get_Is(ksS)

    len_kS = length(ΣS[:,1])
    WHS = (1000 - len_kS - 2) / (1000 - 1) #wishhart factor

    Σ_correctedS = ΣS ./ WHS
    ΓS = sqrt(Σ_correctedS)
    iΓS = inv(ΓS)
    Bk_recastS = iΓS * Bk_dataS

    #B0_terms one for all computation
    termsS = Bℓ(ksS, Pks_modelS, Minv_modelS, IabcS, 0)
    nkS = length(ksS[:,1])
    Bα₂S = ones(nkS)
    termsS = vcat(termsS, [Bα₂S])
    B0_termsS = reduce(vcat, termsS')

    #compute B0 terms with GIC
    if IC
        termsS = Bℓ(ksS, Pks_model_GICS, Minv_modelS, IabcS, 0)
        nkS = length(ksS[:,1])
        Bα₂S = ones(nkS)
        termsS = vcat(termsS, [Bα₂S])
        BGIC_termsS = reduce(vcat, termsS')
    end

    if !IC
        data_model = B_qso(Bk_recastN, Bk_recastS, p, B0_termsN, B0_termsS, fz_modelN, fz_modelS, PsnN, PsnS, iΓN, iΓS)
    else
        data_model = B_qso(Bk_recastN, Bk_recastS, p, B0_termsN, B0_termsS, BGIC_termsN, BGIC_termsS, fz_modelN, fz_modelS, PsnN, PsnS, iΓN, iΓS, WricN, WricS)
    end
end

map = optimize(data_model, MAP())

sampler = NUTS(10, 0.65)

chain = sample(data_model, sampler, MCMCDistributed(), 30, nprocs, init_theta = map.values.array)

describe(chain)

#outputs 
output_folder = get_output_folderB("joint", p, -1, false)

pl = StatsPlots.plot(chain)
savefig(pl, joinpath(output_folder, "traces.png"))

df = DataFrame(chain)
par = df[:,[:fNL, :b₁N, :b₂N, :bₛN, :c₁N, :α₁N, :α₂N, :b₁S, :b₂S, :bₛS, :c₁S, :α₁S, :α₂S]]

c_plot = PairPlots.pairplot(par)
save(joinpath(output_folder, "corner.png"), c_plot)

npzwrite(joinpath(output_folder, "chain.npy"), chain.value.data)
writedlm(joinpath(output_folder, "map.dat"), map.values.array)