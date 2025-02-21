#P+B joint analysis
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
include("fnl-P_plus_B-model.jl")

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
    
    include("../src/fnl-utils.jl")
    include("fnl-bispectrum-utils.jl")
    include("fnl-P_plus_B-model.jl")
end

#model instance
@everywhere begin
    #info run
    #joint
    p = 1.0 #or 1.6 #p for model
    p_weight = 1.0 #weights used for P(k)
    NN_weights = false #or false
    window = true #or false #window convolution can be turned off for B(k) only
    SN = true
    rebin = "-s3rebin" # or ""
    #the IC correction is always applyed right now
    B_ind = SN ? 5 : 4

    #NGC
    GC = "N"
    #load data and models
    #POWER SPECTRUM
    PkP_modelN = get_Pk_model(GC, p_weight, NN_weights)
    αkP_modelN = get_alphak_model(GC, p_weight, NN_weights)
    fzP_modelN = get_fz_model(GC, p_weight, NN_weights)
    kₚPN = get_kₚ(NN_weights)

    QₗPN = get_Ql(GC, p_weight, NN_weights)

    kPN, Pk_dataN = get_Pk_data(GC, p_weight, NN_weights)

    k_startN = findfirst(x -> x ≈ round(kPN[1], digits=8), kₚPN[:,1] .* 10. ./ 2)

    k_WN, W₀N = get_W₀k(GC, p_weight, NN_weights)
    W₀kPN = interpolate_fk(k_WN, W₀N./W₀N[1], kPN)

    WricPN = get_Wric(GC, p_weight, NN_weights)

    #BISPECTRUM

    PkB_modelN = get_Pk_model(GC, window, rebin)
    MinvB_modelN = get_Minv_model(GC, rebin)
    fzB_modelN = get_fz_model(GC, rebin)

    kB1N, kB2N, kB3N, Bk_dataN = get_Bk_data(GC, B_ind, rebin, NN_weights)
    PsnBN = get_S0(GC)
    
    kBsN = [kB1N kB2N kB3N]

    W₀kB1N = interpolate_fk(k_WN, W₀N./W₀N[1], kB1N)
    W₀kB2N = interpolate_fk(k_WN, W₀N./W₀N[1], kB2N)
    W₀kB3N = interpolate_fk(k_WN, W₀N./W₀N[1], kB3N)

    #GIC applyed to Pks_model
    P0B_modelN = [PkB_modelN[1,1] PkB_modelN[1,2] PkB_modelN[1,3]]
    W₀kBsN = [W₀kB1N W₀kB2N W₀kB3N]
    PkB_model_GICN = PkB_modelN .- P0B_modelN .* W₀kBsN
    
    WricBN = get_Wric(GC, SN, rebin)
    
    #covariance and cutting data if required
    ΣN = get_Σ(GC, true, p_weight, SN, rebin) #P+B covariance
    len_kN = length(ΣN[:,1])

    PaB_dimN = get_PaB_dimensions(GC, rebin)
    P_dimN = PaB_dimN[1]
    B_dimN = PaB_dimN[2]

    len_kPN = length(kPN)
    len_kBN = length(kB1N)

    if (P_dimN + B_dimN) == (len_kPN + len_kBN)
        #cut only P(k) model
        kₚPN = cut_ks(P_dimN + k_startN - 1, kₚPN)
        PkP_modelN = cut_ks(P_dimN + k_startN - 1, PkP_modelN)
        αkP_modelN = cut_ks(P_dimN + k_startN - 1, αkP_modelN)
        QₗPN = cut_ks(P_dimN + k_startN - 1, QₗPN)
        WricPN = cut_ks(P_dimN, WricPN)
    else
        #cut: P(k) model, P(k) data, B(k) model, and B(k) data
        Pk_dataN = cut_ks(P_dimN, Pk_dataN)
        W₀kPN = cut_ks(P_dimN, W₀kPN)

        kₚPN = cut_ks(P_dimN + k_startN - 1, kₚPN)
        PkP_modelN = cut_ks(P_dimN + k_startN - 1, PkP_modelN)
        αkP_modelN = cut_ks(P_dimN + k_startN - 1, αkP_modelN)
        QₗPN = cut_ks(P_dimN + k_startN - 1, QₗPN)
        WricPN = cut_ks(P_dimN, WricPN)

        PkB_modelN = cut_ks(B_dimN, PkB_modelN)
        MinvB_modelN = cut_ks(B_dimN, MinvB_modelN)
        Bk_dataN = cut_ks(B_dimN, Bk_dataN)
        kBsN = cut_ks(B_dimN, kBsN)
        PkB_model_GICN = cut_ks(B_dimN, PkB_model_GICN)
        WricBN = cut_ks(B_dimN, WricBN) #dovrebbe essere già della dimensine giusta
    end
    
    kBminN = minimum(kBsN[:,1])
    maskBN = nothing #(kBsN[:,1] .> kBminN)
    
    #apply mask if presents
    if maskBN != nothing
        maskPN = trues(P_dimN)
        maskN = vcat(maskPN, maskBN)
        
        PkB_modelN = mask_ks(maskBN, PkB_modelN)
        MinvB_modelN = mask_ks(maskBN, MinvB_modelN)
        Bk_dataN = mask_ks(maskBN, Bk_dataN)
        kBsN = mask_ks(maskBN, kBsN)
        PkB_model_GICN = mask_ks(maskBN, PkB_model_GICN)
        WricBN = mask_ks(maskBN, WricBN)
        
        ΣN = mask_Σ(maskN, ΣN)
        
        len_kN = length(ΣN[:,1])
    end
    
    WH_PBN = (1000 - len_kN - 2) / (1000 - 1) #wishhart factor
    Σ_correctedN = ΣN ./ WH_PBN
    
    ΓN = sqrt(Σ_correctedN)
    iΓN = inv(ΓN)

    IabcN = get_Is(kBsN)

    #B0_terms one for all computation
    termsN = Bℓ(kBsN, PkB_modelN, MinvB_modelN, IabcN, 0)
    nkN = length(kBsN[:,1])
    Bα₂N = ones(nkN)
    termsN = vcat(termsN, [Bα₂N])
    B0_termsN = reduce(vcat, termsN')

    #compute B0 terms with GIC
    termsN = Bℓ(kBsN, PkB_model_GICN, MinvB_modelN, IabcN, 0)
    nkN = length(kBsN[:,1])
    Bα₂N = ones(nkN)
    termsN = vcat(termsN, [Bα₂N])
    BGIC_termsN = reduce(vcat, termsN')

    #concatenating datavector
    PpB_dataN = vcat(Pk_dataN, Bk_dataN)
    PpB_recastN = iΓN * PpB_dataN

    #SGC
    GC = "S"
    #load data and models
    #POWER SPECTRUM
    PkP_modelS = get_Pk_model(GC, p_weight, NN_weights)
    αkP_modelS = get_alphak_model(GC, p_weight, NN_weights)
    fzP_modelS = get_fz_model(GC, p_weight, NN_weights)
    kₚPS = get_kₚ(NN_weights)

    QₗPS = get_Ql(GC, p_weight, NN_weights)

    kPS, Pk_dataS = get_Pk_data(GC, p_weight, NN_weights)

    k_startS = findfirst(x -> x ≈ round(kPS[1], digits=8), kₚPS[:,1] .* 10. ./ 2)

    k_WS, W₀S = get_W₀k(GC, p_weight, NN_weights)
    W₀kPS = interpolate_fk(k_WS, W₀S./W₀S[1], kPS)

    WricPS = get_Wric(GC, p_weight, NN_weights)

    #BISPECTRUM

    PkB_modelS = get_Pk_model(GC, window, rebin)
    MinvB_modelS = get_Minv_model(GC, rebin)
    fzB_modelS = get_fz_model(GC, rebin)

    kB1S, kB2S, kB3S, Bk_dataS = get_Bk_data(GC, B_ind, rebin, NN_weights)
    PsnBS = get_S0(GC)
    
    kBsS = [kB1S kB2S kB3S]

    W₀kB1S = interpolate_fk(k_WS, W₀S./W₀S[1], kB1S)
    W₀kB2S = interpolate_fk(k_WS, W₀S./W₀S[1], kB2S)
    W₀kB3S = interpolate_fk(k_WS, W₀S./W₀S[1], kB3S)

    #GIC applyed to Pks_model
    P0B_modelS = [PkB_modelS[1,1] PkB_modelS[1,2] PkB_modelS[1,3]]
    W₀kBsS = [W₀kB1S W₀kB2S W₀kB3S]
    PkB_model_GICS = PkB_modelS .- P0B_modelS .* W₀kBsS
    
    WricBS = get_Wric(GC, SN, rebin)

    #covariance and cutting data if required
    ΣS = get_Σ(GC, true, p_weight, SN, rebin) #P+B covariance
    len_kS = length(ΣS[:,1])

    PaB_dimS = get_PaB_dimensions(GC, rebin)
    P_dimS = PaB_dimS[1]
    B_dimS = PaB_dimS[2]

    len_kPS = length(kPS)
    len_kBS = length(kB1S)

    if (P_dimS + B_dimS) == (len_kPS + len_kBS)
        #cut only P(k) model
        kₚPS = cut_ks(P_dimS + k_startS - 1, kₚPS)
        PkP_modelS = cut_ks(P_dimS + k_startS - 1, PkP_modelS)
        αkP_modelS = cut_ks(P_dimS + k_startS - 1, αkP_modelS)
        QₗPS = cut_ks(P_dimS + k_startS - 1, QₗPS)
        WricPS = cut_ks(P_dimS, WricPS)
    else
        #cut: P(k) model, P(k) data, B(k) model, and B(k) data
        Pk_dataS = cut_ks(P_dimS, Pk_dataS)
        W₀kPS = cut_ks(P_dimS, W₀kPS)

        kₚPS = cut_ks(P_dimS + k_startS - 1, kₚPS)
        PkP_modelS = cut_ks(P_dimS + k_startS - 1, PkP_modelS)
        αkP_modelS = cut_ks(P_dimS + k_startS - 1, αkP_modelS)
        QₗPS = cut_ks(P_dimS + k_startS - 1, QₗPS)
        WricPS = cut_ks(P_dimS, WricPS)

        PkB_modelS = cut_ks(B_dimS, PkB_modelS)
        MinvB_modelS = cut_ks(B_dimS, MinvB_modelS)
        Bk_dataS = cut_ks(B_dimS, Bk_dataS)
        kBsS = cut_ks(B_dimS, kBsS)
        PkB_model_GICS = cut_ks(B_dimS, PkB_model_GICS)
        WricBS = cut_ks(B_dimS, WricBS) #dovrebbe essere già della dimensine giusta
    end

    kBminS = minimum(kBsS[:,1])
    maskBS = nothing #(kBsS[:,1] .> kBminS)
    
    #apply mask if presents
    if maskBS != nothing
        maskPS = trues(P_dimS)
        maskS = vcat(maskPS, maskBS)
        
        PkB_modelS = mask_ks(maskBS, PkB_modelS)
        MinvB_modelS = mask_ks(maskBS, MinvB_modelS)
        Bk_dataS = mask_ks(maskBS, Bk_dataS)
        kBsS = mask_ks(maskBS, kBsS)
        PkB_model_GICS = mask_ks(maskBS, PkB_model_GICS)
        WricBS = mask_ks(maskBS, WricBS)
        
        ΣS = mask_Σ(maskS, ΣS)
        
        len_kS = length(ΣS[:,1])
    end
    
    WH_PBS = (1000 - len_kS - 2) / (1000 - 1) #wishhart factor
    Σ_correctedS = ΣS ./ WH_PBS
   
    ΓS = sqrt(Σ_correctedS)
    iΓS = inv(ΓS)
    
    IabcS = get_Is(kBsS)

    #B0_terms one for all computation
    termsS = Bℓ(kBsS, PkB_modelS, MinvB_modelS, IabcS, 0)
    nkS = length(kBsS[:,1])
    Bα₂S = ones(nkS)
    termsS = vcat(termsS, [Bα₂S])
    B0_termsS = reduce(vcat, termsS')

    #compute B0 terms with GIC
    termsS = Bℓ(kBsS, PkB_model_GICS, MinvB_modelS, IabcS, 0)
    nkS = length(kBsS[:,1])
    Bα₂S = ones(nkS)
    termsS = vcat(termsS, [Bα₂S])
    BGIC_termsS = reduce(vcat, termsS')

    #concatenating datavector
    PpB_dataS = vcat(Pk_dataS, Bk_dataS)
    PpB_recastS = iΓS * PpB_dataS

    #model
    data_model = PpB_qso_convolved_IC(PpB_recastN, PpB_recastS, p, iΓN, iΓS, kₚPN, kₚPS, PkP_modelN, PkP_modelS, αkP_modelN, αkP_modelS, B0_termsN, B0_termsS, BGIC_termsN, BGIC_termsS, fzP_modelN, fzP_modelS, fzB_modelN, fzB_modelS, QₗPN, QₗPS, W₀kPN, W₀kPS, k_startN, k_startS, WricPN, WricPS, WricBN, WricBS, PsnBN, PsnBS, true)
end

#map = optimize(data_model, MAP())

sampler = NUTS(1000, 0.65)

chain = sample(data_model, sampler, MCMCDistributed(), 3000, nprocs)#, init_theta = map.values.array)

describe(chain)

#outputs 
if p == p_weight #optimal weights for P
    output_folder = get_output_folderB("joint", p, "data-SNsub", true, rebin, NN_weights)
else #FKP weights for P
    output_folder = get_output_folderB("joint", p, p_weight, "data-SNsub", true, rebin, NN_weights)
end

pl = StatsPlots.plot(chain)
savefig(pl, joinpath(output_folder, "traces.png"))

df = DataFrame(chain)
par = df[:,["fNL", "b₁sN[1]", "σ_fogN", "NN", "b₁sN[2]", "b₂N", "bₛN", "c₁N", "α₁N", "α₂N", "b₁sS[1]", "σ_fogS", "NS", "b₁sS[2]", "b₂S", "bₛS", "c₁S", "α₁S", "α₂S"]]
#par = df[:,["fNL", "b₁sN[1]", "σ_fogN", "NN", "b₁sN[2]", "b₂N", "bₛN", "α₁N", "α₂N", "b₁sS[1]", "σ_fogS", "NS", "b₁sS[2]", "b₂S", "bₛS", "α₁S", "α₂S"]]

c_plot = PairPlots.pairplot(par)
save(joinpath(output_folder, "corner.png"), c_plot)

npzwrite(joinpath(output_folder, "chain.npy"), chain.value.data)
#writedlm(joinpath(output_folder, "map.dat"), map.values.array)

if maskBN != nothing
    writedlm(joinpath(output_folder, "maskN.dat"), maskBN)
end
if maskBS != nothing
    writedlm(joinpath(output_folder, "maskS.dat"), maskBS)
end