#P+B analysis
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
    GC = "N" #or "N"
    p = 1.0 #or 1.6 #p for model
    p_weight = 1.0 #weights used for P(k)
    NN_weights = false #or true
    window = true #or false #window convolution can be turned off for B(k) only
    SN = true
    rebin = "-s3rebin" # or ""
    IC_B = true

    B_ind = SN ? 5 : 4

    #load data and models
    #POWER SPECTRUM
    PkP_model = get_Pk_model(GC, p_weight, NN_weights)
    αkP_model = get_alphak_model(GC, p_weight, NN_weights)
    fzP_model = get_fz_model(GC, p_weight, NN_weights)
    kₚP = get_kₚ(NN_weights)

    QₗP = get_Ql(GC, p_weight, NN_weights)

    kP, Pk_data = get_Pk_data(GC, p_weight, NN_weights)

    k_start = findfirst(x -> x ≈ round(kP[1], digits=8), kₚP[:,1] .* 10. ./ 2)

    k_W, W₀ = get_W₀k(GC, p_weight, NN_weights)
    W₀kP = interpolate_fk(k_W, W₀./W₀[1], kP)

    WricP = get_Wric(GC, p_weight, NN_weights)

    #BISPECTRUM

    PkB_model = get_Pk_model(GC, window, rebin)
    MinvB_model = get_Minv_model(GC, rebin)
    fzB_model = get_fz_model(GC, rebin)

    kB1, kB2, kB3, Bk_data = get_Bk_data(GC, B_ind, rebin, NN_weights)
    PsnB = get_S0(GC) #I'm loading n^-1 nor n^-2!
    
    kBs = [kB1 kB2 kB3]

    W₀kB1 = interpolate_fk(k_W, W₀./W₀[1], kB1)
    W₀kB2 = interpolate_fk(k_W, W₀./W₀[1], kB2)
    W₀kB3 = interpolate_fk(k_W, W₀./W₀[1], kB3)

    #GIC applyed to Pks_model
    if IC_B
        P0B_model = [PkB_model[1,1] PkB_model[1,2] PkB_model[1,3]]
        W₀kBs = [W₀kB1 W₀kB2 W₀kB3]
        PkB_model_GIC = PkB_model .- P0B_model .* W₀kBs
        
        WricB = get_Wric(GC, SN, rebin)
    end

    #covariance and cutting data if required
    Σ = get_Σ(GC, true, p_weight, SN, rebin) #P+B covariance
    len_k = length(Σ[:,1])

    PaB_dim = get_PaB_dimensions(GC, rebin)
    P_dim = PaB_dim[1]
    B_dim = PaB_dim[2]

    len_kP = length(kP)
    len_kB = length(kB1)

    if (P_dim + B_dim) == (len_kP + len_kB)
        #cut only P(k) model
        kₚP = cut_ks(P_dim + k_start - 1, kₚP)
        PkP_model = cut_ks(P_dim + k_start - 1, PkP_model)
        αkP_model = cut_ks(P_dim + k_start - 1, αkP_model)
        QₗP = cut_ks(P_dim + k_start - 1, QₗP)
        WricP = cut_ks(P_dim, WricP)
    else
        #cut: P(k) model, P(k) data, B(k) model, and B(k) data
        Pk_data = cut_ks(P_dim, Pk_data)
        W₀kP = cut_ks(P_dim, W₀kP)

        kₚP = cut_ks(P_dim + k_start - 1, kₚP)
        PkP_model = cut_ks(P_dim + k_start - 1, PkP_model)
        αkP_model = cut_ks(P_dim + k_start - 1, αkP_model)
        QₗP = cut_ks(P_dim + k_start - 1, QₗP)
        WricP = cut_ks(P_dim, WricP)

        PkB_model = cut_ks(B_dim, PkB_model)
        MinvB_model = cut_ks(B_dim, MinvB_model)
        Bk_data = cut_ks(B_dim, Bk_data)
        kBs = cut_ks(B_dim, kBs)
        if IC_B
            PkB_model_GIC = cut_ks(B_dim, PkB_model_GIC)
            WricB = cut_ks(B_dim, WricB) #dovrebbe essere già della dimensine giusta
        end
    end
    
    maskB = nothing #you can use this to set a mask on the triangles
    
    #apply mask if presents
    if maskB != nothing
        maskP = trues(P_dim)
        mask = vcat(maskP, maskB)
        
        PkB_model = mask_ks(maskB, PkB_model)
        MinvB_model = mask_ks(maskB, MinvB_model)
        Bk_data = mask_ks(maskB, Bk_data)
        kBs = mask_ks(maskB, kBs)
        if IC_B
            PkB_model_GIC = mask_ks(maskB, PkB_model_GIC)
            WricB = mask_ks(maskB, WricB)
        end
        
        Σ = mask_Σ(mask, Σ)
        
        len_k = length(Σ[:,1])
    end
    
    WH_PB = (1000 - len_k - 2) / (1000 - 1) #wishhart factor
    Σ_corrected = Σ ./ WH_PB

    Γ = sqrt(Σ_corrected)
    iΓ = inv(Γ)

    Iabc = get_Is(kBs)

    #B0_terms one for all computation
    terms = Bℓ(kBs, PkB_model, MinvB_model, Iabc, 0)
    nk = length(kBs[:,1])
    Bα₂ = ones(nk)
    terms = vcat(terms, [Bα₂])
    B0_terms = reduce(vcat, terms')

    #compute B0 terms with GIC
    if IC_B
        terms = Bℓ(kBs, PkB_model_GIC, MinvB_model, Iabc, 0)
        nk = length(kBs[:,1])
        Bα₂ = ones(nk)
        terms = vcat(terms, [Bα₂])
        BGIC_terms = reduce(vcat, terms')
    end
    
    
    #concatenating datavector
    PpB_data = vcat(Pk_data, Bk_data)
    PpB_recast = iΓ * PpB_data

    #model
    if IC_B
        data_model = PpB_qso_convolved_IC(PpB_recast, p, iΓ, kₚP, PkP_model, αkP_model, B0_terms, BGIC_terms, fzP_model, fzB_model, QₗP, W₀kP, k_start, WricP, WricB, PsnB, true)
    else
        data_model = PpB_qso_convolved(PpB_recast, p, iΓ, kₚP, PkP_model, αkP_model, B0_terms, fzP_model, fzB_model, QₗP, W₀kP, k_start, WricP, PsnB, true)
    end
end

map = optimize(data_model, MAP())

sampler = NUTS(1000, 0.65)#, adtype=AutoZygote())

chain = sample(data_model, sampler, MCMCDistributed(), 3000, nprocs, init_theta = map.values.array)

describe(chain)

#outputs 

if p == p_weight #optimal weights for P
    output_folder = get_output_folderB(GC, p, "data-SNsub", true, rebin, NN_weights)
else #FKP weights for P
    output_folder = get_output_folderB(GC, p, p_weight, "data-SNsub", true, rebin, NN_weights)
end

pl = StatsPlots.plot(chain)
savefig(pl, joinpath(output_folder, "traces.png"))

df = DataFrame(chain)
par = df[:,["fNL", "b₁s[1]", "σ_fog", "N", "b₁s[2]", "b₂", "bₛ", "c₁", "α₁", "α₂"]]
#par = df[:,["fNL", "b₁s[1]", "σ_fog", "N", "b₁s[2]", "b₂", "bₛ", "α₁", "α₂"]]

c_plot = PairPlots.pairplot(par)
save(joinpath(output_folder, "corner.png"), c_plot)

npzwrite(joinpath(output_folder, "chain.npy"), chain.value.data)
writedlm(joinpath(output_folder, "map.dat"), map.values.array)

if maskB != nothing
    writedlm(joinpath(output_folder, "mask.dat"), maskB)
end
