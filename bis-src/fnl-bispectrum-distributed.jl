#B only analysis
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
    GC = "N" #or "S"
    p = 1.0 #or 1.6
    IC = true #or false
    window = true #or false 
    SN = false
    rebin = "" # or "-s2rebin"
    B_ind = SN ? 5 : 4

    #uploading models
    Pks_model = get_Pk_model(GC, window, rebin)
    Minv_model = get_Minv_model(GC, rebin)
    fz_model = get_fz_model(GC, rebin)

    #uploading data
    k1, k2, k3, Bk_data = get_Bk_data(GC, B_ind, rebin)
    Psn = get_S0(GC)
    Σ = get_Σ(GC, false, -1, SN, rebin)

    ks = [k1 k2 k3]

    if IC
        k_W, W₀ = get_W₀k(GC, -1, false)
        W₀k1 = interpolate_fk(k_W, W₀./W₀[1], k1)
        W₀k2 = interpolate_fk(k_W, W₀./W₀[1], k2)
        W₀k3 = interpolate_fk(k_W, W₀./W₀[1], k3)

        #GIC applyed to Pks_model
        P0_model = [Pks_model[1,1] Pks_model[1,2] Pks_model[1,3]]
        W₀ks = [W₀k1 W₀k2 W₀k3]
        Pks_model_GIC = Pks_model .- P0_model .* W₀ks
        
        Wric = get_Wric(GC, SN, rebin)
    end

    #cuts if necessary
    #it should be necessary to cut stuff only if the data are longer than the covariance
    if !check_k_dimension(k1, Σ)
        len_k = length(Σ[:,1])

        Pks_model = cut_ks(len_k, Pks_model)
        Minv_model = cut_ks(len_k, Minv_model)
        Bk_data = cut_ks(len_k, Bk_data)
        ks = cut_ks(len_k, ks)
        if IC
            Pks_model_GIC = cut_ks(len_k, Pks_model_GIC)
            Wric = cut_ks(len_k, Wric) #dovrebbe essere già della dimensine giusta
        end
    end
    
    maskB = (ks[:,3] .<= ks[:,1] .+ ks[:,2])
    
    #apply mask if presents
    if maskB != nothing
        Pks_model = mask_ks(maskB, Pks_model)
        Minv_model = mask_ks(maskB, Minv_model)
        Bk_data = mask_ks(maskB, Bk_data)
        ks = mask_ks(maskB, ks)
        if IC
            Pks_model_GIC = mask_ks(maskB, Pks_model_GIC)
            Wric = mask_ks(maskB, Wric)
        end
        Σ = mask_Σ(maskB, Σ)
        
        len_k = length(Σ[:,1])
    end

    Iabc = get_Is(ks)

    len_k = length(Σ[:,1])
    WH = (1000 - len_k - 2) / (1000 - 1) #wishhart factor
    
    println(len_k)
    println(WH)

    Σ_corrected = Σ ./ WH
    Γ = sqrt(Σ_corrected)
    iΓ = inv(Γ)
    Bk_recast = iΓ * Bk_data

    #B0_terms one for all computation
    println(size(ks), size(Pks_model), size(Minv_model), size(Iabc))
    terms = Bℓ(ks, Pks_model, Minv_model, Iabc, 0)
    nk = length(ks[:,1])
    Bα₂ = ones(nk)
    terms = vcat(terms, [Bα₂])
    B0_terms = reduce(vcat, terms')

    #compute B0 terms with GIC
    if IC
        terms = Bℓ(ks, Pks_model_GIC, Minv_model, Iabc, 0)
        nk = length(ks[:,1])
        Bα₂ = ones(nk)
        terms = vcat(terms, [Bα₂])
        BGIC_terms = reduce(vcat, terms')
    end

    if !IC
        println("No IC")
        data_model = B_qso(Bk_recast, B0_terms, fz_model, Psn, iΓ, p)
    else
        data_model = B_qso(Bk_recast, B0_terms, BGIC_terms, fz_model, Psn, iΓ, p, Wric)
    end
end

map = optimize(data_model, MAP())

sampler = NUTS(1000, 0.65)#, adtype=AutoZygote())

chain = sample(data_model, sampler, MCMCDistributed(), 2000, nprocs, init_theta = map.values.array)

describe(chain)

#outputs 
output_folder = get_output_folderB(GC, p, -1, "name", false, rebin)

pl = StatsPlots.plot(chain)
savefig(pl, joinpath(output_folder, "traces.png"))

df = DataFrame(chain)
#par = df[:,[:fNL, :b₁, :b₂, :bₛ, :c₁, :α₁, :α₂]]
par = df[:,[:fNL, :b₁, :b₂, :bₛ, :α₁, :α₂]]
#par = df[:,[:b₁, :b₂, :bₛ, :c₁, :α₁, :α₂]]
#par = df[:,[:b₁, :b₂, :bₛ, :α₁, :α₂]]

c_plot = PairPlots.pairplot(par)
save(joinpath(output_folder, "corner.png"), c_plot)

npzwrite(joinpath(output_folder, "chain.npy"), chain.value.data)
writedlm(joinpath(output_folder, "map.dat"), map.values.array)

#ll1 = generated_quantities(data_model, chain)
#writedlm(joinpath(output_folder, "loglikelihood.dat"), ll1)

if maskB != nothing
    writedlm(joinpath(output_folder, "mask.dat"), maskB)
end