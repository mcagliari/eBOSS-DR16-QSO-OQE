using NPZ

include("../src/fnl-model.jl")
include("fnl-bispectrum_model.jl")

#single field complete model
function B_qso_sample(θB, p, B0_terms, BGIC_terms, fB, WricB, PsnB)
    
    fNL = θB[1]
    b₁s = θB[2]
    b₂ = θB[3]
    bₛ = θB[4]
    c₁ = θB[5]
    α₁ = θB[6]
    α₂ = θB[7]
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕ = 2 * δ_c * (b₁s - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁s - 1)) - b₁s + 1) #a caso #usiamo Universal relation

    #parameter vectors
    bias = [b₁s, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂] #for B

    #B likelihood
    B0 = compute_B0(B0_terms, bias, fB, fNL, PsnB)
    BGIC = compute_B0(BGIC_terms, bias, fB, fNL, PsnB)

    #RIC contribution (GIC contribution in P0_model)
    predictionB = BGIC .- B0 .* WricB

    return predictionB
end

function B_qso_sample(θB, p, B0_terms, fB, PsnB)
    
    fNL = θB[1]
    b₁s = θB[2]
    b₂ = θB[3]
    bₛ = θB[4]
    c₁ = θB[5]
    α₁ = θB[6]
    α₂ = θB[7]
    #universal relation for B (bϕ is directly implemented in the P model, you do not have to provide it)
    δ_c = 1.686
    bϕ = 2 * δ_c * (b₁s - p) #~ Uniform(0.1, 6) #usiamo Universal relation
    bϕδ = bϕ + 2 * (δ_c * (b₂ - 8 / 21 * (b₁s - 1)) - b₁s + 1) #a caso #usiamo Universal relation

    #parameter vectors
    bias = [b₁s, b₂, bₛ, bϕ, bϕδ, c₁, α₁, α₂] #for B

    #B likelihood
    B0 = compute_B0(B0_terms, bias, fB, fNL, PsnB)

    return B0
end

function B_best_fit(θ, GC, p, p_weight, rebin::String="-s3rebin", window::Bool=true, IC_B::Bool=true)

    NN_weights = false
    SN = true
    B_ind = SN ? 5 : 4

    PkB_model = get_Pk_model(GC, window, rebin)
    MinvB_model = get_Minv_model(GC, rebin)
    fzB_model = get_fz_model(GC, rebin)

    kB1, kB2, kB3, Bk_data = get_Bk_data(GC, B_ind, rebin, NN_weights)
    PsnB = get_S0(GC) #I'm loading n^-1 nor n^-2!
    
    kBs = [kB1 kB2 kB3]

    k_W, W₀ = get_W₀k(GC, p_weight, NN_weights)
    W₀kB1 = interpolate_fk(k_W, W₀./W₀[1], kB1)
    W₀kB2 = interpolate_fk(k_W, W₀./W₀[1], kB2)
    W₀kB3 = interpolate_fk(k_W, W₀./W₀[1], kB3)

    if IC_B
        P0B_model = [PkB_model[1,1] PkB_model[1,2] PkB_model[1,3]]
        W₀kBs = [W₀kB1 W₀kB2 W₀kB3]
        PkB_model_GIC = PkB_model .- P0B_model .* W₀kBs
            
        WricB = get_Wric(GC, SN, rebin)
    end

    #covariance and cutting data if required
    PaB_dim = get_PaB_dimensions(GC, rebin)
    B_dim = PaB_dim[2]

    len_kB = length(kB1)

    if B_dim == len_kB
        PkB_model = cut_ks(B_dim, PkB_model)
        MinvB_model = cut_ks(B_dim, MinvB_model)
        Bk_data = cut_ks(B_dim, Bk_data)
        kBs = cut_ks(B_dim, kBs)
        if IC_B
            PkB_model_GIC = cut_ks(B_dim, PkB_model_GIC)
            WricB = cut_ks(B_dim, WricB) #dovrebbe essere già della dimensine giusta
        end
    end

    Iabc = get_Is(kBs)

    #B0_terms one for all computation
    terms = Bℓ(kBs, PkB_model, MinvB_model, Iabc, 0)
    nk = length(kBs[:,1])
    Bα₂ = ones(nk)
    terms = vcat(terms, [Bα₂])
    B0_terms = reduce(vcat, terms')
    
    if IC_B
        terms = Bℓ(kBs, PkB_model_GIC, MinvB_model, Iabc, 0)
        nk = length(kBs[:,1])
        Bα₂ = ones(nk)
        terms = vcat(terms, [Bα₂])
        BGIC_terms = reduce(vcat, terms')
    end

    #model
    if IC_B
        Bk = B_qso_sample(θ, p, B0_terms, BGIC_terms, fzB_model, WricB, PsnB)
    else
        Bk = B_qso_sample(θ, p, B0_terms, fzB_model, PsnB)
    end
    return [kB1, kB2, kB3, Bk]
end

function B_save_best_fit(input::String, GC::String, name::String, p_weight::Float64, p::Float64, output::String, rebin::String="-s3rebin", window::Bool=true, IC_B::Bool=true)
    log_density = 15
    p_weightf = p_weight == -1 ? "fkp" : p_weight
    p_folder = p === p_weightf ? "$p" : "$(p_weightf)/$p"
    p_folder = p === nothing ? "$(p_weightf)" : p_folder
    chain = joinpath(eboss_folder, input, "$(GC)GC", name, p_folder, "chain.npy")
    chn = npzread(chain)

    best_like = findmax(chn[:,log_density,:])
    nstep = best_like[2][1]
    nchn = best_like[2][2]

    θ = [chn[nstep,1,nchn], chn[nstep,3,nchn], chn[nstep,6,nchn], chn[nstep,7,nchn], chn[nstep,8,nchn], chn[nstep,9,nchn], chn[nstep,10,nchn]]
    Bk = B_best_fit(θ, GC, p, p_weight, rebin, window, IC_B)

    open(output; write=true) do f
        write(f, "# Best fit parameters fNL=$(θ[1]), b1B=$(θ[2]), b2=$(θ[3]), bs=$(θ[4]), c1=$(θ[5]), alpha1=$(θ[6]), alpha2=$(θ[7])")
        writedlm(f, transpose(Bk))
    end
end

function B_save_best_fit(input::String, name::String, p_weight::Float64, p::Float64, outputN::String, outputS::String, rebin::String="-s3rebin", window::Bool=true, IC_B::Bool=true)
    log_density = 24
    p_weightf = p_weight == -1 ? "fkp" : p_weight
    p_folder = p === p_weightf ? "$p" : "$(p_weightf)/$p"
    p_folder = p === nothing ? "$(p_weightf)" : p_folder
    chain = joinpath(eboss_folder, input, "joint", name, p_folder, "chain.npy")
    chn = npzread(chain)

    best_like = findmax(chn[:,log_density,:])
    nstep = best_like[2][1]
    nchn = best_like[2][2]

    θ_N = [chn[nstep,1,nchn], chn[nstep,3,nchn], chn[nstep,10,nchn], chn[nstep,11,nchn], chn[nstep,12,nchn], chn[nstep,13,nchn], chn[nstep,14,nchn]]
    Bk_N =  B_best_fit(θ_N, "N", p, p_weight, rebin, window, IC_B)
    open(outputN; write=true) do f
        write(f, "# Best fit parameters fNL=$(θ_N[1]), b1B=$(θ_N[2]), b2=$(θ_N[3]), bs=$(θ_N[4]), c1=$(θ_N[5]), alpha1=$(θ_N[6]), alpha2=$(θ_N[7])")
        writedlm(f, transpose(Bk_N))
    end

    θ_S = [chn[nstep,1,nchn], chn[nstep,5,nchn], chn[nstep,15,nchn], chn[nstep,16,nchn], chn[nstep,17,nchn], chn[nstep,18,nchn], chn[nstep,19,nchn]]
    Bk_S = B_best_fit(θ_S, "S", p, p_weight, rebin, window, IC_B)
    open(outputS; write=true) do f
        write(f, "# Best fit parameters fNL=$(θ_S[1]), b1B=$(θ_S[2]), b2=$(θ_S[3]), bs=$(θ_S[4]), c1=$(θ_S[5]), alpha1=$(θ_S[6]), alpha2=$(θ_S[7])")
        writedlm(f, transpose(Bk_S))
    end
end

function B_save_best_fit(fNL::Float64, input::String, GC::String, name::String, p_weight::Float64, p::Float64, output::String, rebin::String="-s3rebin", window::Bool=true, IC_B::Bool=true)
    log_density = 15
    p_weightf = p_weight == -1 ? "fkp" : p_weight
    p_folder = p === p_weightf ? "$p" : "$(p_weightf)/$p"
    p_folder = p === nothing ? "$(p_weightf)" : p_folder
    chain = joinpath(eboss_folder, input, "$(GC)GC", name, p_folder, "chain.npy")
    chn = npzread(chain)

    best_like = findmax(chn[:,log_density,:])
    nstep = best_like[2][1]
    nchn = best_like[2][2]

    θ = [fNL, chn[nstep,3,nchn], chn[nstep,6,nchn], chn[nstep,7,nchn], chn[nstep,8,nchn], chn[nstep,9,nchn], chn[nstep,10,nchn]]
    Bk = B_best_fit(θ, GC, p, p_weight, rebin, window, IC_B)

    open(output; write=true) do f
        write(f, "# Best fit parameters fNL=$(θ[1]), b1B=$(θ[2]), b2=$(θ[3]), bs=$(θ[4]), c1=$(θ[5]), alpha1=$(θ[6]), alpha2=$(θ[7])")
        writedlm(f, transpose(Bk))
    end
end


function B_save_best_fit(fNL::Float64, input::String, name::String, p_weight::Float64, p::Float64, outputN::String, outputS::String, rebin::String="-s3rebin", window::Bool=true, IC_B::Bool=true)
    log_density = 24
    p_weightf = p_weight == -1 ? "fkp" : p_weight
    p_folder = p === p_weightf ? "$p" : "$(p_weightf)/$p"
    p_folder = p === nothing ? "$(p_weightf)" : p_folder
    chain = joinpath(eboss_folder, input, "joint", name, p_folder, "chain.npy")
    chn = npzread(chain)

    best_like = findmax(chn[:,log_density,:])
    nstep = best_like[2][1]
    nchn = best_like[2][2]

    θ_N = [fNL, chn[nstep,3,nchn], chn[nstep,10,nchn], chn[nstep,11,nchn], chn[nstep,12,nchn], chn[nstep,13,nchn], chn[nstep,14,nchn]]
    Bk_N =  B_best_fit(θ_N, "N", p, p_weight, rebin, window, IC_B)
    open(outputN; write=true) do f
        write(f, "# Best fit parameters fNL=$(θ_N[1]), b1B=$(θ_N[2]), b2=$(θ_N[3]), bs=$(θ_N[4]), c1=$(θ_N[5]), alpha1=$(θ_N[6]), alpha2=$(θ_N[7])")
        writedlm(f, transpose(Bk_N))
    end

    θ_S = [fNL, chn[nstep,5,nchn], chn[nstep,15,nchn], chn[nstep,16,nchn], chn[nstep,17,nchn], chn[nstep,18,nchn], chn[nstep,19,nchn]]
    Bk_S = B_best_fit(θ_S, "S", p, p_weight, rebin, window, IC_B)
    open(outputS; write=true) do f
        write(f, "# Best fit parameters fNL=$(θ_S[1]), b1B=$(θ_S[2]), b2=$(θ_S[3]), bs=$(θ_S[4]), c1=$(θ_S[5]), alpha1=$(θ_S[6]), alpha2=$(θ_S[7])")
        writedlm(f, transpose(Bk_S))
    end
end