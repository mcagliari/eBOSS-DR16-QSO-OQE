using NPZ
include("post-model.jl")

function confidence_interval(chain, cl)
    d = sort(chain)
    n = length(chain)
    
    n_samples = floor(Int, n * cl)
    
    int_width = d[n_samples+1:end] .- d[1:n-n_samples]
    min_int = argmin(int_width)
    
    return [d[min_int], d[min_int+n_samples]]
end

function best_fit(θ, GC, p_weight, p, NN_weights, RIC)

    Pk_model = get_Pk_model(GC, p_weight, NN_weights)
    αk_model = get_alphak_model(GC, p_weight, NN_weights)
    fz_model = get_fz_model(GC, p_weight, NN_weights)
    kₚ = get_kₚ(NN_weights)

    Qₗ = get_Ql(GC, p_weight, NN_weights)

    k, Pk_data = get_Pk_data(GC, p_weight, NN_weights)
    Σ = get_Σ(GC, p_weight, NN_weights)

    if p_weight == 3.0
        Pk_data .*= -1.
    end

    k_start = findfirst(x -> x ≈ round(k[1], digits=8), kₚ[:,1] .* 10. ./ 2)

    k_W, W₀ = get_W₀k(GC, p_weight, NN_weights)
    W₀k = interpolate_fk(k_W, W₀./W₀[1], k)

    if RIC
        Wric = get_Wric(GC, p_weight, NN_weights)
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

        k = cut_ks(len_k, k)
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


    if RIC
        Pk = P_qso_sample(θ, kₚ, Pk_model, p, αk_model, fz_model, Qₗ, W₀k, k_start, vec(Wric))
    else
        Pk = P_qso_sample(θ, kₚ, Pk_model, p, αk_model, fz_model, Qₗ, W₀k, k_start)
    end

    return [k Pk]
end

function save_best_fit(input::String, GC::String, p_weight::Float64, p, NN_weights::Bool, RIC::Bool, output::String)
    log_density = 9
    p_weight = p_weight == -1 ? "fkp" : p_weight
    p_folder = p === p_weight ? "$p" : "$(p_weight)/$p"
    p_folder = p === nothing ? "$(p_weight)" : p_folder
    chain = joinpath(eboss_folder, input, "$(GC)GC", p_folder, "chain.npy")
    chn = npzread(chain)

    best_like = findmax(chn[:,log_density,:])
    nstep = best_like[2][1]
    nchn = best_like[2][2]

    θ = [chn[nstep,1,nchn], chn[nstep,2,nchn], chn[nstep,3,nchn], chn[nstep,4,nchn]]
    Pk = best_fit(θ, GC, p_weight, p, NN_weights, RIC)

    writedlm(output, Pk)
end

function save_best_fit(input::String, p_weight::Float64, p::Float64, NN_weights::Bool, RIC::Bool, output_N::String, output_S::String)
    log_density = 12
    p_weight = p_weight == -1 ? "fkp" : p_weight
    p_folder = p == p_weight ? "$p" : "$(p_weight)/$p"
    chain = joinpath(eboss_folder, input, "joint", p_folder, "chain.npy")
    chn = npzread(chain)

    best_like = findmax(chn[:,log_density,:])
    nstep = best_like[2][1]
    nchn = best_like[2][2]

    θ_N = [chn[nstep,1,nchn], chn[nstep,2,nchn], chn[nstep,3,nchn], chn[nstep,4,nchn]]
    Pk_N = best_fit(θ_N, "N", p_weight, p, NN_weights, RIC)
    writedlm(output_N, Pk_N)

    θ_S = [chn[nstep,1,nchn], chn[nstep,5,nchn], chn[nstep,6,nchn], chn[nstep,7,nchn]]
    Pk_S = best_fit(θ_S, "S", p_weight, p, NN_weights, RIC)
    writedlm(output_S, Pk_S)
end