using ArgParse
using DelimitedFiles
using LoopVectorization
using Dierckx

eboss_folder = ENV["EBOSS_DIR"]

function parse_commandline()
    s = ArgParseSettings(description = "Define the Turing model, compute the map and the MCMC")

    choise_sample = ["N", "S"]
    choise_p = [0, 1, 1.6, 3]

    @add_arg_table s begin
        "--sample"
            help = "Galaxy cup sample, either N or S"
            arg_type = String
            range_tester = (x->x ∈ choise_sample)
            required = true
        "--p"
            help = "value of p, either 0, 1, 1.6 or 3. 0 is FKP"
            arg_type = Float64
            range_tester = (x->x ∈ choise_p)
            required = true
    end

    return parse_args(s)
end

function get_sample(sample)
    if sample == "N"
        name = "NGC"
    elseif sample == "S"
        name = "SGC"
    elseif sample == "joint"
        name = "joint"
    end
    return name
end

function get_folder_model(sample, p, NN_weights::Bool=false)
    CG = get_sample(sample)
    in_base = NN_weights ? "input_NN_weights" : "input"
    folder = joinpath(eboss_folder, "fits", in_base, "models", CG, "$p")

    return folder
end

function get_folder_data(sample, p, NN_weights::Bool=false)
    CG = get_sample(sample)
    in_base = NN_weights ? "input_NN_weights" : "input"
    folder = joinpath(eboss_folder, "fits", in_base, "data", "spectra", CG, "$p")

    return folder
end

function get_folder_window(sample, p, NN_weights::Bool=false)
    CG = get_sample(sample)
    in_base = NN_weights ? "input_NN_weights" : "input"
    folder = joinpath(eboss_folder, "fits", in_base, "data", "window", "Qkp", CG, "$p")

    return folder
end

function get_Pk_model(sample, p, NN_weights::Bool=false)
    Pk_folder = get_folder_model(sample, p, NN_weights)
    Pk_file = joinpath(Pk_folder, "Pk_model.dat")

    Pk = readdlm(Pk_file, comments=true)
    return Pk[:,2:end]
end

function get_alphak_model(sample, p, NN_weights::Bool=false)
    alphak_folder = get_folder_model(sample, p, NN_weights)
    alphak_file = joinpath(alphak_folder, "alphatildek_model.dat")

    alphak = readdlm(alphak_file, comments=true)
    return alphak[:,2:end]
end

function get_fz_model(sample, p, NN_weights::Bool=false)
    fz_folder = get_folder_model(sample, p, NN_weights)
    fz_file = joinpath(fz_folder, "Dz_fz_model.dat")

    fz = readdlm(fz_file, comments=true)
    return fz[3]
end

function get_kₚ(NN_weights::Bool=false)
    in_base = NN_weights ? "input_NN_weights" : "input"
    k_folder = joinpath(eboss_folder, "fits" , in_base, "data/window/Qkp/")
    k_file = joinpath(k_folder, "table_p.dat")

    k = readdlm(k_file, comments=true)
    return k
end

function get_Pk_data(sample, p, NN_weights::Bool=false)
    Pk_folder = get_folder_data(sample, p, NN_weights)
    Pk_file = joinpath(Pk_folder, "Pk_data.dat")

    Pk = readdlm(Pk_file, comments=true)
    return Pk[:,1], Pk[:,2]
end

function get_Σ(sample, p, NN_weights::Bool=false)
    Σ_folder = get_folder_data(sample, p, NN_weights)
    Σ_file = joinpath(Σ_folder, "covariance.dat")

    Σ = readdlm(Σ_file, comments=true)
    return Σ
end

function get_Ql(sample, p, NN_weights::Bool=false)
    Ql_folder = get_folder_window(sample, p, NN_weights)

    Qls = []
    for name in readdir(Ql_folder)
        push!(Qls, readdlm(joinpath(Ql_folder, name), comments=true))
    end

    Qₗ = Qls[1]
    for i in eachindex(Qls)[2:end]
        Qₗ = [Qₗ;;;Qls[i]]
    end

    return Qₗ
end

function get_W₀k(sample, p, NN_weights::Bool=false)
    GC = get_sample(sample)
    in_base = NN_weights ? "input_NN_weights" : "input"
    folder = joinpath(eboss_folder, "fits", in_base, "data", "window", "Wk2", GC, "$p")
    name = joinpath(folder, "Wk2.dat")

    W₀k = readdlm(name)
    return W₀k[:,1], W₀k[:,2]
end

function get_Wric(sample, p, NN_weights::Bool=false)
    GC = get_sample(sample)
    in_base = NN_weights ? "input_NN_weights" : "input"
    folder = joinpath(eboss_folder, "fits", in_base, "data", "window", "Wkric", GC, "$p")
    name = joinpath(folder, "Wkric.dat")

    Wkric = readdlm(name, comments=true)
    return Wkric
end

function check_k_dimension(k_data, Σ)
    len_data = length(k_data)
    len_cov = length(Σ[:,1])

    return len_data == len_cov
end

function cut_ks(len::Int, to_cut::Vector)
    return to_cut[1:len]
end

function cut_ks(len::Int, to_cut::Matrix)
    return to_cut[1:len,:]
end

function cut_ks(len::Int, to_cut::Array)
    return to_cut[1:len,:,:]
end

function front_cut_ks(start::Int, to_cut::Vector)
    return to_cut[start:end]
end

function get_output_folder(sample, p, NN_weights::Bool=false)
    CG = get_sample(sample)
    out_base = NN_weights ? "output_NN_weights" : "output"
    folder = joinpath(eboss_folder, "fits", out_base, CG, "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "fits"))
        if isdir(out_base)
            cd(out_base)
        else
            mkdir(out_base)
            cd(out_base)
        end
        if isdir(CG)
            cd(CG)
        else
            mkdir(CG)
            cd(CG)
        end
        if isdir("$p")
            cd("$p")
        else
            mkdir("$p")
            cd("$p")
        end
        cd(eboss_folder)
        return folder
    end
end

function get_output_folder(sample, p, p_weight, NN_weights::Bool=false)
    CG = get_sample(sample)
    out_base = NN_weights ? "output_NN_weights" : "output"
    folder = joinpath(eboss_folder, "fits", out_base, CG,  "$p_weight", "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "fits"))
        if isdir(out_base)
            cd(out_base)
        else
            mkdir(out_base)
            cd(out_base)
        end
        if isdir(CG)
            cd(CG)
        else
            mkdir(CG)
            cd(CG)
        end
        if isdir("$p_weight")
            cd("$p_weight")
        else
            mkdir("$p_weight")
            cd("$p_weight")
        end
        if isdir("$p")
            cd("$p")
        else
            mkdir("$p")
            cd("$p")
        end
        cd(eboss_folder)
        return folder
    end
end

function _mygemmavx(A, B, C)
    Dm = zero(eltype(C))
    @turbo for n ∈ axes(A,1)
        Dm += A[n] * B[n] * C[n]
    end
    return Dm
end

function interpolate_fk(k, fk, new_k)
    fk_interp = Spline1D(k, fk)
    return fk_interp(new_k)
end