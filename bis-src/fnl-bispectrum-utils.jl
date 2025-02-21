using DelimitedFiles
using Zygote
using ChainRulesCore
using LinearAlgebra

eboss_folder = ENV["EBOSS_DIR"]

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

function get_folder_model(sample, rebin::String="") #z_eff between liner and NN is basically the same, I use only one model
    GC = get_sample(sample)
    in_base = "input" * rebin
    folder = joinpath(eboss_folder, "bis-fits", in_base, "models", GC)

    return folder
end

function get_folder_data(sample, rebin::String="", NN_weights::Bool=false)
    GC = get_sample(sample)
    in_base = NN_weights ? "input_NN_weights" * rebin : "input" * rebin
    folder = joinpath(eboss_folder, "bis-fits", in_base, "data", "bispectra", GC)

    return folder
end

function get_Pk_model(sample, convolved::Bool=true, rebin::String="")
    Pk_folder = get_folder_model(sample, rebin)
    Pk_name = convolved ? "Pk_window.dat" : "Pk_model.dat"
    Pk_file = joinpath(Pk_folder, Pk_name)

    Pk = readdlm(Pk_file, comments=true)
    return Pk
end

function get_Minv_model(sample, rebin::String="")
    Minv_folder = get_folder_model(sample, rebin)
    Minv_file = joinpath(Minv_folder, "Minv_model.dat")

    Minv = readdlm(Minv_file, comments=true)
    return Minv
end

function get_fz_model(sample, rebin::String="")
    fz_folder = get_folder_model(sample, rebin)
    fz_file = joinpath(fz_folder, "Dz_fz_model.dat")

    fz = readdlm(fz_file, comments=true)
    return fz[3]
end

function get_kB(rebin::String="")
    in_base = "input" * rebin
    k_folder = joinpath(eboss_folder, "bis-fits" , in_base, "models")
    k_file = joinpath(k_folder, "k_eff-bisp.txt")

    k = readdlm(k_file, comments=true)
    return k
end

function get_Bk_data(sample, B_ind::Int64=4, rebin::String="", NN_weights::Bool=false)
    Bk_folder = get_folder_data(sample, rebin, NN_weights)
    Bk_file = joinpath(Bk_folder, "Bk_data.dat") #check name

    Bk = readdlm(Bk_file, comments=true)
    return Bk[:,1], Bk[:,2], Bk[:,3], Bk[:,B_ind] #k1-eff, k2-eff, k3-eff, B0 (if B_ind=4) or B0-Sn (if B_ind=5)
end

function get_S0(sample)
    GC = get_sample(sample)
    S0_folder = joinpath(eboss_folder, "measurements", "bispectra", "data")
    S0_name = "S0_" * GC * "_data.dat" #pensare al caso mock Polynomials
    S0_file = joinpath(S0_folder, S0_name)

    S0 = readdlm(S0_file, comments=true)
    return S0[1]
end

function get_S0B(sample)
    GC = get_sample(sample)
    S0_folder = joinpath(eboss_folder, "measurements", "bispectra", "data")
    S0_name = "S0B_" * GC * "_data.dat" #pensare al caso mock Polynomials
    S0_file = joinpath(S0_folder, S0_name)

    S0 = readdlm(S0_file, comments=true)
    return S0[1,end]
end

function get_Σ(sample, PpB::Bool=false, pP::Number=-1, SN::Bool=false, rebin::String="") #the covariance has linear weights
    Σ_folder = get_folder_data(sample, rebin)
    GC = get_sample(sample)
    SN_name = SN ? "SNsub_" : ""

    if PpB
        if pP == -1
            Σ_name = "covariance_" * SN_name * GC * "_PB.dat"
        else
            Σ_name = "covariance_" * SN_name * GC * "_PB_opt$pP.dat"
        end
    else
        Σ_name = "covariance_" * SN_name * GC * "_B.dat"
    end

    Σ_file = joinpath(Σ_folder, Σ_name)

    Σ = readdlm(Σ_file, comments=true)
    return Σ
end

function get_Wric(sample, SN::Bool=false, rebin::String="") #all the mocks have linear weights
    GC = get_sample(sample)
    in_base = "input" * rebin
    folder = joinpath(eboss_folder, "bis-fits", in_base, "data", "window", "Wkric", GC)
    SN_name = SN ? "_SNsub" : ""
    name = joinpath(folder, "Wkric" * SN_name * ".dat")

    Wkric = readdlm(name, comments=true)
    return Wkric[:,4]
end

function get_PaB_dimensions(sample, rebin::String="")
    B_folder = get_folder_data(sample, rebin)
    PaB_name = joinpath(B_folder, "cut_dimensions.dat")
    
    PaB = readdlm(PaB_name, Int, comments=true) #should be P_dim, B_dim after the cut necessary to build the covariance
    return PaB
end

function mask_ks(mask::BitVector, to_mask::Vector)
    return to_mask[mask]
end

function mask_ks(mask::BitVector, to_mask::Matrix)
    return to_mask[mask,:]
end

function mask_Σ(mask::BitVector, to_mask::Matrix)
    return to_mask[mask,mask]
end

function get_rotation(sample, len_rotation::Int64, cPB::Bool=false)
    GC = get_sample(sample)
    folder = joinpath(eboss_folder, "measurements/bispectra/compression/$GC/")
    file = cPB ? folder * "compression_matrix_PB_GaussCov_$GC.dat" : folder * "compression_matrix_Bonly_GaussCov_$GC.dat"

    rotation = readdlm(file, comments=true)
    return rotation[1:len_rotation,:]
end

function get_Bmocks(sample, SN::Bool)
    GC = get_sample(sample)
    folder = joinpath(eboss_folder, "bis-fits/input/data/prova/$GC/realistic/")

    B_ind = SN ? 5 : 4
    name = folder * "Bispectrum_$(GC)_EZmock_realistic_1.dat"
    Bks = readdlm(name, comments=true)[:,B_ind]

    n = 1000
    for i in 2:n 
        name = folder * "Bispectrum_$(GC)_EZmock_realistic_$i.dat"
        Bk = readdlm(name, comments=true)[:,B_ind]
        Bks = hcat(Bks, Bk)
    end

    return Bks
end

function get_Bmocks(sample, len_compression::Int64)
    GC = get_sample(sample)
    folder = joinpath(eboss_folder, "/measurements/bispectra/compression/$GC/realistic/")

    name = folder * "Bispectrum_$(GC)_EZmock_realistic_BispCompressed_1.dat"
    Bks = readdlm(name, comments=true)[1:len_compression]

    n = 1000
    for i in 2:n 
        name = folder * "Bispectrum_$(GC)_EZmock_realistic_BispCompressed_$i.dat"
        Bk = readdlm(name, comments=true)[1:len_compression]
        Bks = hcat(Bks, Bk)
    end

    return Bks
end

function get_rotate_SN(sample, R::Matrix)
    GC = get_sample(sample)
    folder = joinpath(eboss_folder, "bis-fits/input/data/prova/$GC/realistic/")

    name = folder * "Bispectrum_$(GC)_EZmock_realistic_1.dat"
    SN = readdlm(name, comments=true)[:,6] #speriamo siano tutti della dimensione giusta
    SNrs = R * SN

    n = 1000
    for i in 2:n 
        name = folder * "Bispectrum_$(GC)_EZmock_realistic_$i.dat"
        SN = readdlm(name, comments=true)[:,6]
        SNr = R * SN
        SNrs = hcat(SNrs, SNr)
    end

    return SNrs
end

function get_output_folderB(sample, p, PpB::Bool=false, rebin::String="", NN_weights::Bool=false)
    GC = get_sample(sample)
    out_base = PpB ? "output-PplusB" : "output-B"
    out_base *= NN_weights ? "_NN_weights" : ""
    out_base *= rebin
    folder = joinpath(eboss_folder, "bis-fits", out_base, GC, "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "bis-fits"))
        if isdir(out_base)
            cd(out_base)
        else
            mkdir(out_base)
            cd(out_base)
        end
        if isdir(GC)
            cd(GC)
        else
            mkdir(GC)
            cd(GC)
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

function get_output_folderB(sample, p, p_weight, PpB::Bool=false, rebin::String="", NN_weights::Bool=false)
    CG = get_sample(sample)
    out_base = PpB ? "output-PplusB" : "output-B"
    out_base *= NN_weights ? "_NN_weights" : ""
    out_base *= rebin
    
    p_weight = p_weight == -1 ? "fkp" : p_weight
    folder = joinpath(eboss_folder, "bis-fits", out_base, CG,  "$p_weight", "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "bis-fits"))
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

function get_output_folderB(sample, p, name::String, PpB::Bool=false, rebin::String="", NN_weights::Bool=false)
    GC = get_sample(sample)
    out_base = PpB ? "output-PplusB" : "output-B"
    out_base *= NN_weights ? "_NN_weights" : ""
    out_base *= rebin
    
    folder = joinpath(eboss_folder, "bis-fits", out_base, GC, name, "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "bis-fits"))
        if isdir(out_base)
            cd(out_base)
        else
            mkdir(out_base)
            cd(out_base)
        end
        if isdir(GC)
            cd(GC)
        else
            mkdir(GC)
            cd(GC)
        end
        if isdir(name)
            cd(name)
        else
            mkdir(name)
            cd(name)
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

function get_output_folderB(sample, p, p_weight, name::String, PpB::Bool=false, rebin::String="", NN_weights::Bool=false)
    CG = get_sample(sample)
    out_base = PpB ? "output-PplusB" : "output-B"
    out_base *= NN_weights ? "_NN_weights" : ""
    out_base *= rebin
    
    p_weight = p_weight == -1 ? "fkp" : p_weight
    folder = joinpath(eboss_folder, "bis-fits", out_base, CG, name,  "$p_weight", "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "bis-fits"))
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
        if isdir(name)
            cd(name)
        else
            mkdir(name)
            cd(name)
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

function get_output_foldercB(sample, p, len_compression, name::String, PpB::Bool=true, NN_weights::Bool=false)
    GC = get_sample(sample)
    out_base = PpB ? "output-PplusB-compression" : "output-B-compression"
    out_base *= NN_weights ? "_NN_weights" : ""
    
    folder = joinpath(eboss_folder, "bis-fits", out_base, "compression-$(len_compression)", GC, name, "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "bis-fits"))
        if isdir(out_base)
            cd(out_base)
        else
            mkdir(out_base)
            cd(out_base)
        end
        if isdir("compression-$(len_compression)")
            cd("compression-$(len_compression)")
        else
            mkdir("compression-$(len_compression)")
            cd("compression-$(len_compression)")
        end
        if isdir(GC)
            cd(GC)
        else
            mkdir(GC)
            cd(GC)
        end
        if isdir(name)
            cd(name)
        else
            mkdir(name)
            cd(name)
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

function get_output_foldercB(sample, p, p_weight::Float64, len_compression::Integer, name::String, PpB::Bool=true, NN_weights::Bool=false)
    CG = get_sample(sample)
    out_base = PpB ? "output-PplusB-compression" : "output-B-compression"
    out_base *= NN_weights ? "_NN_weights" : ""
    
    p_weight = p_weight == -1 ? "fkp" : p_weight
    folder = joinpath(eboss_folder, "bis-fits", out_base, "compression-$(len_compression)", CG, name,  "$p_weight", "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "bis-fits"))
        if isdir(out_base)
            cd(out_base)
        else
            mkdir(out_base)
            cd(out_base)
        end
        if isdir("compression-$(len_compression)")
            cd("compression-$(len_compression)")
        else
            mkdir("compression-$(len_compression)")
            cd("compression-$(len_compression)")
        end
        if isdir(CG)
            cd(CG)
        else
            mkdir(CG)
            cd(CG)
        end
        if isdir(name)
            cd(name)
        else
            mkdir(name)
            cd(name)
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

function get_output_foldercB(sample, p, len_compressionN::Integer, len_compressionS::Integer, name::String, PpB::Bool=true, NN_weights::Bool=false)
    GC = get_sample(sample)
    out_base = PpB ? "output-PplusB-compression" : "output-B-compression"
    out_base *= NN_weights ? "_NN_weights" : ""
    
    folder = joinpath(eboss_folder, "bis-fits", out_base, "compression-N$(len_compressionN)-S$(len_compressionS)", GC, name, "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "bis-fits"))
        if isdir(out_base)
            cd(out_base)
        else
            mkdir(out_base)
            cd(out_base)
        end
        if isdir("compression-N$(len_compressionN)-S$(len_compressionS)")
            cd("compression-N$(len_compressionN)-S$(len_compressionS)")
        else
            mkdir("compression-N$(len_compressionN)-S$(len_compressionS)")
            cd("compression-N$(len_compressionN)-S$(len_compressionS)")
        end
        if isdir(GC)
            cd(GC)
        else
            mkdir(GC)
            cd(GC)
        end
        if isdir(name)
            cd(name)
        else
            mkdir(name)
            cd(name)
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

function get_output_foldercB(sample, p, p_weight, len_compressionN, len_compressionS, name::String, PpB::Bool=true, NN_weights::Bool=false)
    CG = get_sample(sample)
    out_base = PpB ? "output-PplusB-compression" : "output-B-compression"
    out_base *= NN_weights ? "_NN_weights" : ""
    
    p_weight = p_weight == -1 ? "fkp" : p_weight
    folder = joinpath(eboss_folder, "bis-fits", out_base, "compression-N$(len_compressionN)-S$(len_compressionS)", CG, name,  "$p_weight", "$p")

    out_base = join([out_base, "/"])
    check = isdir(folder)
    if check
        return folder
    else
        cd(joinpath(eboss_folder, "bis-fits"))
        if isdir(out_base)
            cd(out_base)
        else
            mkdir(out_base)
            cd(out_base)
        end
        if isdir("compression-N$(len_compressionN)-S$(len_compressionS)")
            cd("compression-N$(len_compressionN)-S$(len_compressionS)")
        else
            mkdir("compression-N$(len_compressionN)-S$(len_compressionS)")
            cd("compression-N$(len_compressionN)-S$(len_compressionS)")
        end
        if isdir(CG)
            cd(CG)
        else
            mkdir(CG)
            cd(CG)
        end
        if isdir(name)
            cd(name)
        else
            mkdir(name)
            cd(name)
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