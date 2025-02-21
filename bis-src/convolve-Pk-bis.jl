using DelimitedFiles
using ArgParse
using Tullio
using Dierckx

include("fnl-bispectrum-utils.jl")

function parse_commandline()
    s = ArgParseSettings(description = "Convolve the linear matter power spectrum with the window monopole for the bispectrum model")

    choise_sample = ["N", "S"]

    @add_arg_table s begin
        "--sample"
            help = "Galaxy cup sample, either N or S"
            arg_type = String
            range_tester = (x->x ∈ choise_sample)
            required = true
        "--rebin"
           help = "rebin folder"
           arg_type = String
           default = ""
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
GC = get_sample(parsed_args["sample"])
in_base = "input" * parsed_args["rebin"]

Pkp_file = joinpath(get_folder_model(parsed_args["sample"], parsed_args["rebin"]), "Pkp_model.dat")
Pkp = readdlm(Pkp_file, comments=true)
Qk0_file = joinpath(eboss_folder, "fits", "input", "data", "window", "Qkp", GC, "fkp", "Qkp0_" * GC * "_FKP.dat")
Qkp0 = readdlm(Qk0_file, comments=true)
kp = readdlm(joinpath(eboss_folder, "fits", "input", "data", "window", "Qkp", "table_p.dat"), comments=true)


l = length(Qkp0[:,1])


k = kp[1:l,1] .* 10. ./ 2
#Pkp = Pkp[:,2:end]

ks = get_kB(parsed_args["rebin"])

println(size(Qkp0), size(Pkp))

@tullio Pk[i] := Qkp0[i,k] * Pkp[i,k]

Pk_interp = Spline1D(k, Pk)
Pk1 = Pk_interp(ks[:,1])
Pk2 = Pk_interp(ks[:,2])
Pk3 = Pk_interp(ks[:,3])

Pks = [Pk1 Pk2 Pk3]

model_folder = joinpath(eboss_folder, "bis-fits", in_base, "models", GC)
writedlm(joinpath(model_folder, "Pk_window.dat"), Pks)

println(joinpath(model_folder, "Pk_window.dat", "ready!"))