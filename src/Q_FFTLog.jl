using DelimitedFiles
using Plots
using FFTLog
using ArgParse

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

function get_poles(GC, p)
    folder = joinpath(eboss_folder, "fits", "input", "data", "window", "subsample_mean")
    name = joinpath(folder, "poles_dr16-QSO-$(GC)-p_$p.dat")
    return name
end

function get_output(GC, p)
    folder = joinpath(eboss_folder, "fits", "input", "data", "window", "Wk2")
    cup = GC == "N" ? "NGC" : "SGC"
    folder = joinpath(folder, cup, "$p")
    name = joinpath(folder, "Wk2.dat")
    return name
end

parsed_args = parse_commandline()

Q = readdlm(get_poles(parsed_args["sample"], parsed_args["p"]))
s = Q[:,1]
Q0 = Q[:,2]

HankelTest = FFTLog.HankelPlan(x=s, ν=2.01, n_extrap_low=1, n_extrap_high=1, n_pad=0)
Ell_Hankel = Array([0.5, 1.5])
prepare_Hankel!(HankelTest, Ell_Hankel)
k_Hankel = get_y(HankelTest)

Wk_Hankel = evaluate_Hankel(HankelTest, Q0 .* (s .^(1/2)))

output = get_output(parsed_args["sample"], parsed_args["p"])
writedlm(output, [k_Hankel[1,:] Wk_Hankel[1,:] ./ (k_Hankel[1,:].^(1/2))])