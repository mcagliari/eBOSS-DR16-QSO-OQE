#IMPORTANT NOTE: this code just simplify running the sbatch command, but the real parameters of the runs must be changed in the *-distributed*.jl file
using ArgParse

eboss_folder = ENV["EBOSS_DIR"]

s = ArgParseSettings(description = "Run distributed chains for B or P+B on slurm")

choise_sample = ["N", "S", "j"]
choise_p = [1.0, 1.6]
choise_w = ["f", "o"]
choise_node = ["node1", "node2", "node3", "node4", "node5", "node6", "node7", "node8"] 
choise_analysis = ["B", "PB", "PcB"]

@add_arg_table s begin
    "--sample"
        help = "Galaxy cup sample, either N or S or joint (j)"
        arg_type = String
        range_tester = (x->x ∈ choise_sample)
        required = true
    "--p"
        help = "value of p, either 1, 1.6"
        arg_type = Float64
        range_tester = (x->x ∈ choise_p)
        required = true
    "--w"
        help = "type of weights FKP (f) or Optimal (o)"
        arg_type = String
        range_tester = (x->x ∈ choise_w)
        required = true
    "--node"
        help = "Node to run on"
        arg_type = String
        range_tester = (x->x ∈ choise_node)
        required = true
    "--analysis"
        help = "type of analysis B or P+B or compression"
        arg_type = String
        range_tester = (x->x ∈ choise_analysis)
        required = true
end

parsed_args = parse_args(s)

GC = parsed_args["sample"]
p = parsed_args["p"]
w = parsed_args["w"]
node = parsed_args["node"]
analysis = parsed_args["analysis"]

slurm_name = analysis * GC * "-" * w
slurm_name *= p == 1.0 ? "1" : "6"


if GC == "j"
    if analysis == "B"
        script_name = joinpath(eboss_folder, "bis-src", "fnl-bispectrum-distributed-joint.jl")
    elseif analysis == "PB"
        script_name = joinpath(eboss_folder, "bis-src", "fnl-P_plus_B-distributed-joint.jl")
    elseif analysis == "PcB"
        script_name = joinpath(eboss_folder, "bis-src", "fnl-P_plus_cB-distributed-joint.jl")
    end
else
    if analysis == "B"
        script_name = joinpath(eboss_folder, "bis-src", "fnl-bispectrum-distributed.jl")
    elseif analysis == "PB"
        script_name = joinpath(eboss_folder, "bis-src", "fnl-P_plus_B-distributed.jl")
    elseif analysis == "PcB"
        script_name = joinpath(eboss_folder, "bis-src", "fnl-P_plus_cB-distributed.jl")
    end
end

cmd = `sbatch -psquire8 --time=5-00:00:00 --job-name=$slurm_name -w$node -N1 -c12 --wrap="julia $script_name"`

run(cmd)
