using ArgParse

eboss_folder = ENV["EBOSS_DIR"]

s = ArgParseSettings(description = "Run fnl-distributed.jl on slurm")

choise_sample = ["N", "S", "joint"]
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

parsed_args = parse_args(s)

GC = parsed_args["sample"]
p = parsed_args["p"]

slurm_out = joinpath(eboss_folder, "doraemon-runs", "slurm-out", "slurm-distributed-ezmocks-$GC-p$p.out")
slurm_name = "$GC-p$p"
if GC == "joint"
    script_name = joinpath(eboss_folder, "src", "fnl-distributed-joint.jl")
else
    script_name = joinpath(eboss_folder, "src", "fnl-distributed.jl")
end

cmd = `sbatch -psquire8 --time=4-00:00:00 --output=$slurm_out --job-name=$slurm_name -wnode6 -N1 -c8 --wrap="julia $script_name --sample $GC --p $p"`

run(cmd)
