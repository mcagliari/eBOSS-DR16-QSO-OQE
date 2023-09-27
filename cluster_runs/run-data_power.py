#conda activate fnl
import os
import argparse
import eboss_qso.measurements as eboss

parser = argparse.ArgumentParser(description="Run with sbatch data_power.py")

h = 'the sample, either North or South'
parser.add_argument('--sample', type = str, choices = ['N', 'S'], help = h, required = True)

h = 'the version to load'
parser.add_argument('--version', type = str, choices = eboss.DATA_VERSIONS, help = h, default = 'dr16')

h = 'the redshift bins to compute in the form of ZMIN, ZMAX'
parser.add_argument('--zbins', type = float, nargs = '+', help = h, required = True)

h = 'the value of p to use'
parser.add_argument('--p', type = float, help = h, choices=[0., 1., 1.6, 3.0], required = True)

h = 'the P0 FKP version to use'
parser.add_argument('--P0_FKP', type = float, default = 3e4, help = h)

h = 'whether to use focal plane weights, to be used only for dr14'
parser.add_argument('--focal-weights', action = 'store_true', help = h)

h = 'whether to use the catalogs with NN systot weights, to be used only for dr16'
parser.add_argument('--NN-weights', action='store_true', help=h)

h = 'number of nodes'
parser.add_argument('--nodes', type = int, default = 1)

args = parser.parse_args()

zmin, zmax = args.zbins

slurm_dir = os.path.join(eboss.data_dir, 'doraemon-runs', 'slurm-out')
slurm_name = f"slurm-data-power-{args.sample}GC-v{args.version}-p{args.p}.out"
slurm_out = os.path.join(slurm_dir, slurm_name)

p_path = os.path.join(eboss.data_dir, 'eBOSS-QSO-PNG/measurements/run', 'data_power.py') 

mpi = ''
if args.nodes > 1:
    mpi = f'srun -n {args.nodes}'


if args.focal_weights and args.version == 'dr14':
    cmd = f'''sbatch -psquire8 --mem=0 --output={slurm_out} --wrap="python {p_path} --sample {args.sample} --version {args.version} --zbins {zmin},{zmax} --p {args.p} --P0_FKP {args.P0_FKP} --focal-weights"'''
elif args.NN_weights and args.version == 'dr16':
    cmd = f'''sbatch -psquire8 --mem=0 --output={slurm_out} -N{args.nodes} --wrap="{mpi} python {p_path} --sample {args.sample} --version {args.version} --zbins {zmin},{zmax} --p {args.p} --P0_FKP {args.P0_FKP} --NN-weights"'''
else:
    cmd = f'''sbatch -psquire8 --mem=0 --output={slurm_out} -N{args.nodes} --wrap="{mpi} python {p_path} --sample {args.sample} --version {args.version} --zbins {zmin},{zmax} --p {args.p} --P0_FKP {args.P0_FKP}"'''
    
print(cmd)

os.system(cmd)
 
