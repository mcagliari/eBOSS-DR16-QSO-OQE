#conda activate fnl
import os
import argparse
import eboss_qso.measurements as eboss

parser = argparse.ArgumentParser(description="Run with sbatch my-scr/window-poles.py")

h = 'the sample, either North or South'
parser.add_argument('--sample', type = str, choices = ['N', 'S'], help = h, required = True)

h = 'the version to load'
parser.add_argument('--version', type = str, choices = eboss.DATA_VERSIONS, help = h, default = 'dr16')

h = 'the minimum redshift to include'
parser.add_argument('--zmin', type=float, help=h, default=0.8)

h = 'the maximum redshift to include'
parser.add_argument('--zmax', type=float, help=h, default=2.2)

h = 'the desired collective size to subsample to'
parser.add_argument('--subsample', type=float, help=h)

h = 'the value of p to use'
parser.add_argument('--p', type=float, help=h, choices=[-1., 0., 1., 1.6, 3.0], required=True)

h = 'iteration number of the DD paircount file'
parser.add_argument('--subnumber', type = int, help = h, default = None)

h = 'identifying name on the DD paircount file if not hash named'
parser.add_argument('--name', type = str, help = h, default = None)

h = 'whether to use the catalogs with NN systot weights, to be used only for dr16'
parser.add_argument('--NN-weights', action='store_true', help=h)


args = parser.parse_args()


slurm_dir = os.path.join(eboss.data_dir, 'doraemon-runs', 'slurm-out')
slurm_name = f"slurm-window-poles-{args.sample}GC-v{args.version}-sub{args.subsample}-p{args.p}.out"
slurm_out = os.path.join(slurm_dir, slurm_name)
######
p_path = os.path.join(eboss.data_dir, 'src/', 'window-poles.py') 

if args.NN_weights and args.version == 'dr16':
    cmd = f'''sbatch -psquire8 --time=14-00:00:00 --output={slurm_out} --wrap="python {p_path} --sample {args.sample} --version {args.version} --zmin {args.zmin} --zmax {args.zmax} --subsample {args.subsample} --p {args.p} --subnumber {args.subnumber} --NN-weights"'''
else:
    cmd = f'''sbatch -psquire8 --time=14-00:00:00 --output={slurm_out} --wrap="python {p_path} --sample {args.sample} --version {args.version} --zmin {args.zmin} --zmax {args.zmax} --subsample {args.subsample} --p {args.p} --subnumber {args.subnumber}"'''
    
print(cmd)

os.system(cmd)
