import os
import argparse
import eboss_qso.measurements as eboss
from eboss_qso import EBOSS_FITS

parser = argparse.ArgumentParser(description="Compute the model P(k), T(k), D, f and alphatilde at a given z_eff")
    
h = 'effective redshift'
parser.add_argument('--zeff', type = float, help = h, required = True)

h = 'the sample, either North or South'
parser.add_argument('--sample', type = str, choices = ['N', 'S'], help = h, required = True)

h = 'the value of p to use'
parser.add_argument('--p', type = float, help = h, choices=[0., 1., 1.6, 3.0], required = True)

h = 'path to D(z) and f(z) file'
Df_file = os.path.join(EBOSS_FITS, 'input', 'models', 'Dz_fz_Planck18.txt')
parser.add_argument('--Dffile', type = str, help = h, default = Df_file)

h = 'flag for raw Tk file'
parser.add_argument('--raw', action = 'store_true')

args = parser.parse_args()

slurm_dir = os.path.join(eboss.data_dir, 'doraemon-runs', 'slurm-out')
slurm_name = f"slurm-model-{args.sample}GC-p{args.p}.out"
slurm_out = os.path.join(slurm_dir, slurm_name)
######
p_path = os.path.join(eboss.data_dir, 'src/', 'cosmo_model.py') 

if args.raw:
    cmd = f'''sbatch -psquire8 --output={slurm_out} --wrap="python {p_path} --zeff {args.zeff} --sample {args.sample} --p {args.p} --Dffile {args.Dffile} --raw"'''
else:
    cmd = f'''sbatch -psquire8 --output={slurm_out} --wrap="python {p_path} --zeff {args.zeff} --sample {args.sample} --p {args.p} --Dffile {args.Dffile}"'''  
    
print(cmd)

os.system(cmd)