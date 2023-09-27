#!/bin/bash
#SBATCH --partition squire8
#SBATCH --time=7-00:00:00
#SBATCH --array=1-4%2
#SBATCH --output=slurm-out/slurm-ezmock-power-NGC-v7-p3-%A_%a.out
#SBATCH -N 1
#SBATCH -c 3
#SBATCH --job-name=N-ezpower

#p_path = $EBOSS_DIR/eBOSS-QSO-PNG/measurements/run/window_paircount.py
((st_art = $SLURM_ARRAY_TASK_ID * 250 - 249)) #* 250 - 249
((st_op = $SLURM_ARRAY_TASK_ID * 250 + 1)) #* 250 + 1

srun -n 1 -c 3 python $EBOSS_DIR/eBOSS-QSO-PNG/measurements/run/ezmock_power.py --sample N --p 0.0 --P0_FKP 3e4 --zbins 0.8 2.2 --start $st_art --stop $st_op --step 1
