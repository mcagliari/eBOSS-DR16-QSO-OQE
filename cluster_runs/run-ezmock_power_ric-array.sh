#!/bin/bash
#SBATCH --partition squire8
#SBATCH --time=7-00:00:00
#SBATCH --array=1-4%2
##SBATCH --array=1-8%8
#SBATCH --output=slurm-out/slurm-ezmock-power-SGC-v7-p0.0-%A_%a.out
#SBATCH -N 1
#SBATCH -c 3
#SBATCH --job-name=S-shuf
#SBATCH --mem=0

module load intel impi

#p_path = $EBOSS_DIR/eBOSS-QSO-PNG/measurements/run/window_paircount.py
((st_art = $SLURM_ARRAY_TASK_ID * 250 - 249)) #* 250 - 249
((st_op = $SLURM_ARRAY_TASK_ID * 250 + 1)) #* 250 + 1

#((st_art = $SLURM_ARRAY_TASK_ID * 125 - 124)) #* 250 - 249
#((st_op = $SLURM_ARRAY_TASK_ID * 125 + 1)) #* 250 + 1


srun -n 1 -c 3 python $EBOSS_DIR/src/ezmock_power_ric_rshuffle.py --sample S --p 0.0 --P0_FKP 3e4 --zbins 0.8 2.2 --start $st_art --stop $st_op --step 1

###srun -n 1 -c 3 python $EBOSS_DIR/src/ezmock_power_ric_complete.py --sample S --p 0.0 --P0_FKP 3e4 --zbins 0.8 2.2 --start $st_art --stop $st_op --step 1