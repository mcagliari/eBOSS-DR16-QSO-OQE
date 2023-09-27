#!/bin/bash
#SBATCH --partition squire8
#SBATCH --time=28-00:00:00
#SBATCH --array=1-5%5
#SBATCH --output=slurm-out/slurm-window-paircount-SGC-vdr16-p0.0-%A_%a.out
#SBATCH -w node3
#SBATCH --job-name=S-wind

###p_path = $EBOSS_DIR/eBOSS-QSO-PNG/measurements/run/window_paircount.py

python $EBOSS_DIR/eBOSS-QSO-PNG/measurements/run/window_paircount.py --sample S --version dr16 --zmin 0.8 --zmax 2.2 --subsample 1440000 --p 0.0 --ell 0 --nsub $SLURM_ARRAY_TASK_ID --NN-weights

###python $EBOSS_DIR/eBOSS-QSO-PNG/measurements/run/window_paircount.py --sample N --version dr16 --zmin 0.8 --zmax 2.2 --subsample 2200000 --p 0.0 --ell 0 --nsub $SLURM_ARRAY_TASK_ID --NN-weights
