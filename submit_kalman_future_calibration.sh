#!/bin/bash
#SBATCH -J kf-future-calibration
#SBATCH -t 12:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=12g
#SBATCH --array=0-269
#SBATCH -o logs_kf_calibration/%A_%a.out
set -euo pipefail
module load matlab/R2019a-rjyk3ws
module load r
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cleanpy
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
idx=$SLURM_ARRAY_TASK_ID
if (( idx < 150 )); then
 scenarios=(fut_ESM1-2-LR_SSP126_constVolc fut_ESM1-2-LR_SSP245_constVolc fut_ESM1-2-LR_SSP370_constVolc)
 scenario=${scenarios[$((idx / 50))]}
 member=$((idx % 50))
else
 scenarios=(fut_NorESM_RCP45_Volc fut_NorESM_RCP45_VolcConst)
 scenario=${scenarios[$(((idx-150) / 60))]}
 member=$(((idx-150) % 60))
fi
python3 calibrate_method_uncertainty.py "$scenario" "$member"
