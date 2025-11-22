#!/bin/bash
#SBATCH -J histens
#SBATCH -t 2:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=20g
#SBATCH --array=0-19  # Sweeps ENSEMBLE_RUN
#SBATCH -o logs3/histens3e-%A_%a.out

i=$((${SLURM_ARRAY_TASK_ID}*5))
eval "$(conda shell.bash hook)"
conda activate cleanpy
# Define scenariosi


scenarios=(
    "histens_satcal_real"
    "histens_rpi_real"
    "histens_satcal_anthro"
    "histens_rpi_anthro"
)

# Loop over scenarios and ranges in parallel
for scenario in "${scenarios[@]}"; do
    echo $scenario
    echo $i
    python3 hist_heads_tails_evaluation_script.py "$scenario" "$i"
done
