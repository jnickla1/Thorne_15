#!/bin/bash
#SBATCH -J histens
#SBATCH -t 2:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=20g
#SBATCH --array=0-19  # Sweeps ENSEMBLE_RUN
#SBATCH -o logs3/histens3e-%A_%a.out

#i=$((${SLURM_ARRAY_TASK_ID}*5))
#eval "$(conda shell.bash hook)"
#TO RUN
#export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
#seq 0 5 100 | xargs -P 16 -n 1 ./run_histens_local.sh

i=$1
mkdir -p logs_histens
exec > "logs_histens/run_${i}.log" 2>&1   # one log per i so output doesn't interleave

# conda activate doesn't work in a non-interactive script until conda is initialised
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cleanPy 


# Define scenarios:q


scenarios=(
    "histens_satcal_real"
    "histens_rpi_real"
    "histens_satcal_anthro"
    "histens_rpi_anthro"
)

#for i in $(seq 0 5 100); do
# Loop over scenarios and ranges in parallel



for scenario in "${scenarios[@]}"; do
    echo $scenario
    echo $i
    python3 hist_heads_tails_evaluation_script.py "$scenario" "$i"
done

#done
