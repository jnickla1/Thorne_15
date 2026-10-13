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
# Temporary coordination for this correction's parallel dispatchers.
if [[ -f logs_histens/.kalman_refresh_active ]]; then
    refresh_root="logs_histens/.kalman_refresh"
    mkdir -p "$refresh_root"
    [[ -f "$refresh_root/$i.done" ]] && exit 0
    mkdir "$refresh_root/$i.lock" 2>/dev/null || exit 0
    trap 'status=$?; if [[ $status == 0 ]]; then touch "$refresh_root/$i.done"; fi; rmdir "$refresh_root/$i.lock"; exit "$status"' EXIT
fi
set -e
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
    python3 evaluate_heads_tails_ensemble.py "$scenario" "$i"
done

#done
