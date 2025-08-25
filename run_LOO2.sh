#!/bin/bash
#SBATCH -J all_comp2_eval
#SBATCH -t 2:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=20g
#SBATCH --array=25-29  # Sweeps ENSEMBLE_RUN
#SBATCH -o logs3/ee2e-%A_%a.out

i=$((${SLURM_ARRAY_TASK_ID}*2))
eval "$(conda shell.bash hook)"
conda activate cleanpy
# Define scenarios
# Loop over scenarios and ranges in parallel
for scenario in `seq 0 4`; do
    #for i in `seq 0 10 40`; do
    echo $scenario
    echo $i
    python3 fut_2_script.py "$scenario" "$i"
done
