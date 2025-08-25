#!/bin/bash
#SBATCH -J all_comb16_eval
#SBATCH -t 4:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=20g
#SBATCH --array=0-59  # Sweeps ENSEMBLE_RUN
#SBATCH -o logs3/dg-%A_%a.out

i=$((${SLURM_ARRAY_TASK_ID}*1))
eval "$(conda shell.bash hook)"
conda activate cleanpy
# Define scenarios
# Loop over scenarios and ranges in parallel
for scenario in `seq 0 4`; do
    #for i in `seq 0 10 40`; do
    echo $scenario
    echo $i
    python3 create_combination_method_LLOffcompute.py 16 "$scenario" "$i"
done
