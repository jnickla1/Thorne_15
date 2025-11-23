#!/bin/bash
#SBATCH -J all_comb30_eval
#SBATCH -t 4:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=20g
#SBATCH --array=26  # Sweeps ENSEMBLE_RUN
#SBATCH -o logs3/dee-%A_%a.out

i=$((${SLURM_ARRAY_TASK_ID}*1))
eval "$(conda shell.bash hook)"
conda activate cleanpy
# Define scenarios
# Loop over scenarios and ranges in parallel
i=4
    #for i in `seq 0 10 40`; do
echo $scenario
echo $i
python3 create_combination_method_LLOffcompute.py 30 "$scenario" "$i"
