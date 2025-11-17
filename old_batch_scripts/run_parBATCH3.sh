#!/bin/bash
#SBATCH -J all_crossing_eval
#SBATCH -t 8:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=5g
#SBATCH --array=0,1,2  #0-30  # Sweeps ENSEMBLE_RUN
#SBATCH -o logs4/rab8-%A_%a.out

i=$((${SLURM_ARRAY_TASK_ID}*2))
eval "$(conda shell.bash hook)"
conda activate cleanpy
# Define scenarios
scenarios=(
    "fut_ESM1-2-LR_SSP126_constVolc"
    "fut_ESM1-2-LR_SSP245_constVolc"
    "fut_ESM1-2-LR_SSP370_constVolc"
)

# Loop over scenarios and ranges in parallel
for scenario in "${scenarios[@]}"; do
    #for i in `seq 0 10 40`; do
    if [ "$i" -lt 50 ]; then
        python3 fut_evaluation_script.py "$scenario" "$i" &
    fi
done

scenarios=(
    "fut_NorESM_RCP45_Volc"
    "fut_NorESM_RCP45_VolcConst"
)

# Loop over scenarios and ranges in parallel
for scenario in "${scenarios[@]}"; do
    #for i in `seq 0 10 50`; do
        python3 fut_evaluation_script.py "$scenario" "$i" &
    #done
done


# Wait for all background processes to finish
wait

echo "All processes completed."
