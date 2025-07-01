#!/bin/bash
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
    for i in `seq 0 10 40`; do
        python3 fut_evaluation_script.py "$scenario" "$i" &
    done
done

scenarios=(
    "fut_NorESM_RCP45_Volc"
    "fut_NorESM_RCP45_VolcConst"
)

# Loop over scenarios and ranges in parallel
for scenario in "${scenarios[@]}"; do
    for i in `seq 0 10 50`; do
        python3 fut_evaluation_script.py "$scenario" "$i" &
    done
done


# Wait for all background processes to finish
wait

echo "All processes completed."
