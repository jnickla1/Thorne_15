#!/bin/bash
source activate cleanpy
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

# Wait for all background processes to finish
wait

echo "All processes completed."
