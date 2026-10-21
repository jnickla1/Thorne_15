#!/bin/zsh
# Local historical rerun after the datum-correction commit (acf2df22). Run from the repo root.
cd /Users/JohnMatthew/Downloads/Thorne_15_codefigurestats
L=thorne_2024v2025_investigation/baseline_1900_check/hist_rerun_logs
export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
PY=/Users/JohnMatthew/opt/anaconda3/envs/cleanPy/bin/python
step() { echo "START $1 $(date)" >> $L/status.txt; shift; "$@" >> $L/$(date +%s).log 2>&1 && echo "OK $(date)" >> $L/status.txt || { echo "FAILED $(date)" >> $L/status.txt; exit 1; } }
step calibrate $PY -u calibrate_method_uncertainty.py historical
step combinations $PY -u calibrate_historical_combinations.py
step evaluate $PY -u evaluate_historical_observations.py
echo "START histens $(date)" >> $L/status.txt
seq 0 5 100 | xargs -P 16 -n 1 ./run_histens_local.sh
echo "DONE histens $(date)" >> $L/status.txt
