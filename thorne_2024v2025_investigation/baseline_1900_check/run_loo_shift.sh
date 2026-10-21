#!/bin/zsh
# Runs the read-only shift-check worker for every case not yet done, 8 at a time. From the repo root.
cd /Users/JohnMatthew/Downloads/Thorne_15_codefigurestats
B=thorne_2024v2025_investigation/baseline_1900_check
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=. PYTHONDONTWRITEBYTECODE=1
cat $B/loo_cases.txt | xargs -P 8 -L 1 sh -c '[ -f '$B'/loo_shift/nmes$0run_exp$1_r$2.npz ] || /Users/JohnMatthew/opt/anaconda3/envs/cleanPy/bin/python -u '$B'/LLOffcompute_shift_worker.py $0 $1 $2 > '$B'/loo_logs/$0_$1_$2.log 2>&1'
echo "ALL DONE $(date)" >> $B/loo_logs/_status.txt
