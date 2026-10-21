#!/bin/bash
#SBATCH -J thorne-loo-b1900
#SBATCH -t 12:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=24g
#SBATCH --array=0-809%150
#SBATCH -o thorne_2024v2025_investigation/baseline_1900_check/loo_logs/b1900-%A_%a.out
# Read-only baseline check (1850-1900 vs 1850-1899): reruns the Table S4 leave-one-out fits for modes 7, 16, 32
# and evaluates the first 1.5 °C crossing at 1.5 and at 1.5 + c. Writes only to baseline_1900_check/loo_shift/.
# Submit from the repo root:  sbatch thorne_2024v2025_investigation/baseline_1900_check/run_loo_shift_oscar.sh
set -euo pipefail
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONBREAKPOINT=0
PY=/users/jnickla1/.conda/envs/cleanPy/bin/python
B=thorne_2024v2025_investigation/baseline_1900_check
idx=${SLURM_ARRAY_TASK_ID:?Submit as an array job}
modes=(7 16 32)
mode=${modes[$((idx / 270))]}; r=$((idx % 270))
if (( r < 150 )); then scenario=$((r / 50)); member=$((r % 50))
else scenario=$((3 + (r - 150) / 60)); member=$(((r - 150) % 60)); fi
out=$B/loo_shift/nmes${mode}run_exp${scenario}_r${member}.npz
if [[ -f "$out" ]]; then echo "already done: $out"; exit 0; fi
export MPLCONFIGDIR="${TMPDIR:-/tmp}/thorne-mpl-${SLURM_JOB_ID}-${idx}"; mkdir -p "$MPLCONFIGDIR"
trap 'rm -rf -- "$MPLCONFIGDIR"' EXIT
echo "mode=$mode scenario=$scenario member=$member"
PYTHONPATH=. "$PY" -u $B/LLOffcompute_shift_worker.py "$mode" "$scenario" "$member"
test -f "$out" && echo "DONE $out"
