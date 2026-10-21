#!/bin/bash
#SBATCH -J thorne-loo-fix
#SBATCH -t 6:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=24g
#SBATCH --array=0-1349%150
#SBATCH -o Results3_fix/logs/loo-%A_%a.out
# Submit from the repo root, after oscar_patch_newsamples.sh succeeded and the hashes matched:
#   sbatch fut_LOO/run_LOO_fix.sh
# One task per (selection mode, scenario, member): modes 7, 16, 32, 8 x 270 held-out members (tasks 0-1079),
# then mode 2 = individual reference methods EBMKF_ta4 / GWI_tot_CGWL via fut2_fix_worker.py (tasks 1080-1349).
# Writes Results3_fix/nmes{mode}run_exp{scenario}_r{member}.npz; the old Results3/ is untouched.
set -euo pipefail
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONBREAKPOINT=0
PY=/users/jnickla1/.conda/envs/cleanPy/bin/python
idx=${SLURM_ARRAY_TASK_ID:?Submit as an array job}
modes=(7 16 32 8 2)
mode=${modes[$((idx / 270))]}; r=$((idx % 270))
if (( r < 150 )); then scenario=$((r / 50)); member=$((r % 50))
else scenario=$((3 + (r - 150) / 60)); member=$(((r - 150) % 60)); fi
out=Results3_fix/nmes${mode}run_exp${scenario}_r${member}.npz
if [[ -f "$out" ]]; then echo "already done: $out"; exit 0; fi
export MPLCONFIGDIR="${TMPDIR:-/tmp}/thorne-mpl-${SLURM_JOB_ID}-${idx}"; mkdir -p "$MPLCONFIGDIR"
trap 'rm -rf -- "$MPLCONFIGDIR"' EXIT
echo "mode=$mode scenario=$scenario member=$member"
if (( mode == 2 )); then
    PYTHONPATH=. "$PY" -u fut_LOO/fut2_fix_worker.py "$scenario" "$member"
else
    PYTHONPATH=. "$PY" -u fut_LOO/LLOffcompute_fix_worker.py "$mode" "$scenario" "$member"
fi
test -f "$out" && echo "DONE $out"
