#!/bin/bash
#SBATCH -J thorne-fut-b1900
#SBATCH -t 12:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=12g
#SBATCH --array=0-419%150
#SBATCH -o Results2/rerun_b1900/logs/fut-%A_%a.out
# Future-ensemble rerun with the 1850-1900 (inclusive) preindustrial datum (2026-10-21).
# Same protocol as the 18 Oct 8b61c6a1 run: 270 calibrated members (MPI SSP126/245/370 x 50, NorESM Volc/VolcConst x 60)
# plus 150 evaluation-only MPI VolcConst members. Completion records go to Results2/rerun_b1900/ (old records untouched).
# Submit from the repo root after: mkdir -p Results2/rerun_b1900/logs, and moving the future Bayes-CP / GAM caches aside.
set -euo pipefail
module load matlab/R2025a-bcih
export MATLAB_BIN=matlab
module load r/4.6.0-xnxx
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONBREAKPOINT=0
unset THORNE_REFRESH_METHODS
# GAM_AR0/AR1 refit in R (their future caches were moved aside): fail fast if a required R package is missing.
Rscript --vanilla -e 'library(mgcv); library(this.path)' > /dev/null
PY=/users/jnickla1/.conda/envs/cleanPy/bin/python
# Fail before any writes unless the 1850-1900 future datum and the regenerated ensemble-mean targets are present.
grep -q "this_hsim_yr\[0:51\]" evaluate_future_scenarios.py
grep -q "this_hsim_yr\[0:51\]" calibrate_method_uncertainty.py
grep -q "this_hsim_yr\[0:51\]" fut_evaluation_gen_ensemble.py
grep -q "^    reference_stop = 51$" calibrate_method_uncertainty.py
grep -q "reference_stop = 51 # 1850–1900 inclusive for every experiment" Methods/42_Temp_Alone/5_Kalman/Kalman_RW_method.py
for e in fut_ESM1-2-LR_SSP126_constVolc fut_ESM1-2-LR_SSP245_constVolc fut_ESM1-2-LR_SSP370_constVolc fut_NorESM_RCP45_Volc fut_NorESM_RCP45_VolcConst; do
    test -s "Results/ensemble_mean_${e}.csv"
done
idx=${SLURM_ARRAY_TASK_ID:?Submit as an array job}
if (( idx < 150 )); then
    scenarios=(fut_ESM1-2-LR_SSP126_constVolc fut_ESM1-2-LR_SSP245_constVolc fut_ESM1-2-LR_SSP370_constVolc)
    scenario=${scenarios[$((idx / 50))]}; member=$((idx % 50)); calibration=1
elif (( idx < 270 )); then
    scenarios=(fut_NorESM_RCP45_Volc fut_NorESM_RCP45_VolcConst)
    scenario=${scenarios[$(((idx-150) / 60))]}; member=$(((idx-150) % 60)); calibration=1
elif (( idx < 420 )); then
    scenarios=(fut_ESM1-2-LR_SSP126_VolcConst fut_ESM1-2-LR_SSP245_VolcConst fut_ESM1-2-LR_SSP370_VolcConst)
    scenario=${scenarios[$(((idx-270) / 50))]}; member=$(((idx-270) % 50)); calibration=0
else
    echo "Unexpected array index: $idx" >&2; exit 2
fi
export MPLCONFIGDIR="${TMPDIR:-/tmp}/thorne-mpl-${SLURM_JOB_ID}-${idx}"
mkdir -p "$MPLCONFIGDIR" Results Results2/rerun_b1900
trap 'rm -rf -- "$MPLCONFIGDIR"' EXIT
printf 'Source commit: '; git rev-parse HEAD
printf 'Scenario: %s; member: %s; recalibrate: %s\n' "$scenario" "$member" "$calibration"
if (( calibration )); then
    "$PY" -u calibrate_method_uncertainty.py "$scenario" "$member"
fi
# Force the full evaluation set in memory; do not modify the historical source.
"$PY" -u - "$scenario" "$member" <<'RUN_COMPLETE_EVALUATION'
import runpy,sys
import evaluate_historical_observations as historical
historical.running_subset=('Methods/42_Temp_Alone','Methods/43_Forcing_Based','Methods/44_EarthModel_CGWL')
sys.argv=['evaluate_future_scenarios.py',sys.argv[1],sys.argv[2]]
runpy.run_path('evaluate_future_scenarios.py',run_name='__main__')
RUN_COMPLETE_EVALUATION
"$PY" - "$scenario" "$member" "$calibration" <<'VERIFY_CASE'
from pathlib import Path
import sys,json,os,time,numpy as np,pandas as pd,pickle
sc,member,cal=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]);prefix=sc+str(member)
if cal:
    table=pd.read_csv(f'Results2/{prefix}_names_var_scale.csv')
    assert table.method_name.is_unique
    for key in ['standard','standard_se','central_est','hermguass_lls','newsamples']:
        a=np.load(f'Results2/{prefix}_{key}.npy',mmap_mode='r');assert not np.isinf(a).any(),key
    central=np.load(f'Results2/{prefix}_central_est.npy',mmap_mode='r')
    assert len(table)==len(central)
    for method in ['EBMKF_ta','EBMKF_ta2','EBMKF_ta4']:
        assert np.isfinite(central[table.method_name.tolist().index(method)]).sum()>100
with open(f'Results/results{prefix}.pickle','rb') as f:result=pickle.load(f)
for method in ['EBMKF_ta','EBMKF_ta2','EBMKF_ta4','removeMEI_volc_cons','removeMEI_volc_refit']:
    assert method in result,method
    assert np.isfinite(result[method]['LT_trend'][0]).sum()>100,method
for kind in ['current','all']:
    table=pd.read_csv(f'Results/{kind}_fut_statistics_{prefix}.csv')
    assert {'EBMKF_ta','EBMKF_ta2','EBMKF_ta4'}<=set(table.method_name)
Path(f'Results2/rerun_b1900/{prefix}.json').write_text(json.dumps({'scenario':sc,'member':member,'calibration':bool(cal),'evaluation':True,'baseline':'1850-1900','finished_at':time.time(),'job_id':os.environ['SLURM_JOB_ID']},indent=2))
print('VERIFIED COMPLETE:',prefix)
VERIFY_CASE
