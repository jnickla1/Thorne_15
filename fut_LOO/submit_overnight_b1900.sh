#!/bin/bash
# One command for the overnight 1850-1900 rerun on Oscar. Run on the login node from the repo root, after `git pull`:
#   bash fut_LOO/submit_overnight_b1900.sh
# It moves the stale future Bayes-CP / GAM caches aside (fitted levels on the old datum), submits the 420-task
# future ensemble, and queues fut_LOO/overnight_gate.sh, which verifies the ensemble (retrying failures once),
# regenerates averaged_runs*.csv, submits the 1,350-task LOO array and verifies that too (retrying once).
# Progress: Results2/rerun_b1900/overnight_status.txt; success marker: Results2/rerun_b1900/OVERNIGHT_COMPLETE.txt.
set -euo pipefail
test -f run_future_ensemble_1850-1900.sh && test -f fut_LOO/run_LOO_fix.sh && test -f datum_shift.py
grep -q "this_hsim_yr\[0:51\]" evaluate_future_scenarios.py
mkdir -p Results2/rerun_b1900/logs Results3_fix/logs
BK=cache_backup_baseline1850-1899; mkdir -p "$BK/BSC" "$BK/GAM"
shopt -s nullglob
bsc=(Methods/42_Temp_Alone/2_LT_Fits/Bayes_Sequential_Change_Point/output/fut_*)
gam=(Methods/42_Temp_Alone/4_GAM_AR1/GAM_AR_Stephenson/output/gamAR*_fut_*)
shopt -u nullglob
if (( ${#bsc[@]} )); then mv "${bsc[@]}" "$BK/BSC/"; fi
if (( ${#gam[@]} )); then mv "${gam[@]}" "$BK/GAM/"; fi
echo "moved ${#bsc[@]} Bayes-CP and ${#gam[@]} GAM future cache files to $BK/"
jid=$(sbatch --parsable run_future_ensemble_1850-1900.sh)
gid=$(sbatch --parsable --dependency=afterany:"$jid" fut_LOO/overnight_gate.sh ensemble 0)
echo "$(date '+%F %T') submitted future ensemble $jid (420 tasks); gate $gid queued (afterany)" | tee -a Results2/rerun_b1900/overnight_status.txt
