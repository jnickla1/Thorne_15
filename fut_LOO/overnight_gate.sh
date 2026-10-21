#!/bin/bash
#SBATCH -J thorne-gate-b1900
#SBATCH -t 2:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=8g
#SBATCH -o Results2/rerun_b1900/logs/gate-%j.out
# Automatic chaining for the overnight 1850-1900 rerun (queued by fut_LOO/submit_overnight_b1900.sh).
#   overnight_gate.sh ensemble RETRY  -- runs after the 420-task future ensemble finishes (success or not):
#       missing members -> resubmit them once with a 16 h limit and re-queue this gate;
#       all 420 complete -> regenerate averaged_runs*.csv (LOO method selection), move old Results3_fix outputs aside,
#                           submit fut_LOO/run_LOO_fix.sh and queue the LOO gate.
#   overnight_gate.sh loo RETRY       -- runs after the 1,350-task LOO array: missing cases -> resubmit once with 24 h;
#                                        all present -> write Results2/rerun_b1900/OVERNIGHT_COMPLETE.txt.
# Status lines are appended to Results2/rerun_b1900/overnight_status.txt.
set -euo pipefail
STAGE=${1:?ensemble or loo}; RETRY=${2:-0}
PY=/users/jnickla1/.conda/envs/cleanPy/bin/python
STATUS=Results2/rerun_b1900/overnight_status.txt
log() { echo "$(date '+%F %T') [$STAGE retry=$RETRY] $*" | tee -a "$STATUS"; }

if [[ "$STAGE" == ensemble ]]; then
    missing=$("$PY" - <<'EOF'
from pathlib import Path
miss = []
for idx in range(420):
    if idx < 150:
        sc = ['fut_ESM1-2-LR_SSP126_constVolc', 'fut_ESM1-2-LR_SSP245_constVolc', 'fut_ESM1-2-LR_SSP370_constVolc'][idx // 50]; m = idx % 50
    elif idx < 270:
        sc = ['fut_NorESM_RCP45_Volc', 'fut_NorESM_RCP45_VolcConst'][(idx - 150) // 60]; m = (idx - 150) % 60
    else:
        sc = ['fut_ESM1-2-LR_SSP126_VolcConst', 'fut_ESM1-2-LR_SSP245_VolcConst', 'fut_ESM1-2-LR_SSP370_VolcConst'][(idx - 270) // 50]; m = (idx - 270) % 50
    if not Path(f'Results2/rerun_b1900/{sc}{m}.json').exists():
        miss.append(str(idx))
print(','.join(miss))
EOF
)
    if [[ -n "$missing" ]]; then
        if (( RETRY < 1 )); then
            jid=$(sbatch --parsable --array="$missing" -t 16:00:00 run_future_ensemble_1850-1900.sh)
            gid=$(sbatch --parsable --dependency=afterany:"$jid" fut_LOO/overnight_gate.sh ensemble 1)
            log "missing ensemble tasks [$missing]: resubmitted as $jid (16 h); gate re-queued as $gid"
            exit 0
        fi
        log "FAILED: ensemble tasks still missing after one retry: [$missing]. LOO not submitted."
        exit 1
    fi
    log "all 420 ensemble members complete"
    MPLBACKEND=Agg "$PY" -u plot_fut_results.py > Results2/rerun_b1900/logs/plot_fut_results.log 2>&1
    log "regenerated averaged_runs*.csv and combined7histogramTrue.png"
    old=Results3_fix/old_before_b1900_$(date +%Y%m%d_%H%M%S)
    mkdir -p Results3_fix/logs "$old"
    shopt -s nullglob; prev=(Results3_fix/nmes*.npz); shopt -u nullglob
    if (( ${#prev[@]} )); then mv "${prev[@]}" "$old/"; fi
    log "moved ${#prev[@]} earlier LOO outputs to $old"
    jid=$(sbatch --parsable fut_LOO/run_LOO_fix.sh)
    gid=$(sbatch --parsable --dependency=afterany:"$jid" fut_LOO/overnight_gate.sh loo 0)
    log "submitted LOO array $jid (1,350 tasks); LOO gate queued as $gid"
    exit 0
fi

if [[ "$STAGE" == loo ]]; then
    missing=$("$PY" - <<'EOF'
from pathlib import Path
miss = []
for idx in range(1350):
    mode = [7, 16, 32, 8, 2][idx // 270]; r = idx % 270
    s, m = (r // 50, r % 50) if r < 150 else (3 + (r - 150) // 60, (r - 150) % 60)
    if not Path(f'Results3_fix/nmes{mode}run_exp{s}_r{m}.npz').exists():
        miss.append(str(idx))
print(','.join(miss))
EOF
)
    if [[ -n "$missing" ]]; then
        if (( RETRY < 1 )); then
            jid=$(sbatch --parsable --array="$missing" -t 24:00:00 fut_LOO/run_LOO_fix.sh)
            gid=$(sbatch --parsable --dependency=afterany:"$jid" fut_LOO/overnight_gate.sh loo 1)
            log "missing LOO tasks [$missing]: resubmitted as $jid (24 h); gate re-queued as $gid"
            exit 0
        fi
        log "FAILED: LOO tasks still missing after one retry: [$missing]"
        exit 1
    fi
    n=$(ls Results3_fix/nmes*.npz | wc -l)
    log "all 1,350 LOO cases present ($n files). Copy Results3_fix/*.npz back to the Mac for aggregation."
    date > Results2/rerun_b1900/OVERNIGHT_COMPLETE.txt
    exit 0
fi
echo "unknown stage $STAGE" >&2; exit 2
