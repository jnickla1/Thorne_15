#!/bin/bash
#SBATCH -J thorne-patch-newsamples
#SBATCH -t 2:00:00
#SBATCH --account=epscor-condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=8g
#SBATCH -o Results3_fix/logs/patch-%j.out
# Submit from the repo root:  sbatch fut_LOO/newsamples_patch/oscar_patch_newsamples.sh
# Applies the exact newsamples scale fix (c + (old-c)*s^2 on sampled-method rows) to the 270 future sets,
# then hashes all leave-one-out inputs for comparison with the Mac. Safe to resubmit: patched files are skipped
# via Results2/newsamples_scale_patch_manifest.json (never copy the Mac's manifest to Oscar).
set -euo pipefail
PY=/users/jnickla1/.conda/envs/cleanPy/bin/python
xargs "$PY" fut_LOO/newsamples_patch/patch_newsamples_scale.py < fut_LOO/newsamples_patch/future_270_names.txt
"$PY" -c "import json;d=json.load(open('Results2/newsamples_scale_patch_manifest.json'));n=[l.strip() for l in open('fut_LOO/newsamples_patch/future_270_names.txt')];assert all(x in d for x in n),'not all patched';print('all 270 patched')"
bash fut_LOO/newsamples_patch/hash_inputs.sh Results3_fix/oscar_input_hashes.txt
