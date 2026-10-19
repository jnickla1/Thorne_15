#!/bin/bash
# Usage (from the repo root): bash thorne_2024v2025_investigation/oscar_loo_fix/hash_inputs.sh OUTFILE
# SHA-256 of every input the leave-one-out reads, in a fixed order, so the Mac and Oscar lists can be diffed.
set -euo pipefail
out=${1:?output file}
if command -v sha256sum >/dev/null; then H="sha256sum"; else H="shasum -a 256"; fi
names=thorne_2024v2025_investigation/future_270_names.txt
{
  for f in averaged_runs126.csv averaged_runs245.csv averaged_runs370.csv averaged_runsVolc.csv averaged_runsVolcConst.csv \
           261018_HadCRUT5.1_var_unadjusted/current_methods_statistics_261018True.csv custom_stack_utils.py \
           thorne_2024v2025_investigation/oscar_loo_fix/LLOffcompute_fix_worker.py; do $H "$f"; done
  while read -r n; do
    for k in names_var_scale.csv standard.npy standard_se.npy central_est.npy hermguass_lls.npy newsamples.npy; do
      $H "Results2/${n}_${k}"
    done
  done < "$names"
} > "$out"
echo "wrote $(wc -l < "$out") hashes to $out"
