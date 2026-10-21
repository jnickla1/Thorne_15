# Rerunning the Table S4 leave-one-out on Oscar after the newsamples fix

Everything below uses these two variables. Set `REPO` to your Oscar checkout, which is the directory you ran the 8b61c6a1 ensemble from.

```bash
OSCAR=jnickla1@ssh.ccv.brown.edu
REPO=/path/to/Thorne_15_codefigurestats
```

## 1. Copy the code and small inputs to Oscar (run on the Mac)

`-R` keeps each relative path. `--backup` saves any file it replaces on Oscar as `*.pre_newsamples_fix`.

```bash
cd /Users/JohnMatthew/Downloads/Thorne_15_codefigurestats
rsync -avR --backup --suffix=.pre_newsamples_fix \
  calibrate_method_uncertainty.py \
  create_combination_method_LLOffcompute.py \
  create_comb_method_LLOffcompute_percase.py \
  calibrate_historical_combinations.py \
  custom_stack_utils.py \
  Methods/43_Forcing_Based/3_Human_Induced/eROF_anthro_method.py \
  Methods/43_Forcing_Based/3_Human_Induced/eROF_tot_method.py \
  averaged_runs126.csv averaged_runs245.csv averaged_runs370.csv averaged_runsVolc.csv averaged_runsVolcConst.csv \
  261018_HadCRUT5.1_var_unadjusted/current_methods_statistics_261018True.csv \
  thorne_2024v2025_investigation/patch_newsamples_scale.py \
  thorne_2024v2025_investigation/future_270_names.txt \
  thorne_2024v2025_investigation/oscar_loo_fix/ \
  "$OSCAR:$REPO/"
```

Do **not** copy `Results2/newsamples_scale_patch_manifest.json`. On Oscar it has to start empty; otherwise the patch would skip every file.

## 2. Patch the 270 future newsamples on Oscar, then compare hashes

On Oscar:

```bash
cd $REPO
mkdir -p Results3_fix/logs
sbatch thorne_2024v2025_investigation/oscar_loo_fix/oscar_patch_newsamples.sh
```

The job patches each file exactly once and checks that all 270 are done. It then writes SHA-256 hashes of every leave-one-out input to `Results3_fix/oscar_input_hashes.txt`: the 270 runs × 6 files, plus the selection CSVs and the worker. The log is `Results3_fix/logs/patch-<jobid>.out`. The log should end with `all 270 patched` and `wrote 1628 hashes`.

When the job has finished, run this on the Mac. The Mac's list is already written to `thorne_2024v2025_investigation/oscar_loo_fix/mac_input_hashes.txt`.

```bash
cd /Users/JohnMatthew/Downloads/Thorne_15_codefigurestats
rsync -av "$OSCAR:$REPO/Results3_fix/oscar_input_hashes.txt" thorne_2024v2025_investigation/oscar_loo_fix/
diff thorne_2024v2025_investigation/oscar_loo_fix/mac_input_hashes.txt thorne_2024v2025_investigation/oscar_loo_fix/oscar_input_hashes.txt && echo "IDENTICAL INPUTS"
```

Only go on to step 3 if it prints `IDENTICAL INPUTS`. If any lines differ, send me the diff before submitting anything.

## 3. Submit the array (on Oscar)

```bash
cd $REPO
sbatch thorne_2024v2025_investigation/oscar_loo_fix/run_LOO_fix.sh
```

The array is 1,080 tasks, covering modes 7, 16, 32 and 8 × 270 held-out members, with up to 150 running at once and 24 GB each. Each task writes `Results3_fix/nmes{mode}run_exp{s}_r{m}.npz` and skips any case that already has an output, so resubmitting after failures is safe. The old `Results3/` is left untouched.

Mode 2, the individual-method comparison, isn't rerun. It doesn't use newsamples, so the earlier results still hold.

## 4. Bring the results back and build Table S4 (on the Mac)

```bash
cd /Users/JohnMatthew/Downloads/Thorne_15_codefigurestats
rsync -av --include='*.npz' --exclude='*' "$OSCAR:$REPO/Results3_fix/" Results3_fix_oscar/
/Users/JohnMatthew/opt/anaconda3/envs/cleanPy/bin/python thorne_2024v2025_investigation/oscar_loo_fix/check_and_aggregate.py Results3_fix_oscar thorne_2024v2025_investigation/oscar_loo_fix/table_S4_out
```

The script does three things:

- **Coverage:** checks that 1080/1080 files are present and all values are finite.
- **Cross-check:** compares the 321 cases the stopped local run had already finished (7-method ×270, 18-method ×50, 8-method ×1). It prints the maximum absolute differences, which should be zero or at floating-point level.
- **Aggregation:** writes `table-S4.md`, the numeric and formatted CSVs, and the `final_combined_*` CSVs.

## What was verified on the Mac

- `LLOffcompute_fix_worker.py` is `create_combination_method_LLOffcompute.py` with three changes: it reads the 261018 statistics file, writes to `Results3_fix/`, and saves `n_valid_kl`.
- Run standalone on the Mac, it reproduced two local pipeline cases (7-method exp3 r17, 18-method exp0 r5) with zero difference in every field.
