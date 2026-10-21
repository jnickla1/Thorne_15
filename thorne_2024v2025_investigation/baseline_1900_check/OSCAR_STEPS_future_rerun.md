# Future-ensemble rerun on the 1850–1900 baseline (option A, 2026-10-21)

## Code changes (Mac, uncommitted)

**Future datum, `[0:50]` → `[0:51]`:**
- `evaluate_future_scenarios.py`: observation baseline (NorESM alignment) and MPI `offset_sim`.
- `fut_evaluation_gen_ensemble.py`: same two places.
- `evaluate_future_scenario_subsets.py`: same two places (Figures 22/23).
- `calibrate_method_uncertainty.py`: MPI `offset_sim`, plus the `reference_stop = 51 if historical else 50` gate, now always 51.
- `fut_LOO/fut_2_script.py`: observation baseline and `offset_sim`. `fut2_fix_worker.py` was regenerated from it.

**Method gates that kept the future runs on 1850–1899, now always 51:**
- `min_month_proj`, `Kal_flexLin_ocn`, `Kalman_RW` and `Kalman_RW_ocn` (`reference_stop`).
- `lfca_TAS`: the ESM `preind_end_idx`.

**Deliberately left at `[0:50]`:**
- EBMKF internals, as John asked.
- Kal_flexLin's internal centering, which is added back and is also unchanged for historical.
- removeMEI's `given_preind_base`, which is unused.
- The ocean-record alignment in Kal_flexLin_ocn and _RW_ocn.
- CGWL's `histens`-only branch.

**Ensemble-mean targets:** `Results/ensemble_mean_fut_*.csv` were regenerated on the Mac with the new datum, shifting by −0.00113 °C (MPI) and −0.00244 °C (NorESM). The old files are in `Results/ensemble_mean_backup_baseline1850-1899/`.

**Not synced:** the uncommitted submodule edits (lfca, KCC, AWA, forced-patterns, elastic-net) predate this work and are left as they are on Oscar.

## 0. On Oscar: don't run the shift-check array

The full rerun supersedes it. If `run_loo_shift_oscar.sh` was submitted, `scancel` it.

## 1. On Oscar: pull from GitHub (needs commit acf2df22, "Datum correction for anchored methods", or later)

Earlier rsyncs left untracked copies of the same files on Oscar (fut_LOO/, baseline_1900_check/ and so on). Move them aside, set any local edits to tracked files aside, then fast-forward.

```bash
cd ~/data/jnickla1/Thorne_15_codefigurestats && git fetch origin && git status --short | grep -v '^??' | head -20
```

```bash
mkdir -p pre_pull_backup && for f in $(git diff --name-only HEAD origin/main); do if [ -e "$f" ] && ! git ls-files --error-unmatch "$f" >/dev/null 2>&1; then mkdir -p "pre_pull_backup/$(dirname "$f")"; mv "$f" "pre_pull_backup/$f"; fi; done
```

```bash
git stash push -m "Oscar local edits before the 1850-1900 rerun" ; git pull --ff-only origin main && git log --oneline -1
```

The last command should print `acf2df22 Datum correction for anchored methods under the 1850-1900 baseline`, or a later commit. `git stash list` keeps the set-aside local edits.

## 2. On Oscar: start the whole overnight chain with one command

```bash
bash fut_LOO/submit_overnight_b1900.sh
```

The command does the following:
1. Moves the stale future Bayes change-point and GAM caches into `cache_backup_baseline1850-1899/`. Those caches hold fitted levels on the old datum.
2. Submits the 420-task future ensemble (`run_future_ensemble_1850-1900.sh`).
3. Queues `fut_LOO/overnight_gate.sh` to run once the ensemble finishes, whether or not every task succeeded. The gate then:
   - checks all 420 completion records and resubmits any missing members once, with a 16 h limit;
   - once all 420 are present, regenerates `averaged_runs*.csv`, which the leave-one-out uses for method selection;
   - moves the earlier `Results3_fix/*.npz` aside;
   - submits the 1,350-task leave-one-out array, covering modes 7, 16, 32, 8 and 2 (the individual-method rows);
   - checks those results the same way, resubmitting any missing cases once with a 24 h limit.

Progress goes to `Results2/rerun_b1900/overnight_status.txt`. When everything has succeeded, `Results2/rerun_b1900/OVERNIGHT_COMPLETE.txt` appears.

## 3. In the morning

```bash
cat Results2/rerun_b1900/overnight_status.txt
```

Then copy the results back to the Mac for aggregation and comparison:

```bash
rsync -av --include='*.npz' --exclude='*' "jnickla1@ssh.ccv.brown.edu:/users/jnickla1/data/jnickla1/Thorne_15_codefigurestats/Results3_fix/" fut_LOO/results/Results3_fix_b1900/
```

Also copy the small future outputs (`Results/*fut*` CSVs, the `averaged_runs*.csv` files and `combined7histogramTrue.png`). These rebuild Figure 24, Table 3 and Table S4.
