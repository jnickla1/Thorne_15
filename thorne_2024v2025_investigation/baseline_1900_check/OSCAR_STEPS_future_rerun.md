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

## 1. On Oscar: pull commit 216d3fde ("Consistent 1850-1900 baseline") from GitHub

Earlier rsyncs left untracked copies of the same files on Oscar (fut_LOO/, baseline_1900_check/ and so on). Move them aside, set any local edits to tracked files aside, then fast-forward.

```bash
cd ~/data/jnickla1/Thorne_15_codefigurestats && git fetch origin && git status --short | grep -v '^??' | head -20
```

```bash
mkdir -p pre_pull_backup && for f in $(git diff --name-only HEAD origin/main); do if [ -e "$f" ] && ! git ls-files --error-unmatch "$f" >/dev/null 2>&1; then mkdir -p "pre_pull_backup/$(dirname "$f")"; mv "$f" "pre_pull_backup/$f"; fi; done
```

```bash
git stash push -m "Oscar local edits before 216d3fde" ; git pull --ff-only origin main && git log --oneline -1
```

The last command should print `216d3fde Consistent 1850-1900 baseline`. `git stash list` keeps the set-aside local edits.

## 2. On Oscar: move the stale future caches aside and create the log folder

The Bayes change-point and GAM caches store fitted temperature levels on the old datum. The quartic, cubic-spline, Butterworth and LFCA caches don't depend on the datum, so they stay.

```bash
cd ~/data/jnickla1/Thorne_15_codefigurestats && BK=cache_backup_baseline1850-1899 && mkdir -p $BK/BSC $BK/GAM Results2/rerun_b1900/logs
```

```bash
mv Methods/42_Temp_Alone/2_LT_Fits/Bayes_Sequential_Change_Point/output/fut_* $BK/BSC/ ; mv Methods/42_Temp_Alone/4_GAM_AR1/GAM_AR_Stephenson/output/gamAR*_fut_* $BK/GAM/ ; ls $BK/BSC | wc -l ; ls $BK/GAM | wc -l
```

## 3. On Oscar: submit

```bash
sbatch run_future_ensemble_1850-1900.sh
```

The script stops before writing anything unless the 1850–1900 datum and the five new ensemble-mean target files are present.

## 4. When it finishes, on Oscar

The first command should print 420.

```bash
ls Results2/rerun_b1900/*.json | wc -l
```

```bash
sacct -j <jobid> -X --format=JobID%18,State,Elapsed | grep -v COMPLETED
```

Resubmit any failed tasks with `sbatch --array=<ids> run_future_ensemble_1850-1900.sh`.

## 5. Next, decided after step 4

1. Regenerate `averaged_runs*.csv` and the figures.
2. Choose the historical statistics table for method selection. The new `current_methods_statistics_261020True.csv` moves GWI_tot_AR6 just across the RMS 0.06 cut (0.05998 → 0.06004), which drops it from the 18- and 33-method sets.
3. Move `Results3_fix/*.npz` aside.
4. Run `fut_LOO/run_LOO_fix.sh`: 1,350 tasks, including the individual-method rows.
5. Compare Figure 24, Table 3 and Table S4 with the current numbers.
