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

## 1. On the Mac: copy the code, the new targets and the array script to Oscar

Files that Oscar replaces are kept as `*.pre_b1900`.

```bash
cd /Users/JohnMatthew/Downloads/Thorne_15_codefigurestats && rsync -avR --backup --suffix=.pre_b1900 calibrate_method_uncertainty.py evaluate_future_scenarios.py evaluate_future_scenario_subsets.py fut_evaluation_gen_ensemble.py evaluate_historical_observations.py run_future_ensemble_1850-1900.sh Methods/42_Temp_Alone/3_ST_Fits/min_month_proj_method.py Methods/42_Temp_Alone/5_Kalman/Kal_flexLin_ocn_method.py Methods/42_Temp_Alone/5_Kalman/Kalman_RW_method.py Methods/42_Temp_Alone/5_Kalman/Kalman_RW_ocn_method.py Methods/42_Temp_Alone/6_Remove_IV/lfca_SST_method.py Methods/42_Temp_Alone/6_Remove_IV/lfca_TAS_method.py Methods/43_Forcing_Based/1_ERF_FaIR/FaIR_all_unB_method.py Methods/43_Forcing_Based/1_ERF_FaIR/FaIR_anthro_unB_method.py Methods/43_Forcing_Based/1_ERF_FaIR/FaIR_comb_unB_method.py Methods/43_Forcing_Based/1_ERF_FaIR/FaIR_nonat_unB_method.py Methods/43_Forcing_Based/2_Kalman/EBMKF_ta4_attibutable.py Methods/43_Forcing_Based/3_Human_Induced/GWI_anthro_AR6_method.py Methods/43_Forcing_Based/3_Human_Induced/GWI_anthro_CGWL_method.py Methods/43_Forcing_Based/3_Human_Induced/GWI_anthro_SR15_method.py Methods/43_Forcing_Based/3_Human_Induced/GWI_anthro_method.py Methods/43_Forcing_Based/3_Human_Induced/GWI_tot_AR6_method.py Methods/43_Forcing_Based/3_Human_Induced/GWI_tot_CGWL_method.py Methods/43_Forcing_Based/3_Human_Induced/GWI_tot_SR15_method.py Methods/43_Forcing_Based/3_Human_Induced/GWI_tot_method.py Methods/43_Forcing_Based/3_Human_Induced/KCC_all_method.py fut_LOO/fut_2_script.py fut_LOO/fut2_fix_worker.py fut_LOO/LLOffcompute_fix_worker.py fut_LOO/create_combination_method_LLOffcompute.py fut_LOO/make_oscar_worker.py fut_LOO/run_LOO_fix.sh fut_LOO/check_and_aggregate.py fut_LOO/selection.json Results/ensemble_mean_fut_ESM1-2-LR_SSP126_constVolc.csv Results/ensemble_mean_fut_ESM1-2-LR_SSP245_constVolc.csv Results/ensemble_mean_fut_ESM1-2-LR_SSP370_constVolc.csv Results/ensemble_mean_fut_NorESM_RCP45_Volc.csv Results/ensemble_mean_fut_NorESM_RCP45_VolcConst.csv "jnickla1@ssh.ccv.brown.edu:/users/jnickla1/data/jnickla1/Thorne_15_codefigurestats/"
```

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
