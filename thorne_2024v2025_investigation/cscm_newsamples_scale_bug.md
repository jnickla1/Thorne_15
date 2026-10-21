# CSCM: `newsamples` for sampled methods use the wrong scale (1/s instead of s)

Date: 2026-10-19. This investigation was read-only: no main code was changed and nothing was rerun.

## Verdict: this is a bug

In `calibrate_method_uncertainty.py`, the sampled (dict-returning, non-Gaussian) branch of the newsamples code does the following:

```python
# L377-379: likelihood (used for best_alter_scale and the hermgauss lls)
resc = c + (y - c)/s ;  ll = log p(resc) - log(s)        # density of  c + s*(Y - c)
# L380: samples
newsamples = (resample - c)/s + c                         # draws of  c + (Y - c)/s   <-- inverse
```

Change of variables gives p_s(y) = p(c + (y−c)/s)/s, which is the density of **c + s·(Y−c)**. The likelihood, the `best_alter_scale` fit, the stacking lls (`*_hermguass_lls.npy`) and `evaluate_heads_tails_ensemble.py:386-389` all use this convention. The Gaussian branch uses it too, with `se*s`. The only line that inverts it is the sampled-method newsamples line. As a result, the sampled-method spreads are off by a factor of **s²**: methods with s<1 get wider when they should get narrower, and methods with s>1 get narrower.

### Empirical confirmation from the stored `Results2/historical_newsamples.npy`

I used `cscm_newsamples_diag.py` (output in `cscm_newsamples_diag.log`) over the fit window, the last 100 years excluding the final 10. It computes the RMS of (standard − central)/SD(newsamples):

- **All Gaussian methods:** rms_z = 1.00 ± 0.01, which is calibrated as intended.
- **Sampled methods:** rms_z ≈ s², for example:

| method | s | s² | rms_z | 90 % coverage |
|---|---|---|---|---|
| GWI_tot_CGWL | 0.32 | 0.10 | 0.10 | 1.00 |
| CGWL10y_pUKCP | 0.27 | 0.07 | 0.07 | 1.00 |
| GWI_anthro_CGWL | 0.46 | 0.21 | 0.21 | 1.00 |
| FaIR_comb_unB | 0.77 | 0.59 | 0.59 | 1.00 |
| CGWL10y_sfUKCP | 0.80 | 0.64 | 0.62 | 0.99 |
| FaIR_anthro | 1.22 | 1.49 | 1.51 | 0.65 |
| KCC_all | 1.37 | 1.88 | 1.95 | 0.69 |
| CGWL10y_for_halfU, GWI_tot_AR6 | ≈1.0 | ≈1.0 | 1.0 | ok, essentially unaffected |

The affected families are FaIR_*, GWI_*, KCC_*, eROF_* and CGWL10y_*.

This matches the dashboard's read-only refit (`climate-dashboard/data/supplementary/dashboard_combination_calibration.json`) for the 7-method CSCM:

| samples | nsharp | KL | mean SD |
|---|---|---|---|
| stored | 2.00 | 119.4 | 0.092 |
| rescaled by s² | 1.03 | 78.4 | 0.032 |

The 0.032 mean SD is about the same as the PIVW calibrated scale (0.0323). The stored samples are over-dispersed because the top-tier sampled methods (FaIR_comb_unB, GWI_tot_CGWL, CGWL10y_sfUKCP) all have s<1. `nsharp` (which is bounded to [1, 10]) partly made up for that by averaging about 2 draws, but it could only do this uniformly across methods. That leaves the CSCM mixture too wide and gives a worse KL.

**Side observation, not part of this bug:** eROF_anthro and eROF_tot have 90 % coverage of 0.00 and 0.02 despite rms_z < 1. That suggests their resample draws are offset from their `mean`, or are degenerate. It's worth a separate look.

## Outputs this touches

- **`Results2/historical_newsamples.npy`** and all **`Results2/fut_*_newsamples.npy`** (272 files). The rows for sampled methods are wrong. The Gaussian rows, `*_central_est.npy`, `*_hermguass_lls.npy` and `*_names_var_scale.csv` are all **correct**.
- **Historical CSCM** (`calibrate_historical_combinations.py`): `best_nsharp`, `sharpened_blended`, `best_stack_KL` and the CSCM crossing-year and abs-error numbers. The stacking weights are unaffected because they come from the lls. PIVW is unaffected (`comb_*_scales.npy`, `best_scale_iv`).
- **Future CSCM, i.e. Table S4 leave-one-out** (`create_combination_method_LLOffcompute.py` and `create_comb_method_LLOffcompute_percase.py` → `Results3/nmes*run_exp*_r*.npz` → `create_combination_method_LLOffcombine.py`): this covers the `sharpened_blended` column only. The `inverse_variance` column is unaffected.
- **Any CSCM curves or histograms that are built from those**, e.g. the stacked variant of `plot_fut_results.py` → `combined7histogram*.png`, if it reads the sharpened results.
- **Figure 26: not affected.** It is the 2024/2025 GMST-ensemble paragraph (see `u250_ensemble_check/REPORT.md`) and doesn't use newsamples.
- **Not affected:** the per-method exceedance probabilities in `evaluate_heads_tails_ensemble.py`, which apply the scale correctly through `pvalue`; the per-method statistics; and PIVW everywhere.
- **Inherited from earlier code:** the same line is in `old_combined75/hist_fitprob75.py:388`, so this goes back to before the 12 Oct rename.

## Proposed fix (needs John's approval)

This is a one-line change in `calibrate_method_uncertainty.py:380`:

```python
# before
newsamples = (result['resample'](years, nsamples,k) - central_est[:, np.newaxis])/best_alter_scale + central_est[:, np.newaxis]
# after
newsamples = (result['resample'](years, nsamples,k) - central_est[:, np.newaxis])*best_alter_scale + central_est[:, np.newaxis]
```

**Cheaper alternative** that avoids rerunning the methods: post-process the existing npy files with `c + (old − c)·s²` on the sampled-method rows only. This is algebraically identical to the fix apart from RNG. It needs a list of which rows are sampled. That list follows from `isinstance(result, dict)`, or it can be checked against rms_z ≈ s² in the log. Note that the future files use per-run s, taken from each run's `_names_var_scale.csv`.

## Rerun list (after the fix)

1. `calibrate_method_uncertainty.py` for the historical run and for all future runs (the Oscar ensemble), **or** the post-processing patch above, to regenerate `Results2/*_newsamples.npy`.
2. `calibrate_historical_combinations.py` ('hh') to get the new historical CSCM nsharp, KL and crossing stats. Any paper text quoting these numbers needs updating.
3. `create_combination_method_LLOffcompute.py`, plus `_percase` if it's used, then `create_combination_method_LLOffcombine.py` for **Table S4**. This can be folded into the pending future leave-one-out rerun.
4. Regenerate any CSCM figures or histograms (`plot_fut_results.py`, stacked variant) from the new Results3.

**Expected effect:** CSCM nsharp drops from about 2 to about 1, CSCM KL improves (in the dashboard test, 119 → 78 for the 7-method set), and the CSCM spread becomes comparable to PIVW. Because `nsharp ≥ 1`, CSCM can't widen a blend that is too narrow. Check that the refit doesn't hit the lower bound of 1 in the future leave-one-out cases.

---

## Update 2026-10-19: John approved the fix

### Done

- **Code fix applied.** `calibrate_method_uncertainty.py:380` now uses `*best_alter_scale`. This is the only main-code change.
- **`Results2/historical_newsamples.npy` patched** with `patch_newsamples_scale.py`, which applies the exact transform c + (old − c)·s² to the 26 sampled-method rows.
  - The patch is idempotent: `Results2/newsamples_scale_patch_manifest.json` records which files are done.
  - Acceptance test: every sampled row now has rms_z ≈ 1.0 and 90 % coverage of about 0.9. eROF is the exception because of its separate offset problem.
  - PIVW `comb_*_scales.npy` is byte-identical before and after.
- **The future files have NOT been patched yet.** The Codex-managed controller (pid 32699, `.../future-final-refresh-2026-10-18`) is mid-run. Its workers read `Results2/fut_*_newsamples.npy` live, and it checks those files' size and mtime against `evaluation-validation.json`.
  - Status: 504/1350 jobs done. All 270 mode-2 jobs are done; these are individual methods and don't use newsamples. Mode 7 has 234/270 done, and those used the buggy samples.
  - Modes 8, 16 and 32 haven't started.

### Two further problems found while rerunning the historical CSCM

These come from an instrumented scratch copy; the main script is unchanged. Logs are in the session scratchpad.

**1. The nsharp optimizer never moves.**
`minimize(rescale_KL_stack, x0=2, bounds=[(1,10)])` returns 2.000 every time. That happens in `calibrate_historical_combinations.py:245`, `create_combination_method_LLOffcompute.py:269`, `_percase:265`, and the controller's frozen kernel. The objective is discontinuous at integer nsharp, because the number of draws jumps, so L-BFGS-B's finite-difference gradient is meaningless.

7-method historical CSCM KL. The KL of PIVW is 79.2 in every case.

| samples | nsharp = 2 (as coded) | grid / bounded-scalar optimum |
|---|---|---|
| buggy | 83–136, depending on weight noise | 132.1 at nsharp 1.17 |
| fixed | 85.6 | **78.6 at nsharp 1.11** |

2-method (anthropogenic) CSCM KL. The KL of PIVW is 87.2.

| samples | nsharp = 2 | bounded optimum |
|---|---|---|
| buggy | 136 | 123 |
| fixed | 156 | **97.5 at nsharp 1.01** |

So the scale fix only helps CSCM if the optimizer is fixed as well. With nsharp stuck at 2, correctly scaled samples get over-sharpened: the mean SD is 0.025, against 0.032 for PIVW.

Proposed fix, in each of the four places above:

```python
best_nsharp = minimize_scalar(lambda v: rescale_KL_stack([v]), bounds=(1, 10), method='bounded', options={'xatol': 0.01}).x
```

This needs `from scipy.optimize import minimize_scalar`.

**2. The historical stacking weights aren't reproducible.**
`infill_log_likelihoods` uses unseeded `np.random.choice`, and `calibrate_historical_combinations.py` sets no seed. With identical inputs, the 7-method weights moved between runs. GWI_tot_CGWL, for example, got 0.0007, 0.0036, 0.004 and 0.49. The future leave-one-out code already seeds (`np.random.seed(base_seed)`), so this affects the historical CSCM only.

Proposed fix: add `np.random.seed(<fixed>)` before `infill_log_likelihoods` in `calibrate_historical_combinations.py`.

### Answer to "Is S4 affected? Anything else?"

- **Table S4's CSCM (`sharpened_blended`) columns: yes**, for every mode that uses sampled methods (7, 8, 16 and 32). PIVW and mode 2 are not affected.
- **Historical CSCM numbers** in the text and supplement.
- **The CSCM histogram or curves** if they come from the sharpened results.
- **The dashboard refit** (`climate-dashboard/source/calibrate_dashboard_combinations.py`) applies an s² rescale itself. That must be removed once the npy files are fixed, or the correction is applied twice.
- **Not affected:** Figure 26, PIVW, per-method probabilities and statistics, and the mode-2 individual reference rows of S4.

### Remaining steps, in order, once John decides

1. Stop the controller, or let it finish if he'd rather.
2. Run `patch_newsamples_scale.py` on all 270 primary `fut_*` sets.
3. If approved, apply the optimizer fix and the historical seed fix.
4. Relaunch the 1,350-case combination run. The controller's validation record has to be refreshed because the newsamples mtimes will have changed.
5. Rerun the historical CSCM for the 7- and 2-method sets.
6. Regenerate Table S4 and the CSCM figures.

---

## Update 2026-10-19 (afternoon): John approved everything; it's done or running

### Code changes (uncommitted)

| file | change |
|---|---|
| `calibrate_method_uncertainty.py:380` | `/best_alter_scale` → `*best_alter_scale` |
| `calibrate_historical_combinations.py` | `minimize_scalar` bounded search for nsharp, plus `np.random.seed(20261019 + nmethods)` before `infill_log_likelihoods` |
| `create_combination_method_LLOffcompute.py`, `create_comb_method_LLOffcompute_percase.py` | the same `minimize_scalar` change (they already seed) |

### Data

- **All 271 `Results2/*_newsamples.npy` patched** (historical plus the 270 future sets). The manifest is `Results2/newsamples_scale_patch_manifest.json` and the log is `patch_future_newsamples.log`.
- Spot checks on three future runs give rms_z ≈ 1 for every sampled method.
- The future method set has 18 sampled methods; historical has 26.

### Historical CSCM rerun

- Official stdout: `hist_cscm_rerun_official.log`.
- Instrumented stats: `hist_final_nm7_a.log` and `hist_final_nm2_a.log`.
- Stacking weights are now identical across repeat runs.

| set | before (buggy samples, nsharp stuck at 2) | after (fixed samples + optimizer) | PIVW |
|---|---|---|---|
| 7-method | nsharp 2.00, KL 83–136 (weight noise), SD 0.032–0.11 | **nsharp 1.07, KL 80.0, SD 0.032** | KL 79.2, scale 0.0323 |
| 2-method (anthropogenic) | nsharp 2.00, KL 136, SD 0.085 | **nsharp 1.005, KL 98.6, SD 0.037** | KL 87.2, scale 0.0354 |

Crossing-year statistics hardly move: 7-method all-threshold RMSE 7.32 → 7.36 yr, and the 1.0 °C error stays about 0.8 yr.

### Future / Table S4

- The Codex controller (pid 32699) was stopped at 522/1350 jobs. A marker, `STOPPED_BY_CLAUDE_newsamples_fix.txt`, is in its output folder. Its results used unpatched samples, so its CSCM values are invalid.
- **The relaunch is running detached:**
  - script: `~/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/work/future_refresh_newsamples_fix_2026_10_19/run_all.sh`
  - output: `outputs/future-refresh-newsamples-fix-2026-10-19/` (`run_all.log`, `combination-run/progress.json`)
  - It uses the same pipeline. The readiness check now accepts the patched newsamples via the manifest, and every other input must still be unchanged.
  - Equivalence check: PASS.
  - Mode order is 7, 16, 32, 8, 2, so the S4 modes finish first. It aggregates automatically to `combination-run/table-S4.md` and its CSVs when done.
- Throughput is about 26 s per job with 2 workers, so roughly 10–20 h in total.
- Early mode-7 cases: the PIVW columns are identical to the stopped run, as expected. The CSCM KL moves by a few units per case, in both directions.

### Still to do

- When the run completes, compare `table-S4.md` with the current supplement and update the CSCM text and numbers.
- Remove the s² workaround from the dashboard refit (`climate-dashboard/source/calibrate_dashboard_combinations.py`).
- The eROF offset issue noted above.

---

## Update 2026-10-19: eROF resample offset fixed (approved)

### The problem

`ppf_resample` in `eROF_tot_method.py` and `eROF_anthro_method.py` didn't force `k=1`, unlike the mean, SE, p-value and likelihood functions. With `k=0`, it drew from `cdists[yr-1950]`. But `cdists` is `rdists`, indexed from 1850, so the draws came from the year 100 years earlier.

### What changed

- Added `k=1` at the top of `ppf_resample` in both files.
- Regenerated only the eROF rows (51 and 58) of `Results2/historical_newsamples.npy` with `regen_eROF_newsamples.py`.
  - The script asserts that the central estimate is identical to the stored one.
  - It uses the stored scale and the corrected c + s·(r − c) form, with seed 20261019.
- Result: eROF_anthro now has rms_z 1.00 and coverage 0.88; eROF_tot has rms_z 1.00 and coverage 0.91.

### Impact: historical only

- eROF `run_method` returns NaN for any non-historical experiment, so eROF isn't in any of the 270 future method tables (60 methods each). It therefore can't enter any Table S4 section, including the largest selection (mode 32, 33 methods).
- It isn't in the historical 7- or 2-method CSCM either.
- Nothing else reads the eROF sample rows.
- The running S4 job is unaffected.

---

## Update 2026-10-19: dashboard refit corrected (approved)

### What changed in `climate-dashboard/source/calibrate_dashboard_combinations.py`

- Removed the `scale**2` rescale. The stored samples are now already correct, so keeping it would apply the correction twice.
- Removed the split into separate "stored" and "corrected" fits.
- nsharp now uses the same bounded `minimize_scalar` search as the main scripts.
- The infill seed is now 20261019 + 7, the same as `calibrate_historical_combinations.py`.

### Rerun results (`data/supplementary/dashboard_combination_calibration.json`)

- The 7-method CSCM weights now exactly reproduce the paper script: FaIR_comb_unB 0.697, EBMKF_ta4 0.256, GWI_tot_CGWL 0.009, CGWL10y_sfUKCP 0.037.
- CSCM: nsharp 1.088, KL 80.0, mean SD 0.032.
- PIVW: scale 0.0323, KL 79.2. These are unchanged and still pass the `comb_real_scales.npy` assertion.

### Before (previous JSON)

- Weights: FaIR 0.475, EBMKF 0.20, GWI 0.32.
- nsharp 1.031, KL 78.4.
- The weight difference comes only from the different infill seed, which shows how sensitive the stacking weights are to the random infilling. Note that they collapse onto three forcing-based methods in both cases.

### Downstream

`source/export_dashboard_climate_metrics.py` reads `cscm_stack_weights` and `cscm_nsharp`; both keys are kept. The exported dashboard data has **not** been regenerated yet.

---

## Update 2026-10-19: dashboard export rerun (approved)

`climate-dashboard/source/export_dashboard_climate_metrics.py all` ran with exit 0.

### Changes

A full diff against the pre-run exports, which are snapshotted in the scratchpad, shows that **only `CSCM_comb`** changed: its paths and `supplementary_calibration`. Nothing changed in the PIVW, individual-method or anthropogenic exports.

rpi_real CSCM:

| | before | after | PIVW |
|---|---|---|---|
| 2024 level | 1.357 [1.226, 1.520] | 1.369 [1.242, 1.535] | 1.366 [1.230, 1.532] |
| mean SE, 1975–2024 | 0.0384 | 0.0328 | 0.0323 |

- Crossing comparisons move by about 1 year or less.
- At the 1.5 °C comparison, the number of members that crossed rose from 2 to 7.

### Not yet done

`render_dashboard_figures.py` and `assemble_dashboard_site.py` have not been run, so the site and figures still show the old CSCM.

---

## Update 2026-10-19: dashboard figures and site rebuilt (approved)

`render_dashboard_figures.py` and `assemble_dashboard_site.py` both ran with exit 0 (96 variants).

- `site/dist/assets/data/rpi_real.json` is byte-identical to the new export (nsharp 1.088).
- The CSCM figure shows 1.37 [1.24, 1.54] °C for 2025, with 7/100 members crossing 1.5 °C.
- `climate-dashboard-website.zip` and `climate-dashboard-figures.zip` (dated 11 Oct) are **not** rebuilt by these scripts and are stale.

---

## Update 2026-10-19 18:25: S4 run resumed

- The Mac rebooted at 17:55, and the last result was written at 17:14, so the run stopped with 309/1350 cases done: all of mode 7, and 39 of mode 16.
- It was resumed with `work/future_refresh_newsamples_fix_2026_10_19/resume_run.sh`. That script runs the `run` step only, skipping verified cases, under `caffeinate -i`.
- Measured times, wall clock with 2 workers: mode 7 about 25 s per job, mode 16 (18 methods) about 131 s per job.
- Estimate: about 8 h for the rest of mode 16, plus about 15–25 h for mode 32 (33 methods, extrapolated), plus about 2–3 h for modes 8 and 2. That makes roughly 25–35 h in total, finishing late on Tue 20 Oct to midday Wed 21 Oct.
- To get Table S4 earlier, run `aggregate --modes 7 16 32` once mode 32 finishes.

---

## Update 2026-10-19 (evening): PIVW window made consistent (approved)

- **The issue:** `calibrate_method_uncertainty.py:363` computes `err_var100` over `[lhund:]` with `lhund=-101` on the 251-year MPI series. That gives 2000–2090 (91 finite years), while the predictive SD and KL use 2000–2089. NorESM already used 2000–2089 for both.
- **The fix:** `create_combination_method_LLOffcompute.py` now recomputes the PIVW error variances inside the leave-one-out over `[lhund:-10]`, i.e. 2000–2089, from `central_estfut` and `standardsf`. The Oscar worker was regenerated from it.
- **Checks:**
  - The NorESM recomputed variances equal the stored ones to 1e-13; the MPI ones change by up to 6%.
  - The effect on pooled PIVW weights is at most 0.0006.
  - On two test cases, the CSCM outputs are bit-identical and PIVW moves by about 0.03 in KL.
- **Still to do:** the Oscar array is to be cancelled and resubmitted with the new worker, and `check_and_aggregate.py` now expects exact CSCM agreement and small PIVW differences.

---

## Update 2026-10-20: Oscar Table S4 run complete

### Run

- Oscar array 7637718: 1080/1080 cases present and finite. Results are in `Results3_fix_oscar/`.
- Inputs were hash-identical to the Mac. The window fix (2000–2089) was in.
- Task 749 (mode 32, NorESM Volc member 59) hit the SLSQP 100-iteration limit, but its weights had been stable since about iteration 50.

### Cross-check against the 321 local cases

- Differences are cross-platform floating point only.
- CSCM KL: median difference 7e-5, max 0.14 (0.34%).
- RMSE agrees to about 1e-9 (max 6e-5).
- Crossing counts are identical. Only 2 of 321 cases moved the first-crossing time, by one month.

### Outputs

`oscar_loo_fix/table_S4_out/` holds `table-S4.md`, the numeric and formatted CSVs, `final_combined_{7,16,32,8}methods_wide_ff.csv` and `comparison_old_vs_new.csv` (against the Dec 2025 `final_combined_*`).

### Headline

- **Old:** CSCM KL was worse than PIVW in all 15 mode × ensemble cells.
- **New:** CSCM is better in 9/15 (SSP1-2.6, NorESM Volc and VolcConst, at every method count), essentially tied for SSP2-4.5, and still slightly worse for SSP3-7.0 (by +0.02 to +0.05).
- CSCM KL fell by 10–31%; RMSE barely changed.
- The old values also predate the EBMKF and other corrections, so not all of the change is due to the CSCM fixes.

---

## Update 2026-10-20: reorganized into fut_LOO/

The scripts, the newsamples patch, the 1,080 Oscar results and the Table S4 outputs were moved to the top-level `fut_LOO/` folder. Legacy entry points were moved to `fut_LOO/legacy/`, with their paths updated. The root README points to the new folder.

---

## Update 2026-10-20: held-out scoring limited to 2000–2089 (approved)

- MPI held-out KL and RMSE had included 2090: 91/66 valid years against the 90/65 normalization.
- The kernel now scores KL and RMSE over 2000–2089 using `standard_eval`; crossing metrics are unchanged.
- Tested locally: NorESM results are bit-identical; for MPI, KL and RMSE shift slightly and crossings are unchanged.
- Next: rerun the 600 MPI tasks on Oscar.

---

## Update 2026-10-20 (late): Table S4 rebuilt with 2000–2089 held-out scoring

### Run

- On Oscar (array 7668784), 593 of the 600 MPI tasks finished.
- Seven 33-method MPI cases hit the 6 h limit: SSP1-2.6 members 14 and 33, SSP2-4.5 members 25, 26 and 47, SSP3-7.0 members 11 and 22. They were computed on the Mac (`fut_LOO/results/mac_reruns_2089eval/`) and copied into `Results3_fix_oscar/`.
- All 1,080 cases now use 90/65 valid years.

### Table

- `fut_LOO/output/` was regenerated, with headers KL/RMSE 2000–2089 and 2025–2089.
- MPI KL fell by about 1–2%; NorESM is unchanged.
- CSCM still beats PIVW in 9/15 cells, ties for SSP2-4.5, and is worse for SSP3-7.0.
