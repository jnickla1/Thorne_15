# Climate dashboard correction log

> **Latest decision:** preserve the current pathway. Change only the retrospective lag to the full cut2021 all-forcing trajectory. The previous proposal to change current inputs is superseded; its export script is not needed. The latest entry below records validation with all 101 transferred cut2021 matrices and completion of the figures and website. Earlier numerical entries describe intermediate revisions.

## 2026-10-11 — In-place method corrections

The existing method files in `Thorne_15_codefigurestats/Methods` were edited in place at the user's request. No additional method files were created there. The current scientific outputs are unchanged in the regression checks described below. Backups, diffs, corrected result copies, and validation tables are stored outside `Methods` in the dashboard workspace.

### Applied changes

| Change | Current method outputs altered? | Retrospective or derived products altered? |
|---|---|---|
| Align the 2025 historical/heads-tails retrospective FaIR extraction with the Oscar following-boundary convention in all four `_unB` methods | No | Yes: main source selection begins at row 101 rather than 100; 2025 selects boundary 2026 |
| Apply the `FaIR_comb_unB` biases separately to its retrospective non-natural, all-forcing, and lagged all-forcing components, matching the current correction | No | Yes: net bias changes from −0.0086667°C to +0.0026625°C |
| Give `FaIR_comb_unB` retrospective mixing the current coefficients 0.6125, 0.25, 0.1375 and its four-year lag | No | Yes: weights now sum to 1 rather than 7/6; five-year lag replaced by four-year lag |
| Remove obsolete, unreachable correction alternatives from `FaIR_comb_unB` | No | No additional numerical effect beyond the executed changes above |
| Make `GWI_tot_CGWL` expose its current output as a temporary retrospective proxy | No | Yes: its previously absent retrospective branch now duplicates current means, SEs, probabilities, likelihoods, and resampling |
| Apply the same explicit proxy to `GWI_anthro_CGWL` for the anthropogenic option | No | Yes: same temporary proxy policy |
| Rebuild corrected retrospective result copies and dashboard data | No change to current method curves or current combinations | Yes: retrospective means/rates and the rate headline are refreshed in the copies |

The FaIR row alignment change is restricted to the verified 2025 historical and heads/tails branches. Older and future exports retain their existing boundary selection. The combined retrospective coefficients, bias handling, and four-year lag are corrected wherever that method's executed recipe is used. Full retrospective sample spread is retained; current uncertainty contraction is untouched.

### Why the weights were wrong

The error is inside `FaIR_comb_unB`, not in the seven-method PIVW combination. Its older retrospective expression combined 2/3 of the non-natural component with 1/3 of the all-forcing component, already totaling one, then added another 1/6 lagged all-forcing term without reducing the other weights. The current branch had a different, normalized recipe. The seven-method inverse-variance combination already divides by the weights of available methods; its weights were not the source of the 7/6 problem.

`GWI_tot_CGWL` is the missing retrospective member among the seven realized methods. Its reader explicitly implemented only `k=0`; the `rdists` list was never populated or used to construct a full-window retrospective series. The supplied CGWL CSV supplies successive cutoff estimates for 1950–2025. This is an intentionally absent reader branch, not a failed load in the existing branch. The original rate script already used current GWI output as its fallback; making the proxy explicit preserves that policy. The anthropogenic counterpart also lacked a retrospective branch.

### Current-output verification

- All four FaIR methods, 100 members, and both calibration options: **800 cases**, with exact equality of complete current sample arrays, means, and SEs before and after the edits.
- Current FaIR threshold probabilities were compared in both one-sided and two-sided forms. Current endpoint likelihoods and seeded resamples also matched exactly.
- Both GWI source APIs were compared with their backups for a representative member under each calibration; current results matched exactly, and the new retrospective APIs matched current outputs. Their cached current arrays were preserved for every member.
- All selected current method mean/SE arrays and both current combination curves were checked for all **400** result files.
- PIVW and equal-weight sums were verified at every member/year where outputs exist. Both current and retrospective available-method weights sum to one.
- Maximum current mean change: **0.000000°C**. Maximum current SE change: **0.000000°C**.
- Original `Results` pickle files remain byte-for-byte unchanged. The source edits match the tested staged files by SHA-256.

No hour-long rerun is triggered by changed current outputs in the tested 2025 heads/tails cases. Retrospective products still need updated caches; copies have already been rebuilt. The corrected cache directory is `outputs/climate-dashboard/data/corrected_results` in the dashboard workspace. Existing scripts reading the original project's `Results` directory still see the old cached retrospective values until those results are adopted or regenerated.

### Endpoint changes, mean across 100 members

| Calibration | Method | Current before (°C) | Current after (°C) | Retrospective before (°C) | Retrospective after (°C) |
|---|---|---:|---:|---:|---:|
| rpi | `FaIR_all_unB` | 1.440269 | 1.440269 | 1.420339 | 1.440274 |
| rpi | `FaIR_anthro_unB` | 1.348108 | 1.348108 | 1.323742 | 1.348113 |
| rpi | `FaIR_comb_unB` | 1.367696 | 1.367696 | 1.551741 | 1.365381 |
| rpi | `FaIR_nonat_unB` | 1.337859 | 1.337859 | 1.315439 | 1.337865 |
| satcal | `FaIR_all_unB` | 0.745915 | 0.745915 | 0.725984 | 0.745920 |
| satcal | `FaIR_anthro_unB` | 0.653753 | 0.653753 | 0.629387 | 0.653759 |
| satcal | `FaIR_comb_unB` | 0.673341 | 0.673341 | 0.741661 | 0.671027 |
| satcal | `FaIR_nonat_unB` | 0.643505 | 0.643505 | 0.621084 | 0.643510 |

For the preindustrial realized option, the corrected dashboard headline reads the unchanged current PIVW level, **1.364154°C**, and uses the corrected retrospective history for its recent rate, **0.318513°C/decade**, combined SE **0.054385°C/decade**. The old saved headline level, 1.408236°C, came from retrospective-selected endpoints, and its old rate was 0.335652°C/decade. Thus the reported headline changes even though the scientific current method outputs do not. The regenerated headline agrees with the minimally corrected legacy rate script; its uncertainty includes both between-method and between-member dispersion.

### Earlier issues and remaining distinctions

1. **The 1.0°C definitive rate/uncertainty inconsistency.** The original `histens_rates.py` retained an endpoint-difference rate after calculating an OLS slope, but replaced its uncertainty with the OLS slope uncertainty. The dashboard working copy assigns `rate10a = slope10 * 10`, matching the 0.5°C convention. This changes the definitive rate summary, not a current method curve. The original rate script remains untouched as originally requested.
2. **Headline level versus rate history.** The original rate script used the final retrospective-selected value as its headline level. The corrected dashboard separates the unchanged current level from the retrospective 2011–2025 rate history. This is a reporting correction.
3. **CGWL endpoint definition.** `CGWL10y_sfUKCP` currently averages 10 past observations and 10 forecast years with equal half-weights. Its retrospective taper uses 11 observations and 10 forecast years, divided by 21. The preindustrial endpoint means are 1.397890°C current and 1.386392°C retrospective. This separate difference remains unchanged, pending a decision about the intended retrospective window.
4. **FaIR's remaining lag-input distinction.** The current 2025 mixture takes its lagged input from the combined array's 2021 row, which the Oscar slicer obtained from the `cut2021` full trajectory. The retrospective 2025 mixture takes its 2021 input from the `cut2025` trajectory. The trajectories were both created first on Oscar, but with different observation cutoffs/posterior constraints. The difference appears when they are mixed locally. After the applied corrections, FaIR's retrospective mean is 1.365381°C versus current 1.367696°C. The remaining approximately −0.0023°C difference is mainly this lag-input distinction, with small sample/noise effects.

### Proposed next correction — not applied

For the interpretation that every current estimate is computed from the full record available at that cutoff, obtain both the contemporaneous and four-year-lag all-forcing values from the **same cutoff trajectory**:

```python
# Inside the Oscar reslicer, for the trajectory constrained through `year`:
row_index = year - 1750 + 1
current_all = arr[row_index, :]
current_all_lag4 = arr[row_index - 4, :]
```

The local current mixture would consume this separately exported lagged series instead of `combined_all[:-4]`. This also keeps posterior configurations aligned within each cutoff run. Existing full Oscar trajectories can be resliced; new FaIR simulations should not be needed if those files remain available. The complete corrected current history cannot be reconstructed here from only the final 2025 full trajectory.

**Superseded by the latest instruction to preserve the current branch. This proposal was not applied.** The alternative is to retain the operational mixture of past cutoff estimates and define the retrospective endpoint to use that same archived lag input; this preserves current outputs but changes what the retrospective estimator means. No arbitrary endpoint overwrite or extra offset has been applied to force agreement.

### Audit and recovery files

The dashboard workspace contains `review/inplace_changes_manifest.json` with hashes and backup locations, `review/backups/` with originals, and `review/patches/` with exact diffs. The full regression tables are `review/corrections_current_output_regression.csv`, `review/corrections_endpoint_summary.csv`, and `review/corrections_regression_report.json`. Each corrected result identifies a current proxy explicitly, rather than implying it is a full-window retrospective estimate.

Future corrections should append an entry stating the scientific change, current-output comparison, affected caches/summaries, and any rerun required.

## 2026-10-11 — Superseded proposal: matched current cutoff and lag inputs

The user accepted recalculating current results if needed. Local inspection confirms that no `cut2021` trajectories are present in `resliced_headstails25`: their extracted rows are bundled into `combined_all-headstails25r{run}_all.npy`. The complete `cut2025` trajectories are present in its `retrospective/` subdirectory.

Prepared `source/export_fair_same_cut_lag4.py` for Oscar. It exports both the current and four-year-lagged all-forcing row from the same cutoff trajectory, retaining column order and NaN padding in one compressed NPZ. The current rows are included so the transferred inputs can be compared with the existing combined arrays before applying the local correction. Missing cutoffs cause an error rather than silently shortening the series.

Validation used a local fixture based on the real r1 cut2025 trajectory: 2025 current boundary row 276 and lag row 272 were retained exactly, as were the adjacent 2024 rows, column order, and NaNs. Missing-cutoff rejection also passed. This checks the export logic; the original multi-cutoff files still need to be read on Oscar.

**Current-output effect:** none yet. Export preparation and documentation do not change method outputs. This proposal was subsequently superseded by the retrospective-only correction below; its NPZ transfer is no longer required.

## 2026-10-11 — Applied: retrospective lag from the full cut2021 fit

**Selected behavior:** preserve the existing current pathway. This supersedes the proposed current-input change and the per-cutoff NPZ export.

Updated the existing `Methods/43_Forcing_Based/1_ERF_FaIR/FaIR_comb_unB_method.py` in place. For the 2025 historical and heads/tails branches, it now loads `resliced_headstails25/retrospective/all-headstails25r{run}cut2021_temp_all.npy` and uses `[97:273]` for the four-year-lagged retrospective term. Contemporary retrospective components still come from cut2025. Historical uses r0; ensemble members use r1–r100. No additional method modules were created.

The added matrix receives the satellite calibration offset once and the same all-forcing bias subtraction as the other all-forcing input. Full retrospective spread is retained. Where cutoff fits have different sample counts, source columns are padded or cropped to the final retrospective width, preserving column order and the sample prefix used by the existing current mixture. Preserving retrospective width also preserves the subsequent random draw sizes. Older and future export branches retain their previous lag input.

**Current-output effect: none.** The current loading, transformations, mixing expressions, and random perturbation logic are unchanged. Across 200 member/calibration cases, complete current sample arrays, means, and SEs matched the immediately preceding version exactly. The historical r0 case also matched. The expected retrospective change was checked sample by sample, including the satellite calibration; raw-width padding and truncation were exercised. Missing-file and invalid-shape guards passed. Tests used reader fixtures because the real cut2021 trajectories are not yet present locally; no synthetic inputs were saved into the source data folders or scientific result caches.

**Data and rerun status:** the code is applied and ready for the additional data. It reports a clear missing-input error until the cut2021 files are transferred. Actual cut2021 endpoint consistency and the resulting retrospective curves/rates remain to be validated. Existing corrected result copies and dashboard rate summaries still represent the preceding cut2025-lag revision and must be refreshed after transfer. No current-output cache rebuild is required for this retrospective-only change.

Copy only the full all-forcing cut2021 trajectories (r0–r100), retaining the existing filenames. Run on the Mac:

```sh
scp 'jnickla1@ssh.ccv.brown.edu:/users/jnickla1/data/jnickla1/fair-calibrate/output/all-headstails25/all-headstails25r*cut2021_temp_all.npy' '/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats/Methods/43_Forcing_Based/1_ERF_FaIR/resliced_headstails25/retrospective/'
```

No new Oscar simulations or reslicing are required. Neither nonat nor anthro cut2021 matrices are needed. Validation details: `review/retrospective_cut2021_regression_report.json` and `review/retrospective_cut2021_current_guard.csv`. The immediately preceding method is backed up outside Methods as `review/backups/Methods/43_Forcing_Based/1_ERF_FaIR/FaIR_comb_unB_method.py.before_cut2021`; the original pre-correction backup remains intact. The incremental diff is `review/patches/FaIR_comb_unB_retrospective_cut2021.patch`.


## 2026-10-11 — Verified transferred cut2021 files; dashboard completed

The user transferred all 101 full all-forcing cut2021 trajectories (r0–r100) and changed the combined method's input folder to `resliced_headstails25/retrospective_4/`. This user-selected path is retained. Each trajectory's raw row 272 exactly matches the valid column prefix and NaN padding of the corresponding combined array's 2021 row. The transfer is complete; the previous pending-transfer statements are historical.

**Current-output effect: none.** With the actual transferred matrices, all 800 FaIR method/member/calibration cases pass exact current sample, mean, SE, threshold-probability, likelihood, and seeded-resample comparisons against the original backups. All selected current arrays and both current combinations in 400 result files remain unchanged. Current and retrospective available-method weights normalize to one. The original `Results` caches remain byte-for-byte unchanged.

The corrected combined FaIR retrospective preindustrial endpoint is now **1.367686220°C**, compared with **1.367695715°C** current. The difference is **−0.000009495°C**; the earlier approximately −0.0023°C lag-input discrepancy is resolved. Satellite-calibrated endpoints have the same residual. Remaining sample-spread and noise definitions are preserved. The separate CGWL window-definition distinction remains unchanged.

Rebuilt all **400** corrected retrospective cache copies and all four dashboard datasets. The final preindustrial realized PIVW headline is **1.364154°C** current, with corrected retrospective rate **0.318823°C/decade**, SE **0.054156°C/decade**, and 95% rate interval **0.212677–0.424970°C/decade**. The 2025 saved calibrated probability above 1.5°C is **5.871641%**, classified as very unlikely. These replace the preceding intermediate rate values in this log. No current-pipeline rerun is required for the tested cases; original-cache users must still adopt or regenerate retrospective products.

Generated the main 16:9 mockup in PNG/PDF/SVG and **52** method/calibration/definition/background variants. Definitive crossings use centered 20-year means; their rates use centered 10-year means. The seven requested quartic estimates are flagged in the member data. Light-blue comparison segments extend 15 years before the method-predicted crossing, with exact rate labels below the dark definitive labels. A dark-blue 15-year segment ends at 2025.

The three bottom selected-method density panels reproduce the paper's saved current probabilities and smoothing convention. Those exact probabilities also populate their original assessment windows in the heatmap. Gaussian extensions outside those windows and for the fourth threshold are disclosed in the methodology.

The portable static website supports all 52 variants, interactive zoom/hover plots, image switching, downloadable current-method/summary/member CSVs and complete JSON, and verified method-reference hyperlinks. A probability strip supplies the optional uncertainty indicator. Source originals remain untouched; dashboard copies are under `source/`. Hosting/reproduction instructions and scientific limits are in `README.md`.

Evidence: `review/cut2021_transfer_validation.csv`, `review/corrections_regression_report.json`, `review/corrections_current_output_regression.csv`, `review/corrections_endpoint_summary.csv`, `review/combination_weight_sums.csv`, and updated source hashes/diffs in `review/inplace_changes_manifest.json` and `review/patches/`.


## 2026-10-11 — Display adjustments requested during preview

The top-left warming plot now starts at 1960 in the main figure, all 52 variants, and the interactive website. All three crossing-distribution panels use a fixed **0–0.8 probability/year** y-axis across every method/calibration/definition choice. The largest bar fraction in the complete dataset is 0.72, so the shared range includes all bars and saved density curves. Interactive distribution panels retain x-axis zoom; their y-axis is fixed and the autoscale control is omitted to preserve the comparison. These display adjustments change no scientific values or current outputs.


## 2026-10-11 — Red method estimates, paper histogram palette, and explicit density magnification

The method comparison segments, 2025 method line, current headline values, and smoothed density curves now use red; observational benchmarks retain blue. Benchmark and Method annotation headings are on separate lines from their estimates, with a larger vertical separation from the plotted lines.

The histogram palette and opacity were extracted directly from `hist_evaluation_script.gen_color` and `plot_heads_tails_probs.build_method_colors_for_goal`. PIVW retains the paper's dark grey and gap width; individual methods retain their group colors. The palette is saved in `data/paper_histogram_colors.json` and carried as display metadata in the JSON exports.

The shared distribution ceiling is now **0.3 probability/year**, superseding the preceding 0.8 choice. Bars exceeding that ceiling carry `//` and a label with their full height; hover values also retain the original fractions. The website's Magnify switch explicitly changes only the rightmost threshold to a finer labeled y-scale without renormalization. The Background lines choices are renamed Raw heads/tails ensemble and Smoothed ensemble (method); toggling them only replaces the grey/black ensemble curves. These are display changes and alter no scientific values or current outputs.


### 2026-10-11 — Current exceedance wording

Updated the website probability-card heading to “Probability of currently exceeding [threshold]°C relative to [reference period]” for every selected option. Updated the initial HTML label and the JavaScript version query to avoid stale browser caches. Applied in place to the Downloads GitHub checkout, the secondary dashboard export, and the original workspace export. Wording only: no data, scientific values, or current outputs changed. GitHub publication is left to the user.


### 2026-10-12 — Current warming levels in the rate report

Applied the user-authorized in-place reporting correction to `histens_rates.py`: 2025 individual, equal-weight, and PIVW warming levels now use current endpoints, and `year_cross_sd` uses their member spread. Retrospective histories and existing current proxies still supply the 15-year rate fits. A separate console table compares current and raw retrospective endpoints for each reported method, with means, paired differences, maximum differences, matching counts at absolute tolerance 1e-6°C, and missing endpoints. Combination comparisons retain the same current-proxy fallback as the rate fits. Diagnostics do not add CSV columns and do not assert that every method's endpoints must coincide.

Validated the actual reporting/aggregation code on 100 fixed member snapshots for each of all four experiments. All rate, rate spread, uncertainty, and window columns were exactly unchanged; current levels and their member spreads matched independent calculations; the eight-column CSV schema was unchanged. Checked mismatch/missing/no-pair diagnostics. The dashboard exporter now accepts current-level CSVs and retained legacy retrospective-level CSVs, recording `saved_summary.level_source`; 52 format cases passed and invalid saved levels remain rejected. Matching README text was updated; `README_old.md` was untouched. No source Results pickles or research CSVs were overwritten by this validation; the user's rerun workers remain responsible for regenerating them. Actual method current temperature arrays and the website headline convention did not change.

## 2026-10-12 — Clarify workflow script names

Renamed the 12 approved research/dashboard entry points and updated Python imports, active shell wrappers, and both READMEs. The root README includes an old-to-new migration table and updated workflow commands. This is a naming-only change: method identifiers, scientific calculations, data/result filenames, and website exports are unchanged. Archived scripts/review snapshots, README_old.md, and upstream submodules are preserved.

## 2026-10-12 — Coauthor readability suggestions

Updated both website copies with clearer control labels and grouped method choices, two-decimal headline temperatures, the 2011–2025 endpoint-rate explanation, more prominent year/definition-specific exceedance wording and IPCC category ranges, a key and subtle shading for benchmark crossing windows, and heatmap guidance/axis/legend labels. Scientific outputs and publication figures are unchanged. Larger ensemble/rate controls, early anthropogenic treatment, and the inherited density time-placement convention remain pending approval.

## 2026-10-12 — Rate/level controls, display period and density alignment

Added independent benchmark/method rate switches that flatten segments at their warming levels rather than hiding them. Kept individual ensemble lines and omitted ribbons following the author’s preference. Temperature plots and downloadable figures now start in 1965. Public CSV downloads preserve earlier rows and carry leading # comments explaining the early anthropogenic limitations; analysis CSV headers remain unchanged. Corrected the inherited annual-density time placement in the paper renderer and dashboard exporter: P(y)-P(y-1) belongs to (y-1,y], moving plotted curves one year left while preserving heights, histograms and exceedance probabilities. Regenerated both paper figures and the dashboard library. Static likelihood wording now also includes extremely unlikely.

## 2026-10-12 — Crossing markers, range hover and label spacing

Separated temperature-plot year tick labels from the bottom y-axis tick, moved heatmap benchmark years into the gaps above probability rows, and added hover explanations for grey benchmark 95% crossing windows. Added red downward method crossing ticks and blue upward benchmark crossing ticks at their respective crossing years, in both interactive and static plots. Rebuilt the figure library, also bringing static headline temperature precision to two decimals. Scientific estimates/data are unchanged. Kennedy’s original suggestions were checked individually; ribbons remain omitted by author preference, whole-percentage rounding is rejected, and a separate rates-over-time plot remains deferred.

## 2026-10-12 — Selective member ribbon and matching time axes

Added a dismissible default-view hint and a one-member pale-yellow approximate 95% calibrated uncertainty ribbon in selected-method smoothed mode. Hover is temporary; tap/member selection can retain the ribbon; Clear selection dismisses it. Uncertainty is exported separately, lazy-loaded and checked against the existing dataset’s source SHA/year grid. Added the observational retrospective centred 20-year background (ending in 2015), with static figures for all three backgrounds (78 variants). Matched the warming and exceedance panels to major ticks 1965, 1975, …, 2025 and minor ticks 1970, 1980, …, 2020, using the same horizontal range. Existing scientific outputs/headline probabilities are unchanged.

## 2026-10-12 — Optional rates-over-time overlay

Added an interactive-only overlay opened from the warming-rate card, with lazy-loaded code/data. Blue uses the existing centred-mean OLS benchmark (ending in 2015), solid red uses retrospective trailing 15-year quadratic endpoint rates matching the 2025 headline, and an optional faint dotted red series uses annual-current fits. The retrospective line cannot be toggled away. Added four small JSON exports and 26 explanatory CSV downloads, regenerated automatically during website assembly from the existing validated JSON. Existing scientific exports, headline estimates, probability calculations and static figures are unchanged. Numerical endpoint checks pass for all 26 options (maximum difference <7e-8°C/decade from JSON rounding).


## 2026-10-14 — Future LOO correction and recalculation

Validated the complete Oscar calibration download and applied approved future-only year alignment, held-out-member exclusion and KL-window corrections. Replaced the legacy December-2025 aggregation date check with exact expected member coverage. Full future combination recalculation is in progress; final numbers remain pending verification. Historical/heads-tails/dashboard outputs and website publication are unchanged by this stage.
