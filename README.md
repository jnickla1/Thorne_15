# Thorne_15: long-term global warming methods and climate dashboard

This is the academic working directory for evaluating methods that estimate long-term global warming and when warming thresholds are reached, for the paper *How will we know when 1.5°C of global warming has been reached?* See the [ESSD preprint](https://essd.copernicus.org/preprints/essd-2025-825/) for the methods and scientific context.

The current observational update runs through **2025**, using **HadCRUT5.1**. The directory contains historical evaluations, model-based future evaluations, heads/tails ensemble evaluations, intermediate caches, paper figures, and a static climate dashboard. Older data versions and experiment outputs are retained because some methods and experiments still use them. This is a collection of research workflows rather than an installed Python package.

**Run the research commands from the repository root**, where `config.py` and the top-level evaluation scripts are located. On the author's Mac this is `/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats`; a fresh GitHub clone may instead be named `Thorne_15`.

- [Refresh the current heads/tails results and dashboard](#refresh-the-current-headstails-results-and-dashboard)
- [Installation and external data](#installation-and-external-data)
- [What is in this directory](#what-is-in-this-directory)
- [Add another observational year](#add-another-observational-year)
- [Historical calibration and evaluation](#historical-calibration-and-evaluation)
- [Future simulations and combination experiments](#future-simulations-and-combination-experiments)
- [Interpretation, caches, and troubleshooting](#interpretation-caches-and-troubleshooting)

## Script names and workflow order

The entry-point scripts were renamed on 12 October 2026 to make their roles clearer. Update saved commands, notebooks, and job submissions using the table below. Repository imports and batch wrappers use the new names; the old entry-point filenames are no longer provided. Method identifiers, input files, result filenames, and calculations are unchanged.

| Previous script | Current script |
|---|---|
| `hist_fitprob.py` | `calibrate_method_uncertainty.py` |
| `create_combination_method_HHonly.py` | `calibrate_historical_combinations.py` |
| `hist_evaluation_script.py` | `evaluate_historical_observations.py` |
| `hist_heads_tails_evaluation_script.py` | `evaluate_heads_tails_ensemble.py` |
| `histens_rates.py` | `export_heads_tails_rate_summary.py` |
| `plot_heads_tails_probs.py` | `plot_heads_tails_paper_probabilities.py` |
| `climate-dashboard/source/histens_rates_dashboard.py` | `climate-dashboard/source/export_dashboard_climate_metrics.py` |
| `climate-dashboard/source/plot_heads_tails_probs_dashboard.py` | `climate-dashboard/source/render_dashboard_figures.py` |
| `climate-dashboard/source/prepare_website.py` | `climate-dashboard/source/assemble_dashboard_site.py` |
| `climate-dashboard/source/rebuild_corrected_retrospectives.py` | `climate-dashboard/source/audit_and_rebuild_retrospective_caches.py` |
| `fut_evaluation_script.py` | `evaluate_future_scenarios.py` |
| `fut_evaluation_script_elim.py` | `evaluate_future_scenario_subsets.py` |

For an existing heads/tails dataset, the order is: `evaluate_heads_tails_ensemble.py` (normally through `run_histens_local.sh`), `plot_heads_tails_paper_probabilities.py` and `export_heads_tails_rate_summary.py`, then the dashboard exporter, figure renderer, and site assembler. Historical recalibration starts with `calibrate_method_uncertainty.py` and `calibrate_historical_combinations.py`; see the detailed workflows below. The optional retrospective-cache audit is not required for an ordinary refresh.

## Refresh the current heads/tails results and dashboard

Use this workflow when the required inputs and calibration files already exist and you want to incorporate a method correction, refresh result pickles, or rebuild the dashboard for the existing 2025 dataset. An annual data update also needs the input and calibration work described below.

### 1. Regenerate the heads/tails member results

Activate the Mac environment and run the batches from the repository root:

```sh
conda activate cleanPy
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
seq 0 5 95 | xargs -P 4 -n 1 bash ./run_histens_local.sh
```

`-P 4` is an example of four simultaneous batches; choose concurrency to fit your machine. `run_histens_local.sh` activates `cleanPy` itself. It takes a **starting member index**, and each invocation processes five members for all four experiments:

| Experiment | Warming definition | Calibration / reference |
|---|---|---|
| `histens_rpi_real` | Realized | 1850–1900 |
| `histens_rpi_anthro` | Anthropogenic | 1850–1900 |
| `histens_satcal_real` | Realized | 1981–2010 |
| `histens_satcal_anthro` | Anthropogenic | 1981–2010 |

For example, `bash ./run_histens_local.sh 0` processes members 0–4, not the whole ensemble. The full set of batch starts is **0, 5, …, 95**. The comment in the existing shell script says `seq 0 5 100`; use **95** as the endpoint because batch 100 would attempt nonexistent members 100–104. This README does not change that shell script.

The active setting is `regen = 1` in [evaluate_heads_tails_ensemble.py](evaluate_heads_tails_ensemble.py). It reruns the explicitly selected method lists and overwrites their result pickles. `regen = 0` reads cached results. `regen = 2` runs missing methods and can retain old results for methods already present; it does not refresh a correction to an existing method.

Each experiment/member produces:

- `Results/results{experiment}{member}.pickle`: current and retrospective method arrays.
- `Results/headstails_statistics_{experiment}{member}.csv`: member-level statistics and crossing summaries.
- `Results/headstails_method_fineprobs_{experiment}{member}.npy`: member-level exceedance probabilities.

There should be 100 of each type per experiment: **400 pickles, 400 statistics CSVs, and 400 probability arrays** for this workflow. Logs are in `logs_histens/run_0.log`, `run_5.log`, …, `run_95.log`. Wait for all batches to finish before processing or copying results. Check logs for exceptions and check file timestamps: counts alone do not prove that existing files were successfully refreshed. The wrapper can continue to another scenario after a Python failure, so an apparently finished wrapper is not sufficient evidence that all four scenarios succeeded.

### 2. Refresh the paper probability figures and rate CSVs

After the member jobs have finished, remain at the repository root:

```sh
MPLBACKEND=Agg python3 plot_heads_tails_paper_probabilities.py

for experiment in histens_rpi_real histens_rpi_anthro histens_satcal_real histens_satcal_anthro; do
    MPLBACKEND=Agg python3 export_heads_tails_rate_summary.py "$experiment"
done
```

The probability script reads the refreshed member CSVs and NPYs and writes:

- `baseline_1850-1900_heads_tails_probplot.png`
- `baseline_1981-2010_heads_tails_probplot.png`

The rate script reads the refreshed pickles and writes `rates_histens/rates_histens_{experiment}.csv`. Both are needed for a complete paper/dashboard refresh. `MPLBACKEND=Agg` permits saving figures without opening plot windows; omit it when you want interactive Matplotlib windows.

`hist_evaluation_script_rates.py` is a separate historical-observation rate analysis. It is not a required step for assembling the heads/tails dashboard.

### 3. Export the dashboard data from the refreshed Results directory

```sh
python3 climate-dashboard/source/export_dashboard_climate_metrics.py all \
    --source-root "$PWD" \
    --outdir climate-dashboard/data
```

This is the dashboard analysis entry point. It reads the root `Results/`, the newly updated rate CSVs, the saved paper probability arrays, and the calibration files in `Results2/`. It generates four dataset JSONs plus headline, crossing, rate, probability, and member-level CSV exports.

For this normal refresh, **omit `--results-dir`** so the analysis uses the newly regenerated root `Results/`. The previously created `climate-dashboard/data/corrected_results/` copies were a bridge for testing the FaIR retrospective corrections before the root caches were refreshed; they must not silently override a later full regeneration or new year's data.

The earlier discrepancy between the original 0.335652°C/decade PIVW rate and the dashboard's 0.318823°C/decade came entirely from old versus corrected FaIR retrospective caches. With the same 2025 inputs and calibration weights, the refreshed PIVW rate should be approximately **0.318823°C/decade**, with combined 1σ uncertainty **0.054156°C/decade** and a 95% range of **0.212677–0.424970°C/decade**. These are checks for this dataset and correction, not constants to carry into another year.

### 4. Render the dashboard figure library and assemble the website

```sh
python3 climate-dashboard/source/render_dashboard_figures.py
python3 climate-dashboard/source/assemble_dashboard_site.py
```

The first command reads `climate-dashboard/data/` and writes the main figure in PNG/PDF/SVG plus 78 variant PNGs in `climate-dashboard/figures/`. The second regenerates optional rate-history JSON/CSV exports from the validated dashboard data, then copies the exported figures and data into `climate-dashboard/site/dist/assets/` and creates the selected-method CSV downloads. It checks that the complete 78-figure library exists. Both commands derive their default directories from their own script locations, so they can be called from the repository root.

### 5. Review and publish the website exports

The website is static HTML/JavaScript. Browser controls, zoom, hover values, and downloads use exported data; visitors do not run the research Python pipeline.

Temperature plots start in **1965**, while downloads retain earlier rows. The individual ensemble lines are preserved; only a hovered/selected smoothed member can show a temporary pale-yellow approximate 95% calibrated uncertainty ribbon. The interactive benchmark/method rate switches turn sloping segments into horizontal warming-level segments without changing estimates. Static figure downloads show the default rate segments.

Public website CSVs begin with `#` metadata lines describing the early anthropogenic method limitations (FaIR’s special 1930–1964 treatment and GWI’s absence before 1950). Read them using `pandas.read_csv(path, comment="#")`, or skip the leading comment lines in another reader. Analysis CSVs in `climate-dashboard/data/`, `Results/`, and `rates_histens/` retain their normal headers.

The paper/dashboard density curves place annual increments `P(y) - P(y-1)` within `(y-1, y]`. The time-coordinate correction applied on 12 October 2026 moves the curves one year left relative to the earlier plotting convention; density heights, histogram counts, and current exceedance probabilities are unchanged.

Preview the assembled site from the repository root:

```sh
python3 -m http.server 8787 --directory climate-dashboard/site/dist
```

Open `http://localhost:8787`. Check the main PIVW headline, both calibrations, realized/anthropogenic views, annual, method-smoothed, or centred 20-year backgrounds, figure view, and a data download. Stop the preview with Control-C when finished.

On the author's Mac the separate GitHub Pages checkout is `/Users/JohnMatthew/Downloads/climate-state-dashboard`. Copy the complete assembled website into that checkout:

```sh
cp -R climate-dashboard/site/dist/. /Users/JohnMatthew/Downloads/climate-state-dashboard/
cd /Users/JohnMatthew/Downloads/climate-state-dashboard
git status --short
git add app.js member-ribbon.js rate-history.js index.html styles.css assets vendor .nojekyll
git commit -m "Update climate dashboard results"
git push
```

Inspect the status before committing. If working on another machine, substitute your own website checkout path. The `/.` copy includes hidden `.nojekyll`. GitHub Pages should publish the website repository's `main` branch, `/(root)` folder. See [GitHub Pages publishing settings](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site). The same `site/dist/` contents can be served from a lab web server.

## Installation and external data

### Repository and submodules

Clone [jnickla1/Thorne_15](https://github.com/jnickla1/Thorne_15), enter that checkout, and initialize the external method repositories:

```sh
git submodule update --init --recursive
```

The submodules include upstream implementations and exploratory methods. Being present as a submodule does not mean that a method is included in the current dashboard. The authoritative lists for the heads/tails runner are `sel_methods_list_real` and `sel_methods_list_anthro` in `evaluate_heads_tails_ensemble.py`; combination subsets are defined in `plot_fut_results.py` and the dashboard export script.

### Python, R, and optional MATLAB

The author uses Conda. The original development requirements were Python ≥3.11, R ≥4.4.0, and CMake ≥3.2. Some Bayesian change-point experiments also require MATLAB ≥R2016a; saved intermediates are available for some runs. Existing working environments are preferable for reproducing current results.

The checked-in `cleanpy_environment.yml` is an environment snapshot with Linux build pins and an Oscar-specific prefix. It is not a guaranteed portable Mac environment specification. On a compatible system, create an environment with:

```sh
conda env create --name cleanPy --file cleanpy_environment.yml
conda activate cleanPy
```

The YAML calls its environment `cleanpy`, while the Mac wrapper explicitly activates **`cleanPy`**. The cluster wrappers use lowercase **`cleanpy`**. Environment names and platform-specific packages must match the machine being used. The author previously encountered issues with Mamba for this workflow and uses Conda.

R is also called by Python methods. Install the required R packages in the R installation visible to `Rscript`, for example:

```r
install.packages(c("this.path", "mgcv", "abind"))
```

`Bayes_seq_CP_method.py` currently contains a Mac-specific MATLAB executable path (`/Applications/MATLAB_R2019a.app/bin/matlab`); adjust it for your installation when regenerating its intermediates.

The KCC method needs `library(KCC)` to work. Its source is a submodule; if it is not installed in your R library, install it from the repository root with `R CMD INSTALL Methods/43_Forcing_Based/3_Human_Induced/KCC`. Check the submodule's own requirements for additional dependencies.

### External climate data and local paths

The future-model experiments also require the separate [climate_data repository](https://github.com/jnickla1/climate_data). It contains compressed spatial averages and the preprocessing code, rather than all original ESGF files. See [ESGF](https://esgf.github.io/nodes.html) and [NCO](https://nco.sourceforge.net/) for the original-data/preprocessing context.

Create a local `config.py` at the root of this repository. It is ignored by Git so each machine can retain its paths. For example:

```python
import os

CODEBASE_PATH = os.path.dirname(os.path.abspath(__file__))
CLIMATE_DATA_PATH = os.path.expanduser("~/climate_data")
```

Several methods also use paths relative to their own module or hard-coded export names. `config.py` does not replace those method-specific input requirements. A fresh checkout needs the external FaIR, GWI, forecasting, and other method inputs described below; initializing submodules alone does not generate them.

## What is in this directory

| Location | Role |
|---|---|
| `Common_Data/` | Observations, uncertainty files, forcing/covariate data, and preprocessing inputs. Current top-level observation readers use `HadCRUT.5.1.csv` with `Time`, `Anomaly`, `Lower`, and `Upper` columns. |
| `Hist_heads_tails/` | Heads/tails ensemble inputs. `sst_pseudo.csv` is the preindustrial version; `sst_pseudo_1981-2010.csv` uses satellite-era calibration. The first column is year, followed by 100 members. |
| `Methods/42_Temp_Alone/` | Temperature-only estimators, running means, smoothers, state-space methods, and variability adjustments. |
| `Methods/43_Forcing_Based/` | Forcing-based methods, FaIR readers, energy-balance Kalman methods, and GWI/attribution inputs. |
| `Methods/44_EarthModel_CGWL/` | Methods combining observations with climate-model or decadal-forecast information, including UKCP samples and WMO-derived inputs. |
| `Results/` | Cached per-experiment/per-member trends, heads/tails statistics, and probability arrays. A method-code edit does not update these caches by itself. |
| `Results2/` | Historical calibration and fit-probability outputs, including `historical_names_var_scale.csv` and `comb_real_scales.npy` / `comb_anthro_scales.npy`. |
| `Results3/`, `Results3a/` | Outputs from additional simulation and combination experiments. They are not dashboard deployment directories. |
| `rates_histens/` | Four aggregate heads/tails rate CSVs, plus the older dashboard mockup/presentation material. |
| `climate-dashboard/source/` | Separate dashboard analysis, figure rendering, site assembly, and correction-audit scripts. |
| `climate-dashboard/data/`, `figures/`, `site/dist/` | Dashboard analysis exports, figures, and the deployable static website respectively. |
| `climate-dashboard/review/`, `CLIMATE_DASHBOARD_CHANGELOG.md` | Scientific correction checks, retained backups, and change history. |
| `logs_histens/` | Mac heads/tails batch logs. Other cluster wrappers use their own `logs*` directories. |
| `Illustrative_Figures/` | Additional explanatory/paper figures outside the core dashboard workflow. |

Full-method calibration/evaluation runners use `running_subset` to choose directories and can include exploratory implementations. Review that setting and the required external inputs before a full run.

Many runners discover files ending in **`_method.py`** under `Methods/`. Keep backups and alternative copies outside that discovery tree; adding another file with that suffix can add another method and change experiment outputs.

## Add another observational year

An annual update has three distinct sources of work. Acquire and align the inputs first, then rerun calibration/evaluation and the dashboard. The steps below describe the 2025 update and the places to inspect for the next year; the present code has explicit 2025 branches and is not yet a year-independent pipeline.

### 1. Local observation and preprocessing updates

Update `Common_Data/HadCRUT.5.1.csv` and the corresponding HadCRUT annual summary and ensemble uncertainty files. `cent20y_method.py` now computes the running-mean uncertainty using `HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv` and `HadCRUT.5.1.0.0.analysis.ensemble_series.global.annual.csv`. Retain older versions while methods or experiments still reference them.

Refresh both heads/tails ensemble CSVs through the same final year. Keep their year axes and member ordering aligned with the FaIR and GWI exports; array shape alone does not guarantee matching member identity. The `Hist_heads_tails/global_temperature_merge` submodule contains the upstream ensemble construction workflow.

Review other method inputs used by the selected analysis, such as SST, MEI/ENSO, CO₂, forcing, and ocean-heat-content files. The Git history includes changes to these inputs and method implementations; not every changed file is a mechanical annual update.

For the HadCRUT5.1 update, **John manually deleted the cached GAM, quartic/spline, and Bayesian change-point outputs before rerunning the historical analysis**. That cache invalidation was part of the completed 2025 update, not an additional unfinished step for the present FaIR correction.

For a future observation-version or end-year update, review the relevant method-level caches again:

| Reader | Cached output reused when present |
|---|---|
| `Methods/42_Temp_Alone/4_GAM_AR1/GAM_AR1_method.py` (and `GAM_AR0_method.py` if included) | `GAM_AR_Stephenson/output/gamAR*_fits_{experiment}{member}.csv` and the corresponding SE files; the R fit is invoked only if the fit file is absent. |
| `Methods/42_Temp_Alone/2_LT_Fits/quartic_method.py` | `boot_quartic/quartic_se_{experiment}{member}.npy`. |
| `Methods/42_Temp_Alone/4_GAM_AR1/cubic_spline_method.py` (if included) | `AR_cubic_spline/cubic_spline_se_{experiment}{member}.npy`. The quartic polynomial and cubic spline are separate implementations. |
| `Methods/42_Temp_Alone/2_LT_Fits/Bayes_seq_CP_method.py` | `Bayes_Sequential_Change_Point/output/{experiment}{member}BSCtemperatures.csv` and its uncertainty companion; MATLAB is invoked only when the temperature file is absent. |

Archive/invalidate the affected experiment's old intermediates when changing the input observations, and ensure the required R/MATLAB executable is available to rebuild them. A full root-pickle rerun does not automatically invalidate these caches. This concerns the cached estimator outputs, separately from the dashboard's quartic extrapolation of seven non-crossing members.

### 2. FaIR runs on Oscar, performed by John

The constrained FaIR simulations are run separately on Oscar. The local `_unB` method readers consume exported NumPy arrays; rerunning the local heads/tails evaluator does not rerun FaIR on Oscar.

For the present 2025 update, the relevant local directory is `Methods/43_Forcing_Based/1_ERF_FaIR/resliced_headstails25/`:

- `combined_all-headstails25r{run}_{all,nonat,anthro}.npy`: one annual output from each successive observation-cutoff fit; annual labels begin in 1930.
- `retrospective/all-headstails25r{run}cut2025_temp_{all,nonat,anthro}.npy`: full retrospective trajectories from the 2025-constrained fit.
- `retrospective_4/all-headstails25r{run}cut2021_temp_all.npy`: full all-forcing trajectories from the earlier cutoff used by the four-year lag in `FaIR_comb_unB`.

Run **r0** represents the observational run; local heads/tails member **0–99** maps to Oscar run **r1–r100** for these exports. Check the corresponding mapping for every upstream product.

The verified 2025 exports use time boundaries beginning in 1750. Their annual extraction is the following boundary (`year - 1750 + 1`); confusing that row with the preceding boundary caused an earlier discrepancy. The full lagged matrix must come from the **earlier-cutoff fit**, not an earlier row of the latest-cutoff fit. For a final observation year `T`, review the `cutT` and `cut(T−4)` exports and reader branches together, and validate using FaIR's time coordinate before transferring the data. Also review the manual bias/calibration adjustments when changing the observation product.

### 3. Inputs obtained from collaborators or forecast sources

| Input / contact | Local use and update action |
|---|---|
| **Tristram Walsh — GWI** | Updated total/anthropogenic CSVs are read under `Methods/43_Forcing_Based/3_Human_Induced/heads_tails/individual_members25/`. Request the observational run and the matching heads/tails member outputs, retain their baseline/year metadata, and update the year-specific reader branches and calibration offsets. The current readers explicitly recognize 176 annual values for 1850–2025. |
| **Bjorn Samset — Green's functions** | `removeGreensfx_method.py` reads the externally prepared NetCDF output, currently `filterAna_obs_HadCRUT5_1850_2025_v5.nc`. This is a lower-priority input and is outside the default seven-method realized dashboard subset, but is relevant to the broader historical evaluation. Its reader explicitly raises an exception for a later year until new output is added. |
| **Leon Hermanson / WMO decadal forecasts** | The usual input is `Methods/44_EarthModel_CGWL/GlobalT_WMOLC-ADCPforecast_1991-2020.csv`. In the 2025 update, newer forecast information was reconstructed from WMO figures when a direct collaborator update was unavailable. The source figures, reconstruction code, and validation plot are in `figure-based_2026WMOLC-ADCP/`. |

`figure-based_2026WMOLC-ADCP/wmolc_extend.py` reconstructs pseudo-ensemble values from the plotted per-centre mean and min/max, using earlier member distributions as templates. Those are figure-derived reconstructions, not the original new-year member values. Preserve that provenance when adopting them. The script currently reads figures relative to its working directory and writes to `/mnt/user-data/outputs/`; adapt those output paths before a Mac rerun, inspect the validation plot, and explicitly adopt the resulting CSV into the reader's expected input location.

The selected dashboard method **`CGWL10y_sfUKCP` uses UKCP samples**, rather than the WMO CSV. The WMO refresh concerns forecast-based alternatives such as `CGWL10y_for_halfU` and `CGWL10y_forec`, including broader historical comparisons. Distinguish those data requirements when deciding what must block a particular update.

### Year-dependent settings to review before rerunning

Review the input readers first, then the analysis and display endpoints. Relevant examples include:

- `calibrate_method_uncertainty.py` and `evaluate_historical_observations.py`: observation file version, selected folders, calibration choices, plotting endpoints, and dated comparison tables.
- `evaluate_heads_tails_ensemble.py`: selected method lists, `regen`, evaluation windows, and the exclusive upper endpoints in `evalmaxs` (currently 2026 for the 2025 endpoint).
- `export_heads_tails_rate_summary.py`: `last_year` is currently 2025.
- FaIR/GWI/Green's-function readers: 2025 filenames, array-length branches, bias corrections, member mapping, and time boundaries.
- `climate-dashboard/source/export_dashboard_climate_metrics.py`: probability windows, observations, method lists, and the seven non-crossing members' quartic treatment. A later year's set of non-crossing members can differ.
- `climate-dashboard/source/render_dashboard_figures.py` and `climate-dashboard/site/dist/app.js` / `index.html`: explicit 2025 headlines, line endpoints, threshold windows, and other display labels.

Use a targeted search such as `rg -n '2025|2026|175|176|last_year|evalmaxs|cut2025|individual_members25'` on the relevant source files. Review matches individually: historical dates and fixed reference periods must not be replaced indiscriminately.

Useful recorded update checkpoints are [HadCRUT5.1 input/code changes, September 2026](https://github.com/jnickla1/Thorne_15/commit/f0672785aeb014663bec678df34c6d3820fce2b0), [the finalized 2025 update, October 2026](https://github.com/jnickla1/Thorne_15/commit/51de5173496524f11dea9618ea1df63382470580), and [the dashboard addition](https://github.com/jnickla1/Thorne_15/commit/bca152245bd181291e184f5268ff3e3b9799c849). These checkpoints follow the December 2025 work; they are examples of what was changed, not a universal patch to apply for every future year.

## Historical calibration and evaluation

Run this part after updating historical inputs or deliberately changing calibration. It is not required solely to adopt the tested retrospective-only FaIR correction while current outputs and calibration remain unchanged.

1. Run `python3 calibrate_method_uncertainty.py historical` from the root. `historical_regen=True` is required in that script. It runs the chosen methods and writes the `Results2/historical_*` calibration products, including `historical_names_var_scale.csv` used for method variance corrections and PIVW weights.
2. Review `calibrate_historical_combinations.py`: for historical-to-historical calibration use `comparison_type='hh'`. Its current `nmethods_list=[7]` refreshes the realized combination only. To deliberately recalibrate both dashboard combinations, use `[7, 2]` (seven realized methods and two anthropogenic methods), then run `python3 calibrate_historical_combinations.py` and verify both `Results2/comb_real_scales.npy` and `Results2/comb_anthro_scales.npy`. This is a manual setting to review; it has not been changed by updating this README.
3. Run `python3 evaluate_historical_observations.py` to regenerate the historical-observation evaluation, `hist_results.pickle`, statistics, and paper figures. Review `historical_regen` and `running_subset`; the current historical evaluator has regeneration enabled and adjusted variances enabled, so calibration must exist first.
4. Run the heads/tails and dashboard refresh above using the new inputs and calibration.

`hist_evaluation_script_rates.py` can additionally be run for the separate historical-observation rate figures. Its own regeneration flag controls whether it reuses `hist_results.pickle`. It is not a substitute for either `export_heads_tails_rate_summary.py` or the dashboard exporter.

## Future simulations and combination experiments

These are separate paper-validation workflows. They are not required for updating only the current dashboard exports. They need the external model data, their own experiment settings, and usually Oscar/Slurm.

The existing future-method batch wrapper is `run_parBATCH_future.sh`. It evaluates MPI-ESM1-2-LR SSP126/245/370 and NorESM RCP45 volcanic variants using `evaluate_future_scenarios.py`. Review the cluster account, module loads, environment name, resources, and method-regeneration settings before submitting it:

```sh
sbatch run_parBATCH_future.sh
python3 plot_fut_results.py
```

Run plotting after the submitted jobs have finished. The combination-focused future figures use:

```sh
python3 evaluate_future_scenario_subsets.py futplotcomb_ESM1-2-LR_SSP126_constVolc -11
python3 evaluate_future_scenario_subsets.py futplotcomb_NorESM_RCP45_Volc -11
```

The negative-member convention is script-specific; inspect each runner before assuming the same batch size everywhere. In the heads/tails runner, a negative value `-(member+1)` selects that single member and enables plotting. For `futplotcomb` experiments, the existing NorESM plotting branch overrides some supplied scenario/member choices; the MPI-ESM branch uses the supplied member for highlighting.

For the simulation-based leave-one-out combination experiments, the existing entry points include:

```sh
sbatch run_parBATCH_hist_fitprob_ff.sh
sbatch run_LOO2.sh
sbatch run_LOO16.sh
sbatch run_LOO8.sh
```

After their outputs are complete, the aggregation script is spelled **`create_combination_method_LLOffcombine.py`**:

```sh
python3 create_combination_method_LLOffcombine.py -2
python3 create_combination_method_LLOffcombine.py -16
python3 create_combination_method_LLOffcombine.py -8
```

The wrappers retain their existing research-specific meanings: `run_LOO2.sh` currently calls `fut_2_script.py`, whereas the 8/16 wrappers call `create_combination_method_LLOffcompute.py`. Inspect these settings before extending an experiment. The commands above preserve the existing entry points rather than implying a new uniform API.

Additional figure scripts include `Illustrative_Figures/betts_illustrating_figv2.py` and `Illustrative_Figures/sect42_temp_record.py`. Their paths and assumptions should be reviewed separately when the paper's data year changes.

## Interpretation, caches, and troubleshooting

### Dashboard scientific conventions

The definitive historical crossing benchmark uses the realized observational **centered 20-year mean**. Its rate uses **centered 10-year means**, fitted around the crossing. The dashboard interpolates upward crossings and, for the current preindustrial 1.0°C benchmark, includes the seven requested quartic estimates from the last 25 available centered 10-year means. Those estimates and their limitations are identified in the downloads. Anthropogenic dashboard views retain this realized observational benchmark as an explicitly labeled comparison.

Current dashboard warming levels use **current method outputs**. The current-year rate is a quadratic endpoint derivative from the last **15 retrospective annual outputs**; the GWI readers temporarily use their current outputs as a retrospective proxy. The rate uncertainty combines within-fit variance, weighted between-method disagreement, and between-member spread; the displayed 95% range is rate ±1.96 times this combined 1σ uncertainty. `rate_sd` alone is only the member spread.

In a `export_heads_tails_rate_summary.py` row whose `level` is 2025, the misleadingly named `year_cross` column now contains the **current warming level**, not a crossing year. The 12 October 2026 reporting correction separates this level from the retrospective series used to fit the rate. `year_cross_sd` uses the current member endpoints as well. Older saved CSVs used retrospective levels; the dashboard exporter accepts both formats and marks the source in `saved_summary.level_source`. The historical crossing rows still retain their existing calculation conventions; use the dashboard's final data exports for its definitive crossing benchmarks.

The rate script also prints a separate per-method/current-combination endpoint comparison: current and retrospective means, paired differences, maximum differences, matching counts at absolute tolerance 1e-6°C, and missing endpoints. These are diagnostic checks, including expected CGWL window/FaIR sample-noise differences, not assertions that stop a run. They are **not added to the CSV**. Individual methods compare their raw retrospective endpoints; combination comparisons include the existing current-proxy fallback used by the rate fits.

Satellite-calibrated views use the source pipeline's 0.7°C reference shift. The Temperature series menu changes only the grey/black ensemble lines, not the calculations. Distribution panels retain missing mass from non-crossing members and share a 0–0.3 probability/year scale, with explicit right-panel magnification. Details are in [climate-dashboard/README.md](climate-dashboard/README.md).

### Which dashboard files to use

Use `climate-dashboard/data/*_corrected_headlines.csv` and each dataset JSON's `headline_summary` for the final dashboard levels and rates. `*_legacy_saved_rates.csv`, `rates_fixed/`, and the original/fixed script copies retain intermediate or historical audit work; they are not additional production entry points. Use `*_crossings_rates.csv` and `*_member_crossings_rates.csv` for the definitive and method-comparison crossing calculations.

`climate-dashboard/source/audit_and_rebuild_retrospective_caches.py` is a one-time correction/regression helper, not a routine annual-update command. It compares the edited FaIR readers with retained backups, guards current outputs, and creates separate corrected pickle copies. Its 2025 assumptions and baselines need review before reuse with a later year. Once the root `Results/` has been fully regenerated, the regular dashboard exporter should read that root directly. `export_fair_same_cut_lag4.py` is a superseded diagnostic: the adopted lagged input is the full earlier-cutoff FaIR matrix from Oscar.

### Common causes of mismatches

- **Method files were edited but old pickles are still being read.** Regenerate the affected root caches or explicitly select the verified corrected cache directory during a correction audit.
- **A method has its own stale intermediate files.** Check its reader/generator, including GAM, quartic/spline, and Bayesian change-point outputs, plus externally supplied FaIR/GWI data. The author cleared the first three cache families for the HadCRUT5.1 update.
- **Only some members or scenarios finished.** Check logs and timestamps before regenerating aggregate CSVs or publishing.
- **CSV and pickle vintages differ.** The dashboard exporter checks the saved levels against current outputs or retained legacy retrospective endpoints; regenerate the four rate CSVs after refreshing pickles.
- **The site still has old exports.** Run `assemble_dashboard_site.py`, copy the complete `site/dist/` into the Pages checkout, then commit and push that checkout. For JavaScript changes, update the script version query in `index.html` to refresh browser caches.
- **New-year arrays or forecasts are too short.** Review external exports, explicit reader branches, probability-window endpoints, and method-specific intermediates before rerunning the full evaluation.

### Existing result and note links

The working [historical comparison](https://docs.google.com/spreadsheets/d/10izz9VruI9L1pNT3pwKLlNPVhzrvGRdYk3VxvdQ1es8/edit), [MPI-ESM1-2-LR comparison](https://docs.google.com/spreadsheets/d/1eWAeL1HHHSqyL1YF2IYaQwgMQuh8Y1RJ/edit), [NorESM comparison](https://docs.google.com/spreadsheets/d/1gHNtpZ4MVIw_NYp62kjtuCHZ2kZcWYwo/edit), and [progress/citation notes](https://docs.google.com/spreadsheets/d/1iShljXO2rmPHpPjPkBGbPSwZc_7XVtCU32O3sn6sjTc/edit) are retained for context. They are author-maintained working documents; confirm their data vintage when comparing them with regenerated local outputs.

### Optional rates-over-time panel

The interactive warming-rate card opens a panel on demand. The blue series is the benchmark OLS slope through 11 centred 10-year observational means (year − 5 through year + 5), and requires observations through year + 10, so currently ends in 2015. It represents realised observations on both warming views. The red series uses trailing 15-year quadratic endpoint fits to the existing retrospective histories, with the same GWI current-output proxy as the headline; its 2025 endpoint matches the headline. These are revised historical estimates, not vintage predictions. A checkbox adds a faint dotted red comparison fitted to annual current outputs, while the retrospective series stays visible.

`climate-dashboard/source/export_dashboard_rate_history.py` creates four small `*_rate_history.json` files and 26 selected-method `*__rate_history.csv` downloads, using the already exported main JSON rather than rerunning methods. CSVs include the blue benchmark and both red conventions, with units and explanatory comments. It is called automatically by `assemble_dashboard_site.py`, so the normal update sequence is unchanged. To regenerate only these analysis exports, run `python3 climate-dashboard/source/export_dashboard_rate_history.py`; run the assembler afterward to update the website. The browser imports `rate-history.js` and fetches rate data only when the panel opens. There is no static-figure counterpart.

### 13 October 2026 Kalman Jacobian correction

The nine active `Methods/43_Forcing_Based/2_Kalman/EBMKF_Nicklas*.py` modules now use `1+Cs/Cd` in their temperature/ocean Jacobians and omit the spurious standalone `+Cs` from J_HT. This changes covariance propagation and hence assimilation: current as well as retrospective means/uncertainties may change. Physical update formulas are unchanged. `evaluate_historical_observations.py` is temporarily configured for only `1_Run_Means` and `2_Kalman` as the first diagnostic. Its subset run overwrites `hist_results.pickle`; restore the full method list for final paper outputs. Existing variance scales remain provisional until recalibration.

After checking those results: recalibrate historical method uncertainty and combinations, regenerate the full historical evaluation, regenerate heads/tails caches/probability arrays, refresh paper probability figures and rate CSVs, then export/render/assemble the dashboard using the root Results directory. Rebuild Kalman-containing future/leave-one-out evaluations used in the paper. Also regenerate the attribution `Nicklas_GMST_timeseries.csv` and `Methods/43_Forcing_Based/3_Human_Induced/anthropogenic-warming-assessment/multi_method_assessment_newStyle3.py` figure. Its separate attribution producer (`EBMKF_ta4_attibutable.py`) requires review of its legacy 2023 cuts, explicit Kalman initialization and conversion from its flat CSV to the figure reader's multi-index schema before rerunning. Historical pickles alone do not update that attribution input. No website publication or scientific rerun was performed as part of this patch.


### 14 October 2026: future leave-one-out recalculation

All 270 primary Oscar method-calibration sets were downloaded and checked for complete NumPy payloads and consistent ordering of 60 methods. The approved future-only corrections in `create_combination_method_LLOffcompute.py` align all training arrays to a common calendar, exclude the held-out member from trial samples and inverse-variance uncertainty calibration, and evaluate the intended full/last-75 KL windows. The script now reads `current_methods_statistics_261013True.csv`. `create_combination_method_LLOffcombine.py` verifies exact expected member coverage rather than requiring December 2025 file dates. Historical and heads/tails method indexing is unchanged.

The local batch reruns all 270 members for legacy selection modes 7, 8, 16 and 32, plus the existing mode-2 comparison. Selection-mode numbers are identifiers rather than guarantees of the number of available methods. The aggregate `final_combined_*methods_wide_ff.csv` files should only be used after every expected member succeeds and validates; partial runs must not be reported as final. Detailed progress, source backups, the approved patch and before/after comparison artifacts are in the author review workspace's `outputs/future-loo-2026-10-14/`. These changes combine corrected Kalman calibration and corrected LOO validation, so changes in future comparison metrics must not be attributed to the Jacobian alone.


## Kalman v11 and frozen 9c inputs — 17 October 2026

The active model-based Kalman cores now use the v11 equation, capacity, forcing-reference, Jacobian and smoother corrections. Their state arrays follow the frozen input length. This update can change both current and retrospective estimates; the previous calibration products and published figures are not refreshed by changing source files alone.

`Common_Data/toyKFmodelData9_ext.xlsx` is the calculation workbook; `Common_Data/toyKFmodelData9c.csv` is its frozen 176-row, nine-column export for 1850–2025. Python readers formerly using 8c now use 9c. The old workbook is retained in `Common_Data/archive/pre-v11-2026-10-17/`; the old CSV is retained for comparison. Do not substitute the shorter 9a CSV for production 9c. The shared `ebmkf_v11_support.py` reads the frozen CSV without silently replacing its temperature columns from another file.

Column H through 2024 comes from the authoritative AE (ARI + ACI) column in `~/Dropbox/KalmanFilterPython/ERF_best_aggregates_1750-2024.xlsx`. Its source values and derived formulas are internalized in the extension workbook. H2025 is an estimate obtained by extending the **2015–2024 ten-year linear trend**, with an explicit `SLOPE`/`INTERCEPT` formula in H176 and explanatory notes below the data. B/E use HadCRUT5.1 observations/standard errors. F/G preserve the original IAP splice and extend it through 2025; the ta4 OHC ensemble remains separate. D2024–2025 uses GloSSAC v2.24 with the documented 525-nm conversion once. I preserves the source through 2022 and uses complete annual CDR v03 solar values for 2023–2025. C2025 retains the established equivalent-CO2/observed-CO2 trend extension, refreshed through 2024. Columns E/G contain standard errors, which the methods square when constructing variances.

For the first visual review, `evaluate_historical_observations.py` is configured for only `1_Run_Means` and `43_Forcing_Based/2_Kalman`, with regeneration enabled. From the repository root:

```bash
conda activate cleanPy
python3 evaluate_historical_observations.py
```

This subset run overwrites `hist_results.pickle` with the subset. Existing uncertainty calibration remains provisional; restore the full method list and regenerate calibration before final paper numbers. No full scientific rerun or publication was performed for this migration.

Validation of the frozen input found no formula errors and preserved both original workbook charts. All nine active cores pass independent Jacobian and aerosol-unit-response checks, with covariance eigenvalues nonnegative within numerical tolerance. All three historical wrappers run through 2025 with finite outputs, repeatable calls, independent returned arrays, matching final current/retrospective endpoints and correct behavior under a baseline shift. These tests do not establish full future/scenario correctness.

The ta2 feedback scan against the HadCRUT5.1 centered-20-year benchmark over 1860–2015 gives historical-data fit candidates fdbkA≈0.638 for 2a/2c and ≈0.972 for 2b. These are fit candidates for author review, not installed defaults or future-model calibrations. The author-requested starting value 0.2 remains in ta2 and ta4. The distinct MPI heat-capacity/deep-ocean/exchange settings are retained.

After visual review and final parameter selection, rebuild historical uncertainty calibration and combinations, full historical evaluation, heads/tails method caches and combinations, crossing probabilities/rates, attribution input and assessment figure, future-method calibration and evaluations, future combinations/leave-one-out results, and dashboard data/figures. Also review the temperature-only Kalman and MEI consumers whose input source changed to 9c. All outputs depending on those inputs are potentially stale. Keep the paused future controllers stopped until deliberately restarting with the final configuration.

The attribution producer `EBMKF_ta4_attibutable.py` still requires its own input-window/component-ledger review before regenerating `Nicklas_GMST_timeseries.csv`; in particular its 2023 cuts and the `forcing_O3` assignment to the solar column must not be mistaken for validated attribution physics. Attribution and future/scenario smoke tests are outstanding and do not block the deliberately narrow historical visual review.


## Ta4 standalone diagnostics and smoother review — 17 October 2026

From the repository root, run `python3 Methods/43_Forcing_Based/2_Kalman/EBMKF_ta4_method.py` to display three diagnostic figures with `plt.show()`: original/adjusted inputs and parameters, filter/retrospective/blind states, and the blind shortwave/longwave/storage budget. Optional `--save-dir <directory>` saves PNGs, a diagnostic CSV and source/parameter JSON; `--no-show` supports a headless check. This is an opt-in historical run, not a publication/cache rebuild. The blind replay uses the actual wrapper-adjusted forcings. The SW multiplier is labelled as a multiplier rather than a cloud-feedback coefficient.

`precompute` is the local function in `Methods/43_Forcing_Based/2_Kalman/ebmkf_v11_support.py`; `_v11_precompute` is its import alias. `EBMKF_Nicklas4.py` now exposes inclusive `equilibrium_reference_years=(1850,1899)` and `radiation_reference_years=(1999,2005)`. Frozen 9c column B matches local HadCRUT5.1 exactly. The legacy absolute alignment makes the 1960–1989 mean 13.85°C, so the current reference is Teq1850=286.6475103 K (13.4975103°C). This is an initialization/calibration reference, not a solved equilibrium for actual 1850 forcing. The older `HadCRUT5_capacity_reference_5.0.2.0.csv` only determines non-ocean heat capacity and does not set Teq1850 or T02. No temperature datum was changed in this diagnostic work.

The ta4 failure at fdbkA=0.42 was reproduced: the early prediction covariance condition number is about 1.2e14, and a general explicit inverse corrupts the RTS subtraction at initialization. The ta4 RTS gain now uses a symmetric Hermitian pseudoinverse instead of an explicit general inverse/borrowed-next-gain fallback. Covariance checks remain active; no negative variances are hidden or clipped. In the isolated 0.42 comparison the forward states/covariances were identical and the maximum retrospective temperature change was about 4e-8°C. This makes 0.42 numerically usable without endorsing it as the best-fitting parameter. Ta4's existing 0.2 setting was retained.

The historical ta2 0.2 default was inappropriate after the physics migration. The 1860–2015 centered20 comparison gives RMSE 0.07747°C at 0.2, 0.06255°C at 0.42, 0.06099°C at its former 0.45, and 0.05599°C at the newly installed historical 2a fit fdbkA=0.63782768485. The 1885–1909 legacy initialization was also tested and worsened the fit; it was not silently restored. These are calibration results, not independent validation or recovered old-paper statistics. MPI/NorESM 2b/2c feedback calibration remains outstanding and was not replaced by the historical 2a fit.


### Completion of the smoother numerical repair

Testing both feedback settings showed that a Hermitian inverse alone was insufficient at 0.2: cancellation in the standard RTS covariance subtraction remained. Ta4 now uses the Joseph conditional-covariance form `(I-JF)P(I-JF)' + J Qprocess J' + J Phat_next J'`, algebraically equivalent to the standard RTS expression with the optimal gain. This avoids subtracting nearly equal covariances; it does not clip negative eigenvalues or weaken their checks. The inverse uses the symmetric eigensolver. Historical tests at both 0.2 and 0.42 now pass, with smallest smoothed eigenvalues at roundoff (about 1e-23). Both leave forward states/covariances identical to the corresponding pre-repair case. The 0.42 retrospective temperature change remains approximately 4e-8°C. This is a numerical repair rather than a reason to prefer one physical feedback value.

The historical ta2 fit has a material endpoint tradeoff: changing 2a feedback from 0.2 to 0.63782768485 raises the 2025 temperature anomaly from 1.09071937 to 1.30320850°C (about +0.21249°C), while reducing 1860–2015 benchmark RMSE from 0.07747 to 0.05599°C. The definitive centered20 benchmark ends in 2015; the later endpoint is not validated by that fit. The installed historical fit is explicitly a candidate for author review, not a final model selection. Ta4 remains at 0.2 and MPI/NorESM fits remain separate.


## Symmetric covariance solver rollout — 17 October 2026

All nine active model-based `EBMKF_Nicklas*.py` cores now use the local shared `covariance_inverse` helper for innovation and RTS prediction covariance inversions. It checks finite/positive-semidefinite covariance and uses the Hermitian eigensolver pseudoinverse. Generic explicit inverses and borrowed-next-year RTS gain fallbacks have been removed from these cores. Their smoother covariance uses the stable Joseph conditional form and the exact process covariance saved from that variant's actual forward prediction for each year. State dimension, process-noise settings and physical parameters are preserved. Archived implementations and temperature-only Kalman methods were outside this bounded rollout.

Direct historical before/after checks pass for all nine cores; all three historical wrappers also pass repeated-call, independent-array, baseline-shift and terminal-current/retrospective checks. Maximum temperature differences across the nine direct comparisons are 1.75e-09°C current and 8.32e-09°C retrospective. The ta4 historical check at fdbkA=0.42 also passes. Covariance checks remain active; this is not a waiver for physically invalid covariances. No full scientific calibration, future batch, dashboard rebuild or publication was launched.


## Historical ta2 aerosol-forcing percentile sensitivity — 17 October 2026

Set the single `AEROSOL_FORCING_PERCENTILE` parameter near the top of `Methods/43_Forcing_Based/2_Kalman/EBMKF_Nicklas2a.py` (default `50.0`, permitted range 5–95). This knob controls observational ARI+ACI forcing for the historical ta2/2a pathway, independently of its cloud-feedback parameter. The model-specific future pathways retain their own prescribed forcing. No automatic percentile optimization or input-workbook/CSV overwrite is performed.

The three `ERF_{p05,best,p95}_aggregates_1750-2024.csv` files reside beside the method. The matching best CSV was copied from the local KalmanFilterPython folder. Their `aerosol` column already contains combined ARI+ACI; use its supplied aggregate percentiles rather than adding the two marginal component percentiles. Linearly interpolate p05→best→p95, treating best as an approximate percentile-50 anchor (not claiming it is the exact median). Higher settings generally give less-negative aerosol ERF. Use the selected curve consistently for the 1999–2005 centering and radiation-coefficient recalculation. Annual labels ending in .5 denote the same calendar year. Each curve is extended to 2025 using its own 2015–2024 linear trend, matching the author-approved central extension convention.

The downloaded aggregate 5th–95th interval is broad: 1980 best −1.3939 W/m² with range −2.2691 to −0.6024, and 2024 best −1.0662 with range −1.8956 to −0.4263. These are pointwise forcing-percentile ranges, not a simultaneous confidence band for a whole forcing trajectory. Source: https://github.com/ClimateIndicator/forcing-timeseries/tree/main/output .

Checks at 5, 50, 75 and 95 pass, reference centering is consistent, invalid settings are rejected, and returning to 50 restores the original default outputs exactly. The repaired workbook, frozen 9c, ta4 and future/model-specific feedback settings remain unchanged by this addition. Rerun the historical evaluator after choosing a different setting; existing saved results do not update automatically.


The aerosol interpolation/application helpers now live in `Methods/43_Forcing_Based/2_Kalman/ebmkf_v11_support.py`. The single percentile setting lives in `EBMKF_Nicklas2a.py`; its `precompute_coeffs()` applies the selected aerosol curve before recalibrating radiation coefficients. This applies to both standalone runs and the historical ta2 wrapper. For standalone inspection, edit that setting, then run `python3 Methods/43_Forcing_Based/2_Kalman/EBMKF_Nicklas2a.py` from the repository root in cleanPy; the existing plots use `plt.show()`. Existing standalone plots are for the core's unsmoothed annual inputs; the ta2 wrapper additionally applies its established trailing volcanic smoothing. Standalone plot export filenames are retained and are written to the current working directory.


## Historical ta2 eCO2 forcing sensitivity and relative parameter likelihood

`EBMKF_Nicklas2a.py` now places `eCO2_FORCING_PERCENTILE=50.0` beside `AEROSOL_FORCING_PERCENTILE` and `fdbkA`. Shared implementations are in `ebmkf_v11_support.py`; both standalone and historical-wrapper runs apply these settings before radiation calibration. Existing author choices of aerosol percentile 75 and fdbkA=0 were preserved. Higher eCO2 percentiles select stronger non-aerosol anthropogenic forcing; the permitted range is 5–95.

The eCO2 channel follows workbook column Z: CO2 + CH4 + N2O + O3 + contrails + land_use + BC_on_snow + H2O_stratospheric + halogen. It excludes aerosol, solar, volcanic and the duplicated minor subtotal. The workbook conversion is log10(eCO2)=(ERF+31.14922366816877)/12.743331707472763 (V100/V99). The best-component sum was checked against frozen 9c eCO2 through 2024. Gaussian component sigma is (p95−p05)/(2 Phi^-1(0.95)); sum uncertainty is combined in quadrature, explicitly assuming independent components. Select ERF=best+Phi^-1(p/100)*sigma and perturb frozen log10(eCO2) by deltaERF/V99. This is an approximate Gaussian uncertainty construction, not a source-provided joint percentile. The central frozen 2025 workbook extrapolation is preserved exactly; only the uncertainty width is continued using its 2015–2024 trend. Setting 50 reproduces previous states exactly. The workbook and frozen CSV were not edited.

`parameter_choice_density` returns detailed diagnostics but prints a dimensionless relative parameter likelihood, maximum 1 at the assumed Gaussian peak. It evaluates the actual selected aerosol ERF amplitude, the non-aerosol ERF amplitude and fdbkA under independent Gaussian marginals at reference year 2024. Aerosol sigma uses the supplied total-aerosol 5–95 range; the non-aerosol sigma uses component quadrature. Cloud feedback uses AR6 WGI Table 7.10: mean 0.42, 5–95 range −0.10 to 0.94 W m−2 K−1, equivalent Gaussian sigma 0.3161375526. Source: https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-7/ . The printed ratio is exp(−0.5 sum(z_i²)); it is a parameter-density ratio in ERF/feedback coordinates, not the probability of model correctness, a fit-to-observation likelihood, a percentile-space density, or a joint probability for every year. Correlations among forcings or with feedback are not estimated. Aerosol's existing piecewise-percentile mapping remains approximate and is evaluated at its actual forcing amplitude.

At aerosol percentile 75 / eCO2 percentile 50 / fdbkA=0, the ratio is 0.301405 (30.14% of peak). At the three central distribution values it is exactly 1; shifting only cloud feedback by one sigma gives exp(−0.5). Tests at eCO2 percentiles 5/50/95 and restoration to 50 pass with finite results. These diagnostics select no optimized parameter setting.

The missing blue innovation mixture in 2a panel b came from commented-out density calculations. These now use the active temperature/OHC observation variances divided by forecast variance; all n−1 prediction innovations contribute equally. The accumulators reset for each run. Both temperature and OHC curves are finite, integrate to approximately one on the plotting grid and reproduce on repeat calls. Standalone main/plot generation was executed and visually inspected. Future-specific 2b/2c and ta4 physics are unchanged by this addition.


## Version-2 reference-preserving volcanic attenuation — 17 October 2026

The author explicitly requested the v10 reference-preserving sensitivity in all three version-2 cores, superseding the earlier original-volcanic-function choice for these variants only. `EBMKF_Nicklas2a.py`, `2b.py` and `2c.py` now define `volcanic_target_efficiency=20.0` W m−2 per raw AOD and retain an uncertainty descriptor of 5.0. This is a deterministic central setting, not automatic uncertainty sampling or an additional term in the existing three-parameter likelihood score. The magnitude can be set to 15 or 25 for sensitivity checks.

The exact mapping found in Corrigendum_2026-10-13/code/ekf_testv10.py is D(tau)=9.7279+tau_ref+s(tau−tau_ref), where tau_ref is the raw 1999–2005 mean and s=target*(9.7279+tau_ref)/Sref (including any mean reference aerosol anomaly as in v10). With the current normalized reference, s=0.8099629413, approximately 19% attenuation. At tau=tau_ref the transmission remains 9.068/(9.7279+tau_ref)=0.9317642805 exactly; this is not a new hard-coded 0.93 reference or a direct multiplication of raw AOD by 0.8. The local shortwave derivative at T02 and centred reference aerosol forcing is −20 W m−2 per raw AOD. Slopes away from the reference vary with state and forcing.

Historical updates, their Jacobians, core future routines, constant/generated volcanic inputs, and ta2 trailing averages use this same mapping. Reciprocal averaging is inverted back to raw AOD using its exact reference-preserving inverse, avoiding double attenuation. The historical reciprocal mean is refreshed when coefficients are recomputed. The existing 525-to-550-nm input conversion is not repeated. Base ta, ta4 and other model variants retain their previous volcanic choice. Existing aerosol/eCO2 percentile knobs and cloud-feedback choices are preserved.

Numerical checks in all three version-2 cores verify reference transmission, raw-AOD inverse round trips, local slopes at target magnitudes 15/20/25, trailing-average consistency and synthetic constant/generated future paths. All three historical wrappers and all nine core Jacobian/unit-response checks pass. These are targeted checks, not a full MPI/NorESM scientific rerun. Frozen inputs/workbook are unchanged; dependent version-2 results need regeneration when the author selects final parameters.


## Ta4 forcing controls and single-case future tuning — 17 October 2026

The recent 2a forcing controls, relative parameter-likelihood/1D equivalent-percentile report, innovation-mixture reset, and reference-preserving volcanic attenuation now also apply to `EBMKF_Nicklas4.py` and its ta4 wrapper. The augmented four-state structure, OHC ensemble inputs, existing wrapper smoothing/scaling choices and explicit temperature reference periods are retained. The mapping is used consistently in the Jacobian, physical updates, raw-AOD reciprocal smoothing/inversion, volcanic observation corrections and diagnostic budgets. The default volcanic target is 20 ± 5 W m−2 per raw AOD; ±5 remains an uncertainty descriptor rather than stochastic sampling. Ta4’s established extra smoothing and ×1.4 raw-AOD adjustment remain explicit wrapper choices; standalone core and wrapper plots therefore have different adjusted inputs.

Starting settings copied from the author's current 2a into **2b, 2c and 4** are `AEROSOL_FORCING_PERCENTILE=85.0`, `eCO2_FORCING_PERCENTILE=17.0`, and `fdbkA=0.00`. These are starting values for review, not an optimization or a recommendation to regenerate all outputs yet. Edit them independently in each core. The MPI 2b variant retains its fitted total heat capacity 8.5, depth 1841 m and gamma 1.1; 2c retains total capacity 10.8661063, depth 1141 m and gamma 0.67. 2a was not edited in this step.

**Future forcing convention (corrected later on 17 October):** `AEROSOL_FORCING_PERCENTILE` selects both the updated observational aerosol prefix and the matching SSP percentile trajectory through the end of each future run. The SSP selection uses `ERF_<scenario>_1750-2500.csv`, `_pc05.csv` and `_pc95.csv`; sum ARI+ACI within each file, then use the same piecewise-linear 5/50/95 interpolation convention as the historical control. These component-sum curves are an approximation to a combined aerosol percentile, not an inferred joint distribution. The historical source supplies total-aerosol percentile curves directly, so source-release and uncertainty-construction differences can still produce a join mismatch. The eCO2 control remains observational; this correction does not change its model-specific future construction. Future ta4 coefficient recalibration preserves its full-length forcing arrays. Review the input panels, including the handoff, before selecting final parameters.

From the repository root in `cleanPy`, inspect exactly one member, without writing `Results/`, `Results2/` or evaluation tables:

```bash
conda activate cleanPy
# MPI selects EBMKF_Nicklas2b:
python3 Methods/43_Forcing_Based/2_Kalman/EBMKF_ta2_method.py fut_ESM1-2-LR_SSP245_constVolc 0
# NorESM selects EBMKF_Nicklas2c:
python3 Methods/43_Forcing_Based/2_Kalman/EBMKF_ta2_method.py fut_NorESM_RCP45_Volc 0
```

The final integer is the zero-based member number. MPI SSP126/370 and NorESM VolcConst are also supported. These diagnostics reuse `evaluate_future_scenarios.collect_data/select_data` and its uncertainty construction, including its existing observational splice convention; they do not silently update that evaluator's legacy HadCRUT input. The component panels use the active arrays through 2100 (MPI) or 2099 (NorESM), showing temperature/centered-20-year benchmark, OHC, AOD, solar, eCO2, aerosol ERF, radiative fluxes and storage terms, plus innovation histograms/observation mixtures. Parameter-choice likelihood evaluates historical forcing amplitudes and feedbacks, not the likelihood of the simulated future member. The plot explicitly identifies that the aerosol percentile also applies to SSP inputs.

Append `--breakpoint` to show the plots non-blocking and pause inside `run_method`, where `ekf`, `years` and `new_observ` are available. At `(Pdb)`, use `p ekf.fdbkA` (or inspect forcing/state arrays), `c` to finish, or `q` to abort. Edit parameters in 2b/2c and rerun the command for each trial. Append `--save-dir /absolute/output/folder --no-show` for PNGs without GUI windows. Neither option is active in normal batch runs. The existing evaluator can also enable the hook with `EBMKF_DEBUG_PLOTS=1` and optionally `EBMKF_DEBUG_BREAKPOINT=1`.

The **production** one-member evaluator remains:

```bash
python3 evaluate_future_scenarios.py fut_ESM1-2-LR_SSP245_constVolc 0
```

It runs the configured `running_subset` (currently running means and model-based Kalman methods), performs evaluation and **writes** the member result cache and tables. `THORNE_REFRESH_METHODS=EBMKF_ta2` refreshes just ta2 when that member's existing cache is present; use the diagnostic command above for tuning without result writes. Future uncertainty-calibration jobs (`calibrate_method_uncertainty.py`) are a separate later step.

Validation: all nine core Jacobian/unit-aerosol-response checks and all three historical-wrapper repeatability/baseline/endpoint checks pass. The new diagnostics ran real MPI SSP245 and NorESM volcanic member 0; the ta4 MPI member-0 path preserved all 251 forcing values and returned finite results. Volcanic reference transmission, inverse mapping, and −20 reference derivative pass in 2b/2c/4. Historical ta4 diagnostic budget closes to <2e−13 W m−2. These are targeted checks, not a fresh full evaluation/calibration/combination/dashboard rebuild. Saved publication and dashboard outputs still need regeneration after final parameters are selected. Backups and validation products: `/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-v11-migration-2026-10-17/port-to-four/`.


### SSP aerosol percentile correction — 17 October 2026

The future aerosol helper previously always read the central SSP curve even when the historical prefix used another percentile. It now honours the core's existing `AEROSOL_FORCING_PERCENTILE`: linear p05→best between 5 and 50, and best→p95 between 50 and 95, followed by subtraction of the same historical aerosol reference. MPI uses the requested SSP scenario; both NorESM RCP4.5 volcanic variants use SSP245 as before. The shared correction applies to ta2, ta4 and any caller with that control. Cores without a percentile setting retain the previous central estimate exactly. Missing/unordered/nonfinite curves and invalid percentiles fail explicitly; no silent central-value fallback is used for a non-central request.

Author settings were preserved: 2a = (aerosol 85, eCO2 17, fdbkA 0); MPI 2b = (85, 50, 0.40); NorESM 2c = (85, 17, 0); ta4 = (85, 17, 0). No percentile values, frozen inputs or SSP files were edited. No anchoring offset or new smoothing was applied. For SSP245 at aerosol percentile 85, the centred forcing is +0.099949 W/m² in historical 2023, +0.340612 in scenario 2024 (previously −0.344306), and +0.497524 in 2100. A +0.240663 join mismatch remains; consistent percentile selection does not itself harmonize the two releases or their combined-uncertainty construction.

Validation: 15 full-period curve checks (SSP126/245/370 at percentiles 5/17/50/85/95), exact central compatibility without a knob, both NorESM scenario mappings and invalid-input rejection pass. Real MPI SSP245 and NorESM volcanic member-0 diagnostic runs return finite outputs and refreshed component PNGs. Normal evaluation/calibration caches, manuscript figures and the dashboard were not regenerated. Future calibration/evaluation products need rerunning after final settings and any desired join treatment are agreed. Backups, preserved settings and validation are in `/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-v11-migration-2026-10-17/scenario-aerosol-percentile/`.


## Active ta4 diagnostics and isolated inferred correction — 17 October 2026

`EBMKF_ta4_method.py` now provides the same single-case diagnostic entry point as ta2, for historical observations and MPI/NorESM future members. These run the actual ta4 wrapper preprocessing, filter and retrospective smoother; they do not merely plot the unfiltered core. From the repository root in cleanPy:

```bash
# Historical observations (also the default with no positional arguments):
python3 Methods/43_Forcing_Based/2_Kalman/EBMKF_ta4_method.py historical
# One MPI future member:
python3 Methods/43_Forcing_Based/2_Kalman/EBMKF_ta4_method.py fut_ESM1-2-LR_SSP245_constVolc 0
# One NorESM member:
python3 Methods/43_Forcing_Based/2_Kalman/EBMKF_ta4_method.py fut_NorESM_RCP45_Volc 0
```

MPI SSP126/370 and NorESM VolcConst are also available. Append `--breakpoint` to pause inside the actual `run_method` after creating plots; inspect `ekf`, `years`, `new_observ` and `diagnostic_series`, then `c` to continue or `q` to abort. Append `--save-dir /absolute/output/folder --no-show` to save without GUI windows. Normal batches do not activate these views. An evaluator can explicitly opt in using `EBMKF_DEBUG_PLOTS=1`, optionally `EBMKF_DEBUG_BREAKPOINT=1`/`EBMKF_DEBUG_SAVE_DIR`/`EBMKF_DEBUG_NO_SHOW=1`.

Three windows show: (1) active temperature/OHC, posterior/retrospective states, a blind replay with the inferred multiplier fixed to one, active AOD/solar/eCO2/aerosol inputs and SW/LW/storage components; (2) temperature, OHC and TOA innovation distributions with observation mixtures; (3) the isolated inferred SW correction in W m−2 plus the dimensionless multiplier departure. Temperature benchmark uses the original input (with the same baseline alignment), rather than the volcanically adjusted assimilation input. Blind replay uses the active post-wrapper parameters, including any existing model-specific exchange adjustment, and is labelled accordingly.

The separate correction is ΔSW(k)=[a(k−1)−1] × SW_base(T(k−1),forcing(k)), using the filtered state or the retrospective state respectively. It holds the previous state and current-year forcings fixed while changing only the extra multiplier to one. It is neither a temperature difference between separately evolved trajectories nor the prescribed volcanic preprocessing correction. Positive values add absorbed shortwave. Filtered correction corresponds exactly to the contribution in the year's physical prediction; retrospective correction is a reconstruction using the full-record smoothed states. No prediction exists at initialization, so first-year flux/correction fields are NaN; multiplier fields are end-of-year states. These timing conventions are stated on the figure and in the JSON report.

`--save-dir` writes three PNGs, a diagnostic CSV (both corrections, multipliers, radiation, temperature/OHC and active inputs) and a parameter/definition JSON. No Results/Results2 caches, calibration tables, or legacy extras logs/PNGs are overwritten by this diagnostic command. Existing regular-evaluation behavior is preserved. The detailed historical `plot_historical_diagnostics(save_dir,show)` Python function remains available and also includes the isolated correction figure; core `EBMKF_Nicklas4.py` retains its standalone raw-input plotting.

Ta4 defaults were matched to the live author-tuned 2a: aerosol percentile 85, observational eCO2 percentile 17, fdbkA 0.00, volcanic target 20 ± 5 W m−2 per raw AOD. They reside independently in EBMKF_Nicklas4.py; MPI/NorESM ta2 tuning is untouched. The reference-preserving volcanic map is already active in both core and wrapper paths. Existing ta4-specific smoothing, historical ×1.4 AOD adjustment, state structure and observation/noise choices remain in place. The ±5 value remains an uncertainty descriptor, not automatic random sampling. The prior SSP aerosol-percentile fix remains active; no new forcing join offset is applied.

Validation covers actual historical, MPI SSP245, NorESM volcanic and NorESM VolcConst single-case diagnostics. Before/after historical/MPI/NorESM outputs match exactly. Independent calls to the augmented update map with/without the multiplier reproduce the correction to <1e−10 W m−2, the plotted predicted net flux matches xhatminus, and diagnostic reconstruction leaves state/covariance/forcing arrays unchanged. Breakpoint frame/arrays were checked. Defaults match 2a; the local volcanic derivative remains −20 and reference transmission is unchanged. These checks do not regenerate scientific evaluation/calibration outputs or the dashboard. Backups and examples: `/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-v11-migration-2026-10-17/ta4-future-diagnostics/`.


### Ta4 TOA intervals, initialization spike and delayed-update audit — 17 October 2026

The active component plot's bottom-right energy panel and separate TOA view now include the actual adjusted TOA input and two assigned Gaussian 95% ranges: darker nominal input ±1.96 sqrt(9 TOA_var), and lighter effective input ±1.96 sqrt(9 TOA_var + N cM2[2,2]). These are the model's observation-noise assumptions, not published instrument confidence intervals or posterior/forecast ranges. The effective variance matches the forward filter's assigned observation covariance; state/prediction uncertainty is not added. Nominal/effective sigmas and effective 95% endpoints are included in diagnostic CSV/JSON exports.

An audit of the historical OHC innovation mixture identifies the near-zero blue spike as the first prediction, 1851. Its observation sigma is broad (64.0977 ZJ), but dividing by the initial 2922.7681 ZJ innovation/forecast sigma yields width0.02193. It carries 1/175 of the mixture mass. All later standardized component widths are >=0.798945. The original full mixture is retained; this is initialization, not a new tightly constrained ocean observation or an Argo transition.

The filtered temperature increases 0.095740°C between 1999 and 2001. Exact K×innovation decomposition shows the largest observation term is temperature, but its input was first raised by the inherited delayed volcanic preprocessing: +0.290438°C in2000 and +0.234974°C in2001. `erf_trapez[j]=(TOA_correction[j]+2*TOA_correction[j+1])/3` is applied through `temps[9:] -= erf_trapez[:-7]/(heatCp-Cs)/2`, so the target year's temperature receives corrections from 8–9 years earlier (Pinatubo-era 1991–1993 here). This line is present in git HEAD. The updated non-ocean capacity denominator2.6661062814 versus HEAD's5.261 nearly doubles this heuristic's temperature adjustment at fixed TOA correction. Direct OHC/TOA gain contributions to the 2000/2001 temperature updates are much smaller; the TOA source/noise changes in2001 also coincide with changing gains. No scientific/preprocessing change was made during this audit. Review the timing and amplitude of this inherited adjustment before further tuning or a full scientific rerun.

Complete inputs, four-state previous/predicted/posterior values, all observation gains/contributions, source-year mapping and same-prior/gain counterfactuals are in `/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-v11-migration-2026-10-17/ta4-observation-audit/audit-report.md`; CSV/JSON tables and explanatory plots are alongside it. Updated historical diagnostic plots are under `refreshed-historical/`. Interval exports were checked against assigned forward noise, including an independent covariance reconstruction using S−A K S; filtered temperatures remain unchanged. No scientific caches, manuscript figures, parameters, sources or dashboard data were modified.


The redundant standalone Net TOA figure has been removed from active ta4 diagnostics. The combined energy panel retains the TOA input and assigned uncertainty bands. New runs open/export components, innovations and the separate inferred-correction figure. Previously saved standalone TOA PNGs are historical artifacts.


## Ta4 model-specific starting parameters — 17 October 2026

`Methods/43_Forcing_Based/2_Kalman/EBMKF_ta4_method.py` now groups its four future tuning settings in `FUTURE_PARAMETER_SETTINGS` near the top of the file. MPI-ESM copies the current author-tuned `EBMKF_Nicklas2b.py` starts: `gad=1.1`, `fdbkA=0.40`, aerosol percentile85 and eCO2 percentile50. Both NorESM volcanic options copy `EBMKF_Nicklas2c.py`: `gad=0.85`, `fdbkA=0.25`, aerosol percentile65 and eCO2 percentile35. These are explicit copied starts for review; changing ta2 later does not silently change ta4. The core's existing 2020-OHC adjustment can subsequently change `gad`; saved ta4 diagnostic JSON includes both the starting settings and its post-adjustment value.

Historical and heads/tails calls keep the author's core defaults (`gad=0.6`, `fdbkA=0`, aerosol85, eCO2 17). Baseline restoration now restores both percentiles as well as feedback and exchange, preventing a preceding future call from contaminating those historical settings. No historical preprocessing, ensemble OHC/TOA design, model capacities, inferred-SW correction or noise tuning was changed in this update. The aerosol choice affects the observational prefix and matching SSP curve; eCO2 retains the existing observational-reference/model-derived future construction, including its normalization and NorESM's0.93 factor.

Five member-0 future checks (MPI SSP126/245/370 and NorESM Volc/VolcConst) pass. Historical and heads/tails member-0 means, standard errors and forcing arrays are exactly unchanged before/after the patch and after interleaved future runs. These targeted checks are not a complete ensemble recalibration: review the author workspace `outputs/kalman-v11-migration-2026-10-17/ta4-ensemble-settings/` comparison plot/CSV and then freeze the configuration before restarting the full scientific rebuild. The preceding historical calibration and future pilot are provisional because they predate the final configuration; dependent bulk controllers are stopped. No website publication was performed.
