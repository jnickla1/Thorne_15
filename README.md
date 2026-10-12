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

The active setting is `regen = 1` in [hist_heads_tails_evaluation_script.py](hist_heads_tails_evaluation_script.py). It reruns the explicitly selected method lists and overwrites their result pickles. `regen = 0` reads cached results. `regen = 2` runs missing methods and can retain old results for methods already present; it does not refresh a correction to an existing method.

Each experiment/member produces:

- `Results/results{experiment}{member}.pickle`: current and retrospective method arrays.
- `Results/headstails_statistics_{experiment}{member}.csv`: member-level statistics and crossing summaries.
- `Results/headstails_method_fineprobs_{experiment}{member}.npy`: member-level exceedance probabilities.

There should be 100 of each type per experiment: **400 pickles, 400 statistics CSVs, and 400 probability arrays** for this workflow. Logs are in `logs_histens/run_0.log`, `run_5.log`, …, `run_95.log`. Wait for all batches to finish before processing or copying results. Check logs for exceptions and check file timestamps: counts alone do not prove that existing files were successfully refreshed. The wrapper can continue to another scenario after a Python failure, so an apparently finished wrapper is not sufficient evidence that all four scenarios succeeded.

### 2. Refresh the paper probability figures and rate CSVs

After the member jobs have finished, remain at the repository root:

```sh
MPLBACKEND=Agg python3 plot_heads_tails_probs.py

for experiment in histens_rpi_real histens_rpi_anthro histens_satcal_real histens_satcal_anthro; do
    MPLBACKEND=Agg python3 histens_rates.py "$experiment"
done
```

The probability script reads the refreshed member CSVs and NPYs and writes:

- `baseline_1850-1900_heads_tails_probplot.png`
- `baseline_1981-2010_heads_tails_probplot.png`

The rate script reads the refreshed pickles and writes `rates_histens/rates_histens_{experiment}.csv`. Both are needed for a complete paper/dashboard refresh. `MPLBACKEND=Agg` permits saving figures without opening plot windows; omit it when you want interactive Matplotlib windows.

`hist_evaluation_script_rates.py` is a separate historical-observation rate analysis. It is not a required step for assembling the heads/tails dashboard.

### 3. Export the dashboard data from the refreshed Results directory

```sh
python3 climate-dashboard/source/histens_rates_dashboard.py all \
    --source-root "$PWD" \
    --outdir climate-dashboard/data
```

This is the dashboard analysis entry point. It reads the root `Results/`, the newly updated rate CSVs, the saved paper probability arrays, and the calibration files in `Results2/`. It generates four dataset JSONs plus headline, crossing, rate, probability, and member-level CSV exports.

For this normal refresh, **omit `--results-dir`** so the analysis uses the newly regenerated root `Results/`. The previously created `climate-dashboard/data/corrected_results/` copies were a bridge for testing the FaIR retrospective corrections before the root caches were refreshed; they must not silently override a later full regeneration or new year's data.

The earlier discrepancy between the original 0.335652°C/decade PIVW rate and the dashboard's 0.318823°C/decade came entirely from old versus corrected FaIR retrospective caches. With the same 2025 inputs and calibration weights, the refreshed PIVW rate should be approximately **0.318823°C/decade**, with combined 1σ uncertainty **0.054156°C/decade** and a 95% range of **0.212677–0.424970°C/decade**. These are checks for this dataset and correction, not constants to carry into another year.

### 4. Render the dashboard figure library and assemble the website

```sh
python3 climate-dashboard/source/plot_heads_tails_probs_dashboard.py
python3 climate-dashboard/source/prepare_website.py
```

The first command reads `climate-dashboard/data/` and writes the main figure in PNG/PDF/SVG plus 52 variant PNGs in `climate-dashboard/figures/`. The second copies the exported figures and data into `climate-dashboard/site/dist/assets/` and creates the selected-method CSV downloads. It checks that the complete 52-figure library exists. Both commands derive their default directories from their own script locations, so they can be called from the repository root.

### 5. Review and publish the website exports

The website is static HTML/JavaScript. Browser controls, zoom, hover values, and downloads use exported data; visitors do not run the research Python pipeline.

Preview the assembled site from the repository root:

```sh
python3 -m http.server 8787 --directory climate-dashboard/site/dist
```

Open `http://localhost:8787`. Check the main PIVW headline, both calibrations, realized/anthropogenic views, raw/smoothed backgrounds, figure view, and a data download. Stop the preview with Control-C when finished.

On the author's Mac the separate GitHub Pages checkout is `/Users/JohnMatthew/Downloads/climate-state-dashboard`. Copy the complete assembled website into that checkout:

```sh
cp -R climate-dashboard/site/dist/. /Users/JohnMatthew/Downloads/climate-state-dashboard/
cd /Users/JohnMatthew/Downloads/climate-state-dashboard
git status --short
git add app.js index.html styles.css assets vendor .nojekyll
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

The submodules include upstream implementations and exploratory methods. Being present as a submodule does not mean that a method is included in the current dashboard. The authoritative lists for the heads/tails runner are `sel_methods_list_real` and `sel_methods_list_anthro` in `hist_heads_tails_evaluation_script.py`; combination subsets are defined in `plot_fut_results.py` and the dashboard export script.

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

- `hist_fitprob.py` and `hist_evaluation_script.py`: observation file version, selected folders, calibration choices, plotting endpoints, and dated comparison tables.
- `hist_heads_tails_evaluation_script.py`: selected method lists, `regen`, evaluation windows, and the exclusive upper endpoints in `evalmaxs` (currently 2026 for the 2025 endpoint).
- `histens_rates.py`: `last_year` is currently 2025.
- FaIR/GWI/Green's-function readers: 2025 filenames, array-length branches, bias corrections, member mapping, and time boundaries.
- `climate-dashboard/source/histens_rates_dashboard.py`: probability windows, observations, method lists, and the seven non-crossing members' quartic treatment. A later year's set of non-crossing members can differ.
- `climate-dashboard/source/plot_heads_tails_probs_dashboard.py` and `climate-dashboard/site/dist/app.js` / `index.html`: explicit 2025 headlines, line endpoints, threshold windows, and other display labels.

Use a targeted search such as `rg -n '2025|2026|175|176|last_year|evalmaxs|cut2025|individual_members25'` on the relevant source files. Review matches individually: historical dates and fixed reference periods must not be replaced indiscriminately.

Useful recorded update checkpoints are [HadCRUT5.1 input/code changes, September 2026](https://github.com/jnickla1/Thorne_15/commit/f0672785aeb014663bec678df34c6d3820fce2b0), [the finalized 2025 update, October 2026](https://github.com/jnickla1/Thorne_15/commit/51de5173496524f11dea9618ea1df63382470580), and [the dashboard addition](https://github.com/jnickla1/Thorne_15/commit/bca152245bd181291e184f5268ff3e3b9799c849). These checkpoints follow the December 2025 work; they are examples of what was changed, not a universal patch to apply for every future year.

## Historical calibration and evaluation

Run this part after updating historical inputs or deliberately changing calibration. It is not required solely to adopt the tested retrospective-only FaIR correction while current outputs and calibration remain unchanged.

1. Run `python3 hist_fitprob.py historical` from the root. `historical_regen=True` is required in that script. It runs the chosen methods and writes the `Results2/historical_*` calibration products, including `historical_names_var_scale.csv` used for method variance corrections and PIVW weights.
2. Review `create_combination_method_HHonly.py`: for historical-to-historical calibration use `comparison_type='hh'`. Its current `nmethods_list=[7]` refreshes the realized combination only. To deliberately recalibrate both dashboard combinations, use `[7, 2]` (seven realized methods and two anthropogenic methods), then run `python3 create_combination_method_HHonly.py` and verify both `Results2/comb_real_scales.npy` and `Results2/comb_anthro_scales.npy`. This is a manual setting to review; it has not been changed by updating this README.
3. Run `python3 hist_evaluation_script.py` to regenerate the historical-observation evaluation, `hist_results.pickle`, statistics, and paper figures. Review `historical_regen` and `running_subset`; the current historical evaluator has regeneration enabled and adjusted variances enabled, so calibration must exist first.
4. Run the heads/tails and dashboard refresh above using the new inputs and calibration.

`hist_evaluation_script_rates.py` can additionally be run for the separate historical-observation rate figures. Its own regeneration flag controls whether it reuses `hist_results.pickle`. It is not a substitute for either `histens_rates.py` or the dashboard exporter.

## Future simulations and combination experiments

These are separate paper-validation workflows. They are not required for updating only the current dashboard exports. They need the external model data, their own experiment settings, and usually Oscar/Slurm.

The existing future-method batch wrapper is `run_parBATCH_future.sh`. It evaluates MPI-ESM1-2-LR SSP126/245/370 and NorESM RCP45 volcanic variants using `fut_evaluation_script.py`. Review the cluster account, module loads, environment name, resources, and method-regeneration settings before submitting it:

```sh
sbatch run_parBATCH_future.sh
python3 plot_fut_results.py
```

Run plotting after the submitted jobs have finished. The combination-focused future figures use:

```sh
python3 fut_evaluation_script_elim.py futplotcomb_ESM1-2-LR_SSP126_constVolc -11
python3 fut_evaluation_script_elim.py futplotcomb_NorESM_RCP45_Volc -11
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

In a `histens_rates.py` row whose `level` is 2025, the misleadingly named `year_cross` column now contains the **current warming level**, not a crossing year. The 12 October 2026 reporting correction separates this level from the retrospective series used to fit the rate. `year_cross_sd` uses the current member endpoints as well. Older saved CSVs used retrospective levels; the dashboard exporter accepts both formats and marks the source in `saved_summary.level_source`. The historical crossing rows still retain their existing calculation conventions; use the dashboard's final data exports for its definitive crossing benchmarks.

The rate script also prints a separate per-method/current-combination endpoint comparison: current and retrospective means, paired differences, maximum differences, matching counts at absolute tolerance 1e-6°C, and missing endpoints. These are diagnostic checks, including expected CGWL window/FaIR sample-noise differences, not assertions that stop a run. They are **not added to the CSV**. Individual methods compare their raw retrospective endpoints; combination comparisons include the existing current-proxy fallback used by the rate fits.

Satellite-calibrated views use the source pipeline's 0.7°C reference shift. The raw/smoothed background menu changes only the grey/black ensemble lines, not the calculations. Distribution panels retain missing mass from non-crossing members and share a 0–0.3 probability/year scale, with explicit right-panel magnification. Details are in [climate-dashboard/README.md](climate-dashboard/README.md).

### Which dashboard files to use

Use `climate-dashboard/data/*_corrected_headlines.csv` and each dataset JSON's `headline_summary` for the final dashboard levels and rates. `*_legacy_saved_rates.csv`, `rates_fixed/`, and the original/fixed script copies retain intermediate or historical audit work; they are not additional production entry points. Use `*_crossings_rates.csv` and `*_member_crossings_rates.csv` for the definitive and method-comparison crossing calculations.

`climate-dashboard/source/rebuild_corrected_retrospectives.py` is a one-time correction/regression helper, not a routine annual-update command. It compares the edited FaIR readers with retained backups, guards current outputs, and creates separate corrected pickle copies. Its 2025 assumptions and baselines need review before reuse with a later year. Once the root `Results/` has been fully regenerated, the regular dashboard exporter should read that root directly. `export_fair_same_cut_lag4.py` is a superseded diagnostic: the adopted lagged input is the full earlier-cutoff FaIR matrix from Oscar.

### Common causes of mismatches

- **Method files were edited but old pickles are still being read.** Regenerate the affected root caches or explicitly select the verified corrected cache directory during a correction audit.
- **A method has its own stale intermediate files.** Check its reader/generator, including GAM, quartic/spline, and Bayesian change-point outputs, plus externally supplied FaIR/GWI data. The author cleared the first three cache families for the HadCRUT5.1 update.
- **Only some members or scenarios finished.** Check logs and timestamps before regenerating aggregate CSVs or publishing.
- **CSV and pickle vintages differ.** The dashboard exporter checks the saved levels against current outputs or retained legacy retrospective endpoints; regenerate the four rate CSVs after refreshing pickles.
- **The site still has old exports.** Run `prepare_website.py`, copy the complete `site/dist/` into the Pages checkout, then commit and push that checkout. For JavaScript changes, update the script version query in `index.html` to refresh browser caches.
- **New-year arrays or forecasts are too short.** Review external exports, explicit reader branches, probability-window endpoints, and method-specific intermediates before rerunning the full evaluation.

### Existing result and note links

The working [historical comparison](https://docs.google.com/spreadsheets/d/10izz9VruI9L1pNT3pwKLlNPVhzrvGRdYk3VxvdQ1es8/edit), [MPI-ESM1-2-LR comparison](https://docs.google.com/spreadsheets/d/1eWAeL1HHHSqyL1YF2IYaQwgMQuh8Y1RJ/edit), [NorESM comparison](https://docs.google.com/spreadsheets/d/1gHNtpZ4MVIw_NYp62kjtuCHZ2kZcWYwo/edit), and [progress/citation notes](https://docs.google.com/spreadsheets/d/1iShljXO2rmPHpPjPkBGbPSwZc_7XVtCU32O3sn6sjTc/edit) are retained for context. They are author-maintained working documents; confirm their data vintage when comparing them with regenerated local outputs.
