# U250: 2024/2025 GMST statistics for the Figure 26 paragraph

Date: 2026-10-18. Everything was read-only outside `u250_ensemble_check/`.

## Bottom line

- **The 10,000-member ensemble is not stored anywhere I could reach, and it can't be regenerated from public inputs alone.** Details are in §1.
- **The current 2024 text matches the 100-member ensemble used in the analysis exactly**, to 2 decimal places: median 1.53 °C, 95 % range [1.38, 1.72] °C, mean 1.53 °C. That match suggests the published "10,000-member" numbers were actually computed from the 100-member subset, or that the two agree very closely. **The 10k values themselves can't be checked independently.**
- **For 2025, the 100-member ensemble gives** median 1.44 (1.435) °C, 95 % range [1.29, 1.61] °C, mean 1.44 °C.
- **WMO 2025 headline:** 1.43 °C in the full State of the Global Climate 2025 report (23 Mar 2026), and 1.44 ± 0.13 °C in the 14 Jan 2026 press release. The ensemble agrees with both.

## 1. Where the 10,000-member ensemble is (and isn't)

**Local submodule**, `Hist_heads_tails/global_temperature_merge`:
- `git submodule status` shows `-dabb185…`, i.e. not initialised as a submodule. The directory does hold a checkout of `dabb185` ("consistent formatting", 5 Sep 2025).
- It contains code only: no `Output/`, and no data.

**Upstream clone**, `u250_ensemble_check/upstream` (https://github.com/jjk-code-otter/global_temperature_merge):
- `main` is at `5cd1c2b` (6 Jan 2026). It hard-codes `end_year=2024` in `make_meta_ensemble.py`.
- `origin/merger` is at `7dd8aff` (15 Sep 2026, "correct erroneous dataset names in trees"). This is the 2025-capable code: `end_year = 2025`, with family trees revised "for resubmission".
- Neither branch commits ensemble output. The only output file is `Output*/Frechet.txt`.

**Elsewhere on this machine:**
- `DATADIR` is unset.
- There is no `ManagedData` directory.
- Spotlight finds no other `sst_pseudo*` files.

**The pipeline never writes the 10k ensemble to disk anyway.** In `final_clustered` (10000 → 100, balanced k-means, seed 178525), `MetaEnsembleFactory.make_meta_ensemble` replaces `output_dataset` with the clustered result before returning, so `to_csv` writes only the 100 members. The only experiment that writes a 10k file is `basic` (10000 members, no thinning, seed 8794), which would produce `Output/basic/sst_pseudo.csv`. That draw is equivalent but not the identical 10k ensemble.

**`Hist_heads_tails/metaensemblefactory.py`** is identical to the `main`/`dabb185` version. Against `merger`, it differs only in the default `end_year` (2024 vs 2025) and a size-check guard. Which tree and commit produced the local 2025 `sst_pseudo.csv` isn't recorded locally.

**Regeneration is blocked by inputs, not runtime.** The `merger` `sst_pseudo` tree needs:

| Input | Availability |
|---|---|
| Kadow et al. ensemble (`Kadow_Kadow_ensemble`, in tails) | **Request only**, from C. Kadow |
| Cooper et al. 2025 200-member ensemble (`cooper.py`: "received by email") | **Request only** |
| ERA5 ensemble | Request from A. Simmons, or CDS download (needs a CDS account) |
| JRA-3Q | Request from JMA, or the UCAR alternative translator |
| DCENT-I, DCENT_MLE v1.2, HadCRU_MLE v1.4, HadCRUT5, GloSAT, Berkeley Earth Hires, CMST3, CMA-GMST, COBE-STEMP3, NOAA v5.0 ensemble / v5.1 / v6.1, GISTEMP | Public downloads |

The trees on `merger` have also changed since the paper's run (for example, DCENT_I, CMA_GMST and Cooper were added). A regeneration from that branch might therefore not reproduce the paper's ensemble anyway.

**What to ask John Kennedy / Bruce Calvert for:**
- either the pre-clustering 10k array from the 2025 `final_clustered` `sst_pseudo` run;
- or `Output/basic/sst_pseudo.csv` from the same code and data version;
- or simply the 2024/2025 median, mean and 2.5/97.5 percentiles from that 10k ensemble.

## 2. Statistics (1850–1900 baseline)

**Source:** 100-member ensemble `Hist_heads_tails/sst_pseudo.csv`. This is a **fallback**: the 10k ensemble was not available.
- Baseline check: each member's 1850–1900 mean is within 1.5e-5 of 0.
- Percentiles use numpy's default (linear interpolation). With 100 members, each tail percentile rests on about 2–3 members.

| Year | Median | Mean | 2.5 % | 97.5 % | SD |
|---|---|---|---|---|---|
| 2024 | 1.533 | 1.532 | 1.382 | 1.716 | 0.085 |
| 2025 | 1.435 | 1.439 | 1.293 | 1.614 | 0.080 |

- The 2025 median is exactly 1.4350, the average of the two central members. It sits on a rounding boundary: 1.44 with round-half-up, 1.43 with round-half-even.
- These agree with `sst_pseudo_summary.csv`: mean 1.532 / 1.438, SD 0.084 / 0.079.
- `to_2024/sst_pseudo_rpi.csv` and `sst_pseudo_satcal.csv` are other variants, not relevant here.

**Is "1.53 °C [1.38, 1.72], mean 1.53" confirmed for 2024?** Yes, against the 100-member subset, which reproduces it exactly. It is *not* independently confirmed against the 10k ensemble. Balanced k-means is meant to preserve the distribution, so the 10k values are expected to be within a few hundredths, most likely in the tails.

## 3. Suggested replacement text

Option A keeps 2024 for the WMO comparison and adds 2025.

> The median estimate of GMST/GSAT for the 2024 annual anomaly relative to 1850–1900 is 1.53 °C (95 % range [1.38, 1.72] °C), with a mean of 1.53 °C. This is very close to the headline estimate of 1.55 °C in the WMO State of the Global Climate 2024 (WMO, 2025), and the ranges are consistent. For 2025, the final year of the analysis, the median is 1.44 °C (95 % range [1.29, 1.61] °C) and the mean 1.44 °C. This is again consistent with the WMO estimate of about 1.43 °C (WMO State of the Global Climate 2025, WMO, 2026).

Notes:
- If the authors keep "10,000-member" as the source, they should first confirm these numbers with Kennedy/Calvert (§1). Otherwise, say "the 100-member ensemble used in the analysis".
- WMO 2025 sources:
  - [State of the Global Climate 2025 page](https://wmo.int/publication-series/state-of-global-climate/state-of-global-climate-2025), 23 Mar 2026: "about 1.43 °C above the 1850-1900 average".
  - [Press release, 14 Jan 2026](https://wmo.int/news/media-centre/wmo-confirms-2025-was-one-of-warmest-years-record): 1.44 °C ± 0.13 °C from eight datasets.
- Cite the full report (1.43 °C) and add WMO (2026) to the references.

Option B is minimal: keep the 2024 sentence as it is and add "(2024 is used here for comparison with WMO, 2025)". This leaves the 2025 value unstated, which is weaker given that Section 5 runs to 2025.

## 4. Addendum: how different is the 10k ensemble likely to be?

The method is `Dataset.cluster_ensemble` with default `centres=False` (merger `7dd8aff`).

**How the 100 members are chosen:**
- Balanced k-means splits the 10,000 trajectories into exactly 100 clusters of exactly 100 members each.
- The trajectories are the 1981–2010 anomalies for 1850–2025.
- One member is then drawn uniformly at random from each cluster.

**Consequence: the 100 is a stratified random sample of the 10k.**
- Every member has an inclusion probability of exactly 1/100, so the 100-member mean and quantiles are unbiased for the 10k values.
- The sampling error is at most that of a simple random sample of 100.

**The stratification helps for exactly the quantity in the paper.** In the 100-member ensemble:
- 61 % of the squared distance that k-means works on comes from 1850–1900 (91 % from 1850–1950); only 2 % comes from 1996–2025.
- A member's 2024 anomaly relative to 1850–1900 correlates at r = 0.88 with its 1850–1900 baseline offset.
- 77 % of the variance in that 2024 anomaly is explained by the first three principal components of the trajectories (93 % by the first ten).
- So the clusters are strata in nearly the variable that sets the spread of the 2024 and 2025 anomalies relative to 1850–1900.

**Expected 1σ differences between the 100-member and 10k values (2024; 2025 similar):**

| Statistic | Upper bound (simple random sample, n = 100, σ ≈ 0.084) | Likely, with stratification |
|---|---|---|
| Mean | 0.008 | ~0.004 |
| Median | 0.011 | ~0.005–0.007 |
| 2.5 / 97.5 % | 0.022 | ~0.01–0.02 |

**Conclusion:**
- The 10k mean and median for 2024 and 2025 almost certainly round to the same values (1.53; 1.44, or 1.43 for the 2025 median, which is on the rounding boundary).
- The 95 % range bounds could plausibly differ by about 0.01–0.02 °C at 2 dp, for example [1.37–1.39, 1.70–1.73].
- This is an estimate, not a verification.

## 5. Update (2026-10-19): 10,000-member ensemble now available

**Source.** Kennedy & Calvert (2026), figshare, https://doi.org/10.6084/m9.figshare.34640442.v1 (CC BY 4.0).
- The zip `meta_ensemble_20260915.zip` was unpacked locally to `~/Downloads/meta_ensemble_20260915/`. It is kept outside the repo and read-only.
- `Output/final_clustered/sst_pseudo.csv` is **byte-identical** to `Hist_heads_tails/sst_pseudo.csv`. This confirms that the paper's 100 members come from the `final_clustered` run (10,000 → 100, balanced k-means, seed 178525).
- The pre-clustering parent of those 100 is not saved. The 10,000-member file is `Output/basic/sst_pseudo.csv`:
  - It uses the same `sst_pseudo` tree, the same 1981–2010 30-year overlap, the same code version and the same run date, but a different seed (8794).
  - Two independent 10k draws differ by only about 0.001 °C in the mean, so this file represents the 10k population.
  - The 1850–1900 baseline offset is ≤1.8e-5 for every member.

**2024 and 2025 statistics (1850–1900 baseline):**

| Year | Ensemble | Median | Mean | 2.5 % | 97.5 % | SD |
|---|---|---|---|---|---|---|
| 2024 | 10k (`basic`) | 1.528 | 1.536 | 1.388 | 1.720 | 0.085 |
| 2024 | 100 (`final_clustered`) | 1.533 | 1.532 | 1.382 | 1.716 | 0.085 |
| 2025 | 10k (`basic`) | 1.427 | 1.434 | 1.280 | 1.616 | 0.084 |
| 2025 | 100 (`final_clustered`) | 1.435 | 1.439 | 1.293 | 1.614 | 0.080 |

- **Agreement between the 100-member and 10k ensembles.** The 0.013 °C figure applies to **2024 and 2025 only**. Over all 176 years, at the 2.5/5/25/50/75/95/97.5 % percentiles:

  | Statistic | Largest gap | Typical (median) gap |
  |---|---|---|
  | Mean | 0.013 °C (1949) | 0.002 °C |
  | Median | 0.025 °C | 0.004 °C |
  | 25 % and 75 % | up to 0.026 °C | 0.004 °C |
  | 2.5 % and 97.5 % | up to 0.048 °C (97.5 %, 1865) | 0.011 °C |

  - The tails differ by more than 0.02 °C in about 1 year in 5, and by more than 0.03 °C in 7–8 years.
  - Scaled by the sampling error of a simple random sample of 100, the gaps have a median of 0.44σ. Only 0.6 % of the gaps exceed 2σ, against about 5 % expected for a simple random sample. The stratified sampling is therefore doing better than random sampling, as argued in §4.
- **2024 text check.** The text "median 1.53 [1.38, 1.72], mean 1.53" matches the 100-member values. Against the 10k values:
  - the median (1.53) and the upper bound (1.72) still hold;
  - the lower bound rounds to 1.39 (1.388);
  - the mean rounds to 1.54 (1.536).
- The 10k distribution is slightly right-skewed: in 2024 the mean exceeds the median by 0.008 °C.

**|mean − median| over 1850–2025 (176 years):**

| | ≤ 0.01 °C | ≤ 0.02 °C | Max | Max since 1950 |
|---|---|---|---|---|
| 10k | 143 (81 %) | 170 (97 %) | 0.027 (1910) | 0.015 |
| 100 | 122 (69 %) | 169 (96 %) | 0.032 (1854) | 0.017 |

- In the 10k ensemble, the years exceeding 0.02 °C are 1850, 1909–1911, 1914 and 1944.
- "Generally within 0.02 °C" is true for both ensembles.

**Suggested Figure 26 text using the 10k values:**

> The median estimate of GMST/GSAT for the 2024 annual anomaly relative to 1850–1900 is 1.53 °C (95 % range [1.39, 1.72] °C), with a mean of 1.54 °C. This is very close to the headline estimate of 1.55 °C in the WMO State of the Global Climate 2024 (WMO, 2025), and the ranges are consistent. For 2025, the median is 1.43 °C (95 % range [1.28, 1.62] °C) and the mean 1.43 °C, matching the WMO estimate of about 1.43 °C (WMO, 2026).

**Data availability: suggested sentence**

> The full 10,000-member meta-ensemble, the 100-member balanced k-means subset used in the analysis, and the processed input series are available from Kennedy and Calvert (2026; https://doi.org/10.6084/m9.figshare.34640442.v1), and the code used to generate them is at https://github.com/jjk-code-otter/global_temperature_merge.

Reference (ESSD style): Kennedy, J. and Calvert, B.: Meta ensemble of global mean temperatures from Thorne et al. 2026, figshare [data set], https://doi.org/10.6084/m9.figshare.34640442.v1, 2026.

The code reference should point to a specific commit or tag. The 2025 run is on the `merger` branch, at `7dd8aff` (15 Sep 2026), not on `main`. A Zenodo release would give it a DOI.
