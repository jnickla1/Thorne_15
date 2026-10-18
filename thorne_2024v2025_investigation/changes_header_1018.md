# 17 Oct → 18 Oct refresh of the 2024→2025 investigation

The 17 Oct note was built on the 10 Oct paper state (git `f36ba84b` and the Oct-10 heads/tails backup). This refresh swaps in the outputs regenerated on 18 Oct 2026 as the "2025 vintage". The 2024 vintage (git `93793643`, `snap24/`, `hist_results_93793643.pickle`) is unchanged.

## Inputs used (all checked; see `inputs_1018_manifest.csv` for path, mtime and md5 of 1,427 files)

| Input | Count | mtime (18 Oct unless stated) |
|---|---|---|
| `Results/resultshistens_{rpi,satcal}_{real,anthro}{0..99}.pickle` | 4 × 100 | 00:01–01:21 |
| `Results/headstails_method_fineprobs_histens_*`, `headstails_statistics_histens_*` | 4 × 100 each | 00:01–01:21; copied to `snap25_1018/` |
| `hist_results.pickle` | 1 | 11:29; copied to `hist_results_1018.pickle` |
| `261018_HadCRUT5.1_var_adjusted/current_methods_statistics_261018True.csv` | 1 | 11:31 |
| `Results/resultsfut_ESM1-2-LR_SSP245/SSP370_constVolc*`, `NorESM_RCP45_VolcConst*` | 50, 50, 60 | 01:07–01:12 |
| `Results/resultsfut_ESM1-2-LR_SSP126_constVolc*` | 50 | **12:21–12:22** (see below) |
| `rates_histens/rates_histens_histens_*`, `climate-dashboard/data/*_corrected_headlines.csv` | 4 + 4 | 01:21 |
| `Results2/historical_names_var_scale.csv`, `comb_*_scales.npy` | 3 | 17 Oct 23:58–23:59 |

**Completeness and stability.**
- All member counts are complete: 100 per heads/tails experiment, 50 per MPI future and 60 for NorESM.
- The SSP1-2.6 future pickles were rewritten at 12:21, about 11 h after the other futures. That was probably the open IDLE session; the running `future-refresh` combination job only writes to its own folder.
- They have not changed since 12:22. Their 1850–2014 values differ from SSP2-4.5/3-7.0 for FaIR, GWI and CGWL by up to 0.04 °C. This is inherent to those methods: SSP2-4.5 and SSP3-7.0 also differ from each other by up to 0.05 °C. The temperature-only methods agree exactly.
- The SSP1-2.6 column is therefore treated as part of the 18 Oct vintage. If a further SSP1-2.6 rerun is planned, that column should be refreshed again.

**Observational inputs are unchanged, so (a) is unchanged by construction.**
- `Hist_heads_tails/sst_pseudo*.csv` is from 8 Sep, and `HadCRUT.5.1.csv` from 21 Sep.
- `toyKFmodelData9c.csv` (17 Oct 01:23) gives the same eruption-flag years.
- So `compare_inputs.py` is not rerun. Its numbers (0.036 °C cooler 2024, 0.030 °C cooler 2015–2024, SD 0.084 vs 0.093, SD ratio 1.10, HadCRUT5.0.2→5.1 −0.005 °C) all stand.

## What actually changed in the model outputs (`diff_old_new_1018.py`, `diff_old_new_1018.txt`)

**Heads/tails ensemble: only EBMKF_ta4 changed materially.**
- Its current estimate rose by +0.008 °C in 2024 and +0.009 °C in 2025 (ensemble mean), with up to 0.12 °C differences in earlier years.
- FaIR_nonat_unB changed by < 0.0003 °C.
- Every other component method (lowess, GAM, Kal_flexLin, FaIR_comb, GWI_tot/anthro, CGWL10y) is bit-identical, including `cent20y`.

**Single HadCRUT5.1 run.**
- EBMKF_ta4 rose by +0.013 °C (2024) and +0.012 °C (2025). EBMKF_ta and the `removeMEI_volc_*` methods also changed, but these are not in the PIVW.
- FaIR_nonat changed by ≤ 0.0002 °C in 2024/25.

**PIVW weights.** These were recomputed as 1/`err_var100` from `Results2/historical_names_var_scale.csv`. They reproduce the rates-CSV PIVW to 2×10⁻⁶ °C.

| | lowess | GAM | Kal_flexLin | FaIR_comb | **EBMKF_ta4** | GWI_tot | CGWL10y | FaIR_nonat | GWI_anthro |
|---|---|---|---|---|---|---|---|---|---|
| 17 Oct note (paper) | 7.5 | 7.4 | 7.9 | 24.6 | **16.4** | 25.7 | 10.5 | 69 | 31 |
| 18 Oct (Results2 23:58) | 7.7 | 7.6 | 8.1 | 25.3 | **14.1** | 26.4 | 10.8 | 69.0 | 31.0 |

**The anthropogenic side is effectively unchanged.**
- PIVW anthro 2025 is 1.33664 vs 1.33664.
- All anthropogenic decomposition numbers are unchanged to the stated precision.

## Changes to the note's text (rounded values that change)

| Where in note | 17 Oct text | 18 Oct text |
|---|---|---|
| Headline table, realised 2025 | 1.36 [1.23–1.53] °C | **1.37** [1.23–1.53] °C (1.3655) |
| Headline table, realised change | −0.04 | **−0.03** (unrounded −0.041 vs 1.406) |
| Headline: realised P(>1.5 °C) 2025 | 5.9 % | **6.4 %** (still "very unlikely") |
| Headline: dashboard values | 1.364 [1.230, 1.528] | **1.366 [1.229, 1.532]** |
| Decomposition (ii) realised | 1.336, ≈ −0.025 | **1.337, ≈ −0.024** |
| Decomposition (iii) realised | 1.364, +0.029 (min +0.021) | **1.366, +0.028 (min +0.020)** |
| (ii) PIVW spread-adjusted shift | 0.026; 1.336 + 0.026 = 1.362 | **0.022; 1.337 + 0.022 = 1.359** (direct estimate 1.361) |
| (iii) non-FaIR 2025-step percentile | 97th–99th | **97th–>99th** (EBMKF_ta4 now 99.7–99.8) |
| (c) EBMKF_ta4 spread-adjusted shift | +0.017 °C | **−0.006 °C** |
| (c) single-run same-year differences | "All other methods: within ±0.008 °C" | true except **EBMKF_ta4: +0.012 (2015), +0.033 (2020), +0.017 (2023)**, i.e. a (c) method change |
| (b) weight-sensitivity table | 1.364 / 1.364 / 1.364 / 1.373 | **1.366 / 1.365 / 1.365 / 1.374** |
| (b) EBMKF_ta4 weight | "10.9 %–19.0 %" range noted | now **14.1 %** (inside that range) |
| (d) EBMKF_ta4 row | 5 (0) · 4 (0.3) · 41 (23) · 18 (7) · 5 (1) · 14 (4) | **8 (0) · 12 (1) · 42 (29) · 25 (14) · 10 (4) · 19 (8)** |
| (d) PIVW realised row | 11 (0) · 8 (1) · 40 (18) · 13 (4) · 3 (0.5) · 10 (2) | **13 (0) · 9 (1) · 40 (19) · 15 (5) · 3 (0.7) · 10 (2)** |
| (d) PIVW anthro row | 5 (0) · … · 5 (1) | **3 (0)** · … · **6 (1)** |
| (d) FaIR_nonat row, SSP3-7.0 / NorESM | 2 (0.3) / 7 (2) | **1 (0.3) / 8 (2)** |
| (d) Kal_flexLin, NorESM | 21 (10) | **22 (10)** |
| (d) "Futures from git f36ba84b; SSP245 rerun not included" | — | futures now 18 Oct, SSP245 rerun included |
| Caveat: realised inversion | 1.362 | **1.359** |
| Caveat: EBMKF_ta4 weight | ≈16 % | **≈14 %** |
| Thorne reply | "to 1.36 °C in 2025" | **"to 1.37 °C"** (see note) |

**Unchanged at the stated precision:**
- Decomposition (i) for both sides: 1.361 / 1.355, −0.045 / +0.015. The direct estimate is 1.3609 with the 18 Oct weights and 1.3606 with the paper weights.
- All anthropogenic rows.
- FaIR numbers: +0.047 / +0.064 shifts, 30th–35th percentiles, 69 % weight.
- Temperature-only shifts: +0.026 to +0.035.
- GWI/CGWL shifts.
- The cent20y check: 0.5 % of member-years, min −0.003 °C, no HadCRUT5.1 decreases 1966–2015.
- (d) conclusion: no realised member fell in 2025; one FaIR_nonat member at −0.0004.
- The "4–20 %" and "< 0.01 °C" ranges.
- The "95–100 % spurious" range for warming futures (new minimum 94.7 %, old 93.8 %; both already round to about 95).
- The 1.408 legacy-convention check.
- The HadCRUT5.0-run PIVW of 1.360.

**Rounding hazard.** The new realised headline, 1.36555, rounds to **1.37** at 2 dp. If the manuscript still says 1.36 °C, it needs updating to match the 18 Oct outputs, or else the reply should keep three decimals. This is flagged here only; the manuscript is handled in another chat.

## Full table (every number the note uses)

Generated by `build_changes_1018.py`, which also writes `changes_1017_vs_1018.csv`. "Rounded text changes?" compares the old and new values at the precision used in the note.
