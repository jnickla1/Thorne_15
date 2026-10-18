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

### Headline

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| Realised PIVW 2025 level (rates/dashboard) | 1.364 | 1.366 | +0.0013 | **yes** | 17 Oct text 1.364; rounds to 1.37 (was 1.36) at 2 dp |
| Realised PIVW 2025, 2 dp as in table | 1.36 | 1.37 | +0.010 | **yes** | 1.36555 → 1.37: conflicts with manuscript 1.36 |
| Realised 2025 interval lower | 1.230 | 1.229 | -0.0014 | **yes** | [1.23–1.53] at 2 dp unchanged |
| Realised 2025 interval upper | 1.528 | 1.532 | +0.0037 | **yes** |  |
| Realised change vs preprint (2 dp) | -0.04 | -0.03 | +0.010 | **yes** | 1.37 − 1.40; unrounded 1.3655 − 1.406 = −0.041 |
| Anthropogenic PIVW 2025 level | 1.337 | 1.337 | -0.0004 | no | bit-identical inputs except FaIR_nonat (<3e-4) |
| Anthropogenic 2025 interval | [1.209, 1.496] | [1.209, 1.496] |  | no |  |
| Realised mean-member P(>1.5) 2025 (%) | 5.9 | 6.4 | +0.52 | **yes** | likelihood wording unchanged (very unlikely <10 %, extremely unlikely <5 %) |
| Realised P(>1.5) 2024, 2024 vintage (%) | 7.5 | 7.5 | +0.00 | no |  |
| Realised P(>1.5) 2024, new vintage (%) | 2.5 | 2.9 | +0.39 | **yes** |  |
| Anthropogenic mean-member P(>1.5) 2025 (%) | 2.1 | 2.1 | -0.00 | no | likelihood wording unchanged (very unlikely <10 %, extremely unlikely <5 %) |
| Anthropogenic P(>1.5) 2024, 2024 vintage (%) | 7.1 | 7.1 | +0.00 | no |  |
| Anthropogenic P(>1.5) 2024, new vintage (%) | 0.9 | 0.9 | -0.00 | no |  |

### Weights (PIVW %)

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| lowess1dt36wnc | 7.5 | 7.7 | +0.17 | **yes** | new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58) |
| GAM_AR1 | 7.4 | 7.6 | +0.25 | **yes** | new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58) |
| Kal_flexLin | 7.9 | 8.1 | +0.18 | **yes** | new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58) |
| FaIR_comb_unB | 24.6 | 25.3 | +0.67 | **yes** | new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58) |
| EBMKF_ta4 | 16.4 | 14.1 | -2.28 | **yes** | new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58) |
| GWI_tot_CGWL | 25.7 | 26.4 | +0.72 | **yes** | new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58) |
| CGWL10y_sfUKCP | 10.5 | 10.8 | +0.29 | **yes** | new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58) |
| FaIR_nonat_unB | 69.0 | 69.0 | +0.02 | no | new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58) |
| GWI_anthro_CGWL | 31.0 | 31.0 | -0.02 | no | new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58) |

### Decomposition

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| Realised: (i) 2024 vintage, current endpoints | 1.361 | 1.361 | +0.0003 | no | 18 Oct weights; paper weights give 1.3606 (realised) |
| Realised: (i) step vs reported 1.406 | -0.045 | -0.045 | +0.0003 | no |  |
| Realised: (ii) level, new vintage at 2024 | 1.336 | 1.337 | +0.0015 | **yes** |  |
| Realised: (ii) step | -0.025 | -0.024 | +0.0012 | **yes** |  |
| Realised: (iii) 2025 level | 1.364 | 1.366 | +0.0014 | **yes** |  |
| Realised: (iii) one-year step (mean) | 0.029 | 0.028 | -0.0001 | **yes** |  |
| Realised: (iii) min member step | 0.021 | 0.020 | -0.0007 | **yes** |  |
| Realised: (i)+(ii) shortfall | -0.07 | -0.07 | +0.001 | no |  |
| Anthropogenic: (i) 2024 vintage, current endpoints | 1.355 | 1.355 | +0.0001 | no | 18 Oct weights; paper weights give 1.3606 (realised) |
| Anthropogenic: (i) step vs reported 1.34 | 0.015 | 0.015 | +0.0001 | no |  |
| Anthropogenic: (ii) level, new vintage at 2024 | 1.319 | 1.319 | -0.0000 | no |  |
| Anthropogenic: (ii) step | -0.036 | -0.036 | -0.0001 | no |  |
| Anthropogenic: (iii) 2025 level | 1.337 | 1.337 | -0.0000 | no |  |
| Anthropogenic: (iii) one-year step (mean) | 0.017 | 0.017 | +0.0000 | no |  |
| Anthropogenic: (iii) min member step | 0.010 | 0.010 | +0.0000 | no |  |
| Anthropogenic: (i)+(ii) shortfall | -0.02 | -0.02 | -0.000 | no |  |

### Per-method 2025 step (heads/tails mean)

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| lowess1dt36wnc | 0.041 | 0.041 | +0.0000 | no | min member +0.0300; frac<0 0.00 |
| GAM_AR1 | 0.027 | 0.027 | +0.0000 | no | min member +0.0049; frac<0 0.00 |
| Kal_flexLin | 0.041 | 0.041 | +0.0000 | no | min member +0.0254; frac<0 0.00 |
| FaIR_comb_unB | 0.013 | 0.013 | +0.0000 | no | min member +0.0038; frac<0 0.00 |
| EBMKF_ta4 | 0.038 | 0.040 | +0.0014 | **yes** | min member +0.0325; frac<0 0.00 |
| GWI_tot_CGWL | 0.030 | 0.030 | +0.0000 | no | min member +0.0226; frac<0 0.00 |
| CGWL10y_sfUKCP | 0.031 | 0.031 | +0.0000 | no | min member +0.0115; frac<0 0.00 |
| FaIR_nonat_unB | 0.011 | 0.011 | +0.0000 | no | min member -0.0004; frac<0 0.01 |
| GWI_anthro_CGWL | 0.032 | 0.032 | +0.0000 | no | min member +0.0255; frac<0 0.00 |

### Per-method 2024 level, new vintage

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| lowess1dt36wnc | 1.329 | 1.329 | +0.0000 | no |  |
| GAM_AR1 | 1.359 | 1.359 | +0.0000 | no |  |
| Kal_flexLin | 1.352 | 1.352 | +0.0000 | no |  |
| FaIR_comb_unB | 1.355 | 1.355 | +0.0000 | no |  |
| EBMKF_ta4 | 1.323 | 1.330 | +0.0078 | **yes** |  |
| GWI_tot_CGWL | 1.303 | 1.303 | +0.0000 | no |  |
| CGWL10y_sfUKCP | 1.367 | 1.367 | +0.0000 | no |  |
| FaIR_nonat_unB | 1.327 | 1.327 | -0.0000 | no |  |
| GWI_anthro_CGWL | 1.302 | 1.302 | +0.0000 | no |  |

### 2025 step percentile (after annual drop ≤ thr)

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| lowess1dt36wnc, thr -0.08 | 99 | 99 | +0.0 | no | pct_all 84 → 84; median step +0.0071 → +0.0071 |
| lowess1dt36wnc, thr -0.1 | 99 | 99 | +0.0 | no | pct_all 84 → 84; median step +0.0049 → +0.0049 |
| GAM_AR1, thr -0.08 | 98 | 98 | +0.0 | no | pct_all 56 → 56; median step -0.0084 → -0.0084 |
| GAM_AR1, thr -0.1 | 99 | 99 | +0.0 | no | pct_all 56 → 56; median step -0.0109 → -0.0109 |
| Kal_flexLin, thr -0.08 | 98 | 98 | +0.0 | no | pct_all 75 → 75; median step +0.0007 → +0.0007 |
| Kal_flexLin, thr -0.1 | 98 | 98 | +0.0 | no | pct_all 75 → 75; median step -0.0018 → -0.0018 |
| FaIR_comb_unB, thr -0.08 | 30 | 30 | +0.0 | no | pct_all 25 → 25; median step +0.0238 → +0.0238 |
| FaIR_comb_unB, thr -0.1 | 30 | 30 | +0.0 | no | pct_all 25 → 25; median step +0.0233 → +0.0233 |
| EBMKF_ta4, thr -0.08 | 98 | 100 | +1.3 | **yes** | pct_all 88 → 84; median step +0.0126 → +0.0013 |
| EBMKF_ta4, thr -0.1 | 99 | 100 | +0.7 | **yes** | pct_all 88 → 84; median step +0.0114 → -0.0005 |
| GWI_tot_CGWL, thr -0.08 | 98 | 98 | +0.0 | no | pct_all 76 → 76; median step +0.0083 → +0.0083 |
| GWI_tot_CGWL, thr -0.1 | 98 | 98 | +0.0 | no | pct_all 76 → 76; median step +0.0068 → +0.0068 |
| CGWL10y_sfUKCP, thr -0.08 | 97 | 97 | +0.0 | no | pct_all 65 → 65; median step -0.0114 → -0.0114 |
| CGWL10y_sfUKCP, thr -0.1 | 98 | 98 | +0.0 | no | pct_all 65 → 65; median step -0.0172 → -0.0172 |
| FaIR_nonat_unB, thr -0.08 | 34 | 34 | +0.1 | no | pct_all 26 → 26; median step +0.0155 → +0.0155 |
| FaIR_nonat_unB, thr -0.1 | 35 | 35 | +0.1 | no | pct_all 26 → 26; median step +0.0148 → +0.0148 |
| GWI_anthro_CGWL, thr -0.08 | 99 | 99 | +0.0 | no | pct_all 82 → 82; median step +0.0075 → +0.0075 |
| GWI_anthro_CGWL, thr -0.1 | 99 | 99 | +0.0 | no | pct_all 82 → 82; median step +0.0067 → +0.0067 |
| PIVW_real, thr -0.08 | 95 | 96 | +0.5 | **yes** | pct_all 68 → 67; median step +0.0079 → +0.0063 |
| PIVW_real, thr -0.1 | 97 | 97 | +0.2 | no | pct_all 68 → 67; median step +0.0061 → +0.0047 |
| PIVW_anthro, thr -0.08 | 65 | 65 | +0.0 | no | pct_all 43 → 43; median step +0.0125 → +0.0125 |
| PIVW_anthro, thr -0.1 | 71 | 71 | +0.0 | no | pct_all 43 → 43; median step +0.0118 → +0.0118 |

### Spread-adjusted equivalent shift

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| rpi_real lowess1dt36wnc | 0.030 | 0.030 | +0.0000 | no | unadjusted +0.0420 → +0.0420 |
| rpi_real GAM_AR1 | 0.026 | 0.026 | +0.0000 | no | unadjusted +0.0339 → +0.0339 |
| rpi_real Kal_flexLin | 0.035 | 0.035 | +0.0000 | no | unadjusted +0.0434 → +0.0434 |
| rpi_real FaIR_comb_unB | 0.047 | 0.047 | +0.0000 | no | unadjusted +0.0557 → +0.0557 |
| rpi_real EBMKF_ta4 | 0.017 | -0.006 | -0.0231 | **yes** | unadjusted +0.0307 → +0.0086 |
| rpi_real GWI_tot_CGWL | 0.024 | 0.024 | +0.0000 | no | unadjusted +0.0388 → +0.0388 |
| rpi_real CGWL10y_sfUKCP | 0.056 | 0.056 | +0.0000 | no | unadjusted +0.0612 → +0.0612 |
| rpi_real PIVW_comb | 0.026 | 0.022 | -0.0038 | **yes** | unadjusted +0.0384 → +0.0347 |
| rpi_anthro FaIR_nonat_unB | 0.064 | 0.064 | +0.0000 | no | unadjusted +0.0741 → +0.0741 |
| rpi_anthro GWI_anthro_CGWL | 0.021 | 0.021 | +0.0000 | no | unadjusted +0.0345 → +0.0345 |
| rpi_anthro PIVW_comb | 0.051 | 0.051 | +0.0000 | no | unadjusted +0.0620 → +0.0620 |
| Realised PIVW implied 2024-vintage level (caveat) | 1.362 | 1.359 | -0.0024 | **yes** |  |
| Anthro PIVW implied 2024-vintage level (caveat) | 1.370 | 1.370 | -0.0000 | **yes** |  |

### HadCRUT single run (5.1 run − 5.0.2 run, same year)

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| lowess1dt36wnc vintage_eff_2024 | -0.002 | -0.002 | +0.0000 | no |  |
| lowess1dt36wnc vintage_eff_2015 | -0.005 | -0.005 | +0.0000 | no |  |
| lowess1dt36wnc year_step_25 | 0.033 | 0.033 | +0.0000 | no |  |
| GAM_AR1 vintage_eff_2024 | 0.008 | 0.008 | +0.0000 | no |  |
| GAM_AR1 vintage_eff_2015 | 0.005 | 0.005 | +0.0000 | no |  |
| GAM_AR1 year_step_25 | 0.017 | 0.017 | +0.0000 | no |  |
| Kal_flexLin vintage_eff_2024 | -0.001 | -0.001 | +0.0000 | no |  |
| Kal_flexLin vintage_eff_2015 | -0.005 | -0.005 | +0.0000 | no |  |
| Kal_flexLin year_step_25 | 0.031 | 0.031 | +0.0000 | no |  |
| FaIR_comb_unB vintage_eff_2024 | -0.013 | -0.013 | +0.0000 | no |  |
| FaIR_comb_unB vintage_eff_2015 | 0.019 | 0.019 | +0.0000 | no |  |
| FaIR_comb_unB year_step_25 | 0.012 | 0.012 | +0.0000 | no |  |
| EBMKF_ta4 vintage_eff_2024 | — | — |  | no |  |
| EBMKF_ta4 vintage_eff_2015 | -0.009 | 0.012 | +0.0207 | **yes** |  |
| EBMKF_ta4 year_step_25 | 0.033 | 0.032 | -0.0012 | **yes** |  |
| GWI_tot_CGWL vintage_eff_2024 | 0.001 | 0.001 | +0.0000 | no |  |
| GWI_tot_CGWL vintage_eff_2015 | -0.003 | -0.003 | +0.0000 | no |  |
| GWI_tot_CGWL year_step_25 | 0.023 | 0.023 | +0.0000 | no |  |
| CGWL10y_sfUKCP vintage_eff_2024 | -0.001 | -0.001 | +0.0000 | no |  |
| CGWL10y_sfUKCP vintage_eff_2015 | -0.004 | -0.004 | +0.0000 | no |  |
| CGWL10y_sfUKCP year_step_25 | 0.019 | 0.019 | +0.0000 | no |  |
| FaIR_nonat_unB vintage_eff_2024 | -0.038 | -0.038 | -0.0002 | no |  |
| FaIR_nonat_unB vintage_eff_2015 | 0.021 | 0.020 | -0.0002 | **yes** |  |
| FaIR_nonat_unB year_step_25 | 0.010 | 0.010 | +0.0004 | no |  |
| GWI_anthro_CGWL vintage_eff_2024 | -0.005 | -0.005 | +0.0000 | no |  |
| GWI_anthro_CGWL vintage_eff_2015 | -0.007 | -0.007 | +0.0000 | no |  |
| GWI_anthro_CGWL year_step_25 | 0.024 | 0.024 | +0.0000 | no |  |
| PIVW_real vintage_eff_2024 | -0.003 | -0.003 | +0.0000 | no |  |
| PIVW_real vintage_eff_2015 | 0.002 | 0.005 | +0.0032 | **yes** |  |
| PIVW_real year_step_25 | 0.022 | 0.022 | -0.0005 | no |  |
| PIVW_anthro vintage_eff_2024 | -0.028 | -0.028 | -0.0001 | no |  |
| PIVW_anthro vintage_eff_2015 | 0.012 | 0.012 | -0.0001 | no |  |
| PIVW_anthro year_step_25 | 0.014 | 0.014 | +0.0002 | no |  |
| EBMKF_ta4 vintage_eff_2023 (2024 run ends 2023) | 0.007 | 0.017 | +0.0099 | **yes** | breaks "all others within ±0.008" |
| EBMKF_ta4 vintage_eff_2020 | 0.006 | 0.033 | +0.0275 | **yes** |  |
| PIVW_real level 2024, 2025-vintage run | 1.360 | 1.361 | +0.0015 | **yes** |  |
| PIVW_real level 2024, 2024 run (text 1.360) | 1.360 | 1.360 | +0.0000 | no |  |

### (b) Weight sensitivity: 2025 PIVW level

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| Paper weights | 1.364 | 1.366 | +0.0013 | **yes** | EBMKF_ta4 weight 14.1 % |
| Paper weights (17 Oct) on new levels | 1.364 | 1.366 | +0.0015 | **yes** | EBMKF_ta4 weight 16.4 % |
| 1/RMS² Dec 2025 HadCRUT5.0 | 1.364 | 1.365 | +0.0013 | **yes** | EBMKF_ta4 weight 19.5 % |
| 1/RMS² HadCRUT5.1 (Oct 10 False → Oct 18 var_adjusted True) | 1.364 | 1.365 | +0.0013 | **yes** | EBMKF_ta4 weight 17.6 % |
| 1/RMS² HadCRUT5.1 Oct 10 False | 1.364 | 1.365 | +0.0014 | **yes** | EBMKF_ta4 weight 19.9 % |
| Equal weights | 1.373 | 1.374 | +0.0008 | **yes** | EBMKF_ta4 weight 14.3 % |
| Range of 2025 PIVW across non-equal weight sets | 0.000 | 0.000 | +0.0004 | no | "≤0.005 °C" claim holds |

### (d) % of non-volcanic years with Δ<0

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| CGWL10y_sfUKCP | HadCRUT5.1 ≥1960 | 26 | 26 | +0.0 | no |  |
| EBMKF_ta4 | HadCRUT5.1 ≥1960 | 5 | 8 | +2.6 | **yes** |  |
| FaIR_comb_unB | HadCRUT5.1 ≥1960 | 8 | 8 | +0.0 | no |  |
| FaIR_nonat_unB | HadCRUT5.1 ≥1960 | 11 | 11 | +0.0 | no |  |
| GAM_AR1 | HadCRUT5.1 ≥1960 | 21 | 21 | +0.0 | no |  |
| GWI_anthro_CGWL | HadCRUT5.1 ≥1960 | 11 | 11 | +0.0 | no |  |
| GWI_tot_CGWL | HadCRUT5.1 ≥1960 | 16 | 16 | +0.0 | no |  |
| Kal_flexLin | HadCRUT5.1 ≥1960 | 16 | 16 | +0.0 | no |  |
| PIVW_anthro | HadCRUT5.1 ≥1960 | 5 | 3 | -2.6 | **yes** |  |
| PIVW_real | HadCRUT5.1 ≥1960 | 11 | 13 | +2.6 | **yes** |  |
| lowess1dt36wnc | HadCRUT5.1 ≥1960 | 8 | 8 | +0.0 | no |  |
| CGWL10y_sfUKCP | heads/tails ≥1960 | 23 | 23 | +0.0 | no |  |
| EBMKF_ta4 | heads/tails ≥1960 | 4 | 12 | +7.3 | **yes** |  |
| FaIR_comb_unB | heads/tails ≥1960 | 7 | 7 | +0.0 | no |  |
| FaIR_nonat_unB | heads/tails ≥1960 | 7 | 7 | +0.0 | no |  |
| GAM_AR1 | heads/tails ≥1960 | 18 | 18 | +0.0 | no |  |
| GWI_anthro_CGWL | heads/tails ≥1960 | 7 | 7 | +0.0 | no |  |
| GWI_tot_CGWL | heads/tails ≥1960 | 9 | 9 | +0.0 | no |  |
| Kal_flexLin | heads/tails ≥1960 | 13 | 13 | +0.0 | no |  |
| PIVW_anthro | heads/tails ≥1960 | 4 | 4 | -0.0 | no |  |
| PIVW_real | heads/tails ≥1960 | 8 | 9 | +0.9 | **yes** |  |
| lowess1dt36wnc | heads/tails ≥1960 | 8 | 8 | +0.0 | no |  |
| CGWL10y_sfUKCP | MPI_SSP126 | 41 | 41 | +0.0 | no |  |
| EBMKF_ta4 | MPI_SSP126 | 41 | 42 | +1.9 | **yes** |  |
| FaIR_comb_unB | MPI_SSP126 | 42 | 42 | +0.0 | no |  |
| FaIR_nonat_unB | MPI_SSP126 | 42 | 42 | -0.1 | no |  |
| GAM_AR1 | MPI_SSP126 | 43 | 43 | +0.0 | no |  |
| GWI_anthro_CGWL | MPI_SSP126 | 32 | 32 | +0.0 | no |  |
| GWI_tot_CGWL | MPI_SSP126 | 33 | 33 | +0.0 | no |  |
| Kal_flexLin | MPI_SSP126 | 44 | 44 | -0.0 | no |  |
| PIVW_anthro | MPI_SSP126 | 40 | 40 | -0.1 | no |  |
| PIVW_real | MPI_SSP126 | 40 | 40 | +0.4 | no |  |
| lowess1dt36wnc | MPI_SSP126 | 42 | 42 | +0.0 | no |  |
| CGWL10y_sfUKCP | MPI_SSP245 | 28 | 28 | +0.0 | no |  |
| EBMKF_ta4 | MPI_SSP245 | 18 | 25 | +6.7 | **yes** |  |
| FaIR_comb_unB | MPI_SSP245 | 6 | 6 | +0.0 | no |  |
| FaIR_nonat_unB | MPI_SSP245 | 8 | 8 | +0.1 | no |  |
| GAM_AR1 | MPI_SSP245 | 27 | 27 | +0.0 | no |  |
| GWI_anthro_CGWL | MPI_SSP245 | 0 | 0 | +0.0 | no |  |
| GWI_tot_CGWL | MPI_SSP245 | 1 | 1 | +0.0 | no |  |
| Kal_flexLin | MPI_SSP245 | 26 | 26 | +0.0 | no |  |
| PIVW_anthro | MPI_SSP245 | 6 | 6 | +0.0 | no |  |
| PIVW_real | MPI_SSP245 | 13 | 15 | +1.1 | **yes** |  |
| lowess1dt36wnc | MPI_SSP245 | 20 | 20 | +0.0 | no |  |
| CGWL10y_sfUKCP | MPI_SSP370 | 17 | 17 | +0.0 | no |  |
| EBMKF_ta4 | MPI_SSP370 | 5 | 10 | +5.3 | **yes** |  |
| FaIR_comb_unB | MPI_SSP370 | 1 | 1 | +0.0 | no |  |
| FaIR_nonat_unB | MPI_SSP370 | 2 | 1 | -0.0 | **yes** |  |
| GAM_AR1 | MPI_SSP370 | 14 | 14 | +0.0 | no |  |
| GWI_anthro_CGWL | MPI_SSP370 | 0 | 0 | +0.0 | no |  |
| GWI_tot_CGWL | MPI_SSP370 | 0 | 0 | +0.0 | no |  |
| Kal_flexLin | MPI_SSP370 | 12 | 12 | +0.0 | no |  |
| PIVW_anthro | MPI_SSP370 | 1 | 1 | +0.0 | no |  |
| PIVW_real | MPI_SSP370 | 3 | 3 | +0.4 | no |  |
| lowess1dt36wnc | MPI_SSP370 | 6 | 6 | +0.0 | no |  |
| CGWL10y_sfUKCP | NorESM_RCP45_VolcConst | 23 | 23 | -0.0 | no |  |
| EBMKF_ta4 | NorESM_RCP45_VolcConst | 14 | 19 | +5.6 | **yes** |  |
| FaIR_comb_unB | NorESM_RCP45_VolcConst | 7 | 7 | +0.0 | no |  |
| FaIR_nonat_unB | NorESM_RCP45_VolcConst | 7 | 8 | +0.2 | **yes** |  |
| GAM_AR1 | NorESM_RCP45_VolcConst | 23 | 23 | +0.0 | no |  |
| GWI_anthro_CGWL | NorESM_RCP45_VolcConst | 0 | 0 | +0.0 | no |  |
| GWI_tot_CGWL | NorESM_RCP45_VolcConst | 0 | 0 | +0.0 | no |  |
| Kal_flexLin | NorESM_RCP45_VolcConst | 21 | 22 | +0.1 | **yes** |  |
| PIVW_anthro | NorESM_RCP45_VolcConst | 5 | 6 | +0.2 | **yes** |  |
| PIVW_real | NorESM_RCP45_VolcConst | 10 | 10 | +0.6 | no |  |
| lowess1dt36wnc | NorESM_RCP45_VolcConst | 15 | 15 | +0.0 | no |  |

### (d) % of non-volcanic years with Δ<−0.01

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| CGWL10y_sfUKCP | HadCRUT5.1 ≥1960 | 21 | 21 | +0.0 | no |  |
| EBMKF_ta4 | HadCRUT5.1 ≥1960 | 0.0 | 0.0 | +0.00 | no |  |
| FaIR_comb_unB | HadCRUT5.1 ≥1960 | 0.0 | 0.0 | +0.00 | no |  |
| FaIR_nonat_unB | HadCRUT5.1 ≥1960 | 0.0 | 0.0 | +0.00 | no |  |
| GAM_AR1 | HadCRUT5.1 ≥1960 | 16 | 16 | +0.0 | no |  |
| GWI_anthro_CGWL | HadCRUT5.1 ≥1960 | 3 | 3 | +0.0 | no |  |
| GWI_tot_CGWL | HadCRUT5.1 ≥1960 | 3 | 3 | +0.0 | no |  |
| Kal_flexLin | HadCRUT5.1 ≥1960 | 3 | 3 | +0.0 | no |  |
| PIVW_anthro | HadCRUT5.1 ≥1960 | 0.0 | 0.0 | +0.00 | no |  |
| PIVW_real | HadCRUT5.1 ≥1960 | 0.0 | 0.0 | +0.00 | no |  |
| lowess1dt36wnc | HadCRUT5.1 ≥1960 | 3 | 3 | +0.0 | no |  |
| CGWL10y_sfUKCP | heads/tails ≥1960 | 15 | 15 | +0.0 | no |  |
| EBMKF_ta4 | heads/tails ≥1960 | 0 | 1 | +0.7 | **yes** |  |
| FaIR_comb_unB | heads/tails ≥1960 | 0.1 | 0.1 | +0.00 | no |  |
| FaIR_nonat_unB | heads/tails ≥1960 | 0.1 | 0.1 | +0.00 | no |  |
| GAM_AR1 | heads/tails ≥1960 | 10 | 10 | +0.0 | no |  |
| GWI_anthro_CGWL | heads/tails ≥1960 | 0.6 | 0.6 | +0.00 | no |  |
| GWI_tot_CGWL | heads/tails ≥1960 | 0.7 | 0.7 | +0.00 | no |  |
| Kal_flexLin | heads/tails ≥1960 | 5 | 5 | +0.0 | no |  |
| PIVW_anthro | heads/tails ≥1960 | 0.0 | 0.0 | +0.00 | no |  |
| PIVW_real | heads/tails ≥1960 | 0.7 | 0.8 | +0.08 | **yes** |  |
| lowess1dt36wnc | heads/tails ≥1960 | 2 | 2 | +0.0 | no |  |
| CGWL10y_sfUKCP | MPI_SSP126 | 29 | 29 | +0.0 | no |  |
| EBMKF_ta4 | MPI_SSP126 | 23 | 29 | +5.3 | **yes** |  |
| FaIR_comb_unB | MPI_SSP126 | 10 | 10 | +0.0 | no |  |
| FaIR_nonat_unB | MPI_SSP126 | 13 | 13 | +0.1 | no |  |
| GAM_AR1 | MPI_SSP126 | 31 | 31 | +0.0 | no |  |
| GWI_anthro_CGWL | MPI_SSP126 | 0.1 | 0.1 | +0.00 | no |  |
| GWI_tot_CGWL | MPI_SSP126 | 0.1 | 0.1 | +0.00 | no |  |
| Kal_flexLin | MPI_SSP126 | 30 | 30 | +0.0 | no |  |
| PIVW_anthro | MPI_SSP126 | 9 | 9 | +0.0 | no |  |
| PIVW_real | MPI_SSP126 | 18 | 19 | +1.2 | **yes** |  |
| lowess1dt36wnc | MPI_SSP126 | 25 | 25 | +0.0 | no |  |
| CGWL10y_sfUKCP | MPI_SSP245 | 18 | 18 | +0.0 | no |  |
| EBMKF_ta4 | MPI_SSP245 | 7 | 14 | +6.8 | **yes** |  |
| FaIR_comb_unB | MPI_SSP245 | 0.7 | 0.7 | +0.00 | no |  |
| FaIR_nonat_unB | MPI_SSP245 | 1 | 1 | -0.0 | no |  |
| GAM_AR1 | MPI_SSP245 | 16 | 16 | +0.0 | no |  |
| GWI_anthro_CGWL | MPI_SSP245 | 0.0 | 0.0 | +0.00 | no |  |
| GWI_tot_CGWL | MPI_SSP245 | 0.0 | 0.0 | +0.00 | no |  |
| Kal_flexLin | MPI_SSP245 | 16 | 16 | +0.0 | no |  |
| PIVW_anthro | MPI_SSP245 | 0.8 | 0.8 | -0.05 | no |  |
| PIVW_real | MPI_SSP245 | 4 | 5 | +0.6 | **yes** |  |
| lowess1dt36wnc | MPI_SSP245 | 8 | 8 | +0.0 | no |  |
| CGWL10y_sfUKCP | MPI_SSP370 | 11 | 11 | +0.0 | no |  |
| EBMKF_ta4 | MPI_SSP370 | 1 | 4 | +3.4 | **yes** |  |
| FaIR_comb_unB | MPI_SSP370 | 0.2 | 0.2 | +0.00 | no |  |
| FaIR_nonat_unB | MPI_SSP370 | 0.3 | 0.3 | -0.02 | no |  |
| GAM_AR1 | MPI_SSP370 | 8 | 8 | +0.0 | no |  |
| GWI_anthro_CGWL | MPI_SSP370 | 0.0 | 0.0 | +0.00 | no |  |
| GWI_tot_CGWL | MPI_SSP370 | 0.0 | 0.0 | +0.00 | no |  |
| Kal_flexLin | MPI_SSP370 | 5 | 5 | +0.0 | no |  |
| PIVW_anthro | MPI_SSP370 | 0.1 | 0.1 | +0.00 | no |  |
| PIVW_real | MPI_SSP370 | 0.5 | 0.7 | +0.26 | **yes** |  |
| lowess1dt36wnc | MPI_SSP370 | 2 | 2 | +0.0 | no |  |
| CGWL10y_sfUKCP | NorESM_RCP45_VolcConst | 14 | 14 | +0.0 | no |  |
| EBMKF_ta4 | NorESM_RCP45_VolcConst | 4 | 8 | +4.0 | **yes** |  |
| FaIR_comb_unB | NorESM_RCP45_VolcConst | 1.0 | 1.0 | +0.00 | no |  |
| FaIR_nonat_unB | NorESM_RCP45_VolcConst | 1 | 2 | +0.1 | **yes** |  |
| GAM_AR1 | NorESM_RCP45_VolcConst | 12 | 12 | +0.0 | no |  |
| GWI_anthro_CGWL | NorESM_RCP45_VolcConst | 0.0 | 0.0 | +0.00 | no |  |
| GWI_tot_CGWL | NorESM_RCP45_VolcConst | 0.0 | 0.0 | +0.00 | no |  |
| Kal_flexLin | NorESM_RCP45_VolcConst | 10 | 10 | +0.0 | no |  |
| PIVW_anthro | NorESM_RCP45_VolcConst | 0.9 | 0.9 | +0.04 | no |  |
| PIVW_real | NorESM_RCP45_VolcConst | 2 | 2 | +0.3 | no |  |
| lowess1dt36wnc | NorESM_RCP45_VolcConst | 4 | 4 | +0.0 | no |  |

### (d) mean size of decreases

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| CGWL10y_sfUKCP | HadCRUT5.1 ≥1960 | -0.020 | -0.020 | +0.0000 | no | p05 -0.0339 → -0.0339; spurious 89 → 89 % |
| EBMKF_ta4 | HadCRUT5.1 ≥1960 | -0.003 | -0.003 | +0.0004 | no | p05 +0.0031 → -0.0004; spurious 100 → 100 % |
| FaIR_comb_unB | HadCRUT5.1 ≥1960 | -0.004 | -0.004 | +0.0000 | no | p05 -0.0008 → -0.0008; spurious 67 → 67 % |
| FaIR_nonat_unB | HadCRUT5.1 ≥1960 | -0.002 | -0.002 | -0.0001 | no | p05 -0.0008 → -0.0009; spurious 75 → 75 % |
| GAM_AR1 | HadCRUT5.1 ≥1960 | -0.016 | -0.016 | +0.0000 | no | p05 -0.0197 → -0.0197; spurious 86 → 86 % |
| GWI_anthro_CGWL | HadCRUT5.1 ≥1960 | -0.007 | -0.007 | +0.0000 | no | p05 -0.0059 → -0.0059; spurious 100 → 100 % |
| GWI_tot_CGWL | HadCRUT5.1 ≥1960 | -0.006 | -0.006 | +0.0000 | no | p05 -0.0077 → -0.0077; spurious 100 → 100 % |
| Kal_flexLin | HadCRUT5.1 ≥1960 | -0.006 | -0.006 | +0.0000 | no | p05 -0.0083 → -0.0083; spurious 100 → 100 % |
| PIVW_anthro | HadCRUT5.1 ≥1960 | -0.001 | -0.002 | -0.0011 | **yes** | p05 +0.0007 → +0.0003; spurious 50 → 100 % |
| PIVW_real | HadCRUT5.1 ≥1960 | -0.005 | -0.004 | +0.0007 | **yes** | p05 -0.0028 → -0.0030; spurious 75 → 80 % |
| lowess1dt36wnc | HadCRUT5.1 ≥1960 | -0.007 | -0.007 | +0.0000 | no | p05 -0.0021 → -0.0021; spurious 67 → 67 % |
| CGWL10y_sfUKCP | heads/tails ≥1960 | -0.019 | -0.019 | +0.0000 | no | p05 -0.0352 → -0.0352; spurious 93 → 93 % |
| EBMKF_ta4 | heads/tails ≥1960 | -0.004 | -0.005 | -0.0003 | **yes** | p05 +0.0008 → -0.0049; spurious 94 → 99 % |
| FaIR_comb_unB | heads/tails ≥1960 | -0.003 | -0.003 | +0.0000 | no | p05 -0.0014 → -0.0014; spurious 61 → 61 % |
| FaIR_nonat_unB | heads/tails ≥1960 | -0.002 | -0.002 | -0.0000 | no | p05 -0.0009 → -0.0010; spurious 66 → 67 % |
| GAM_AR1 | heads/tails ≥1960 | -0.014 | -0.014 | +0.0000 | no | p05 -0.0170 → -0.0170; spurious 93 → 93 % |
| GWI_anthro_CGWL | heads/tails ≥1960 | -0.005 | -0.005 | +0.0000 | no | p05 -0.0021 → -0.0021; spurious 93 → 93 % |
| GWI_tot_CGWL | heads/tails ≥1960 | -0.005 | -0.005 | +0.0000 | no | p05 -0.0037 → -0.0037; spurious 81 → 81 % |
| Kal_flexLin | heads/tails ≥1960 | -0.008 | -0.008 | +0.0000 | no | p05 -0.0091 → -0.0091; spurious 96 → 96 % |
| PIVW_anthro | heads/tails ≥1960 | -0.002 | -0.002 | -0.0000 | no | p05 +0.0004 → +0.0004; spurious 57 → 57 % |
| PIVW_real | heads/tails ≥1960 | -0.005 | -0.005 | -0.0001 | no | p05 -0.0026 → -0.0033; spurious 83 → 85 % |
| lowess1dt36wnc | heads/tails ≥1960 | -0.006 | -0.006 | +0.0000 | no | p05 -0.0027 → -0.0027; spurious 85 → 85 % |
| CGWL10y_sfUKCP | MPI_SSP126 | -0.023 | -0.023 | -0.0000 | no | p05 -0.0458 → -0.0458; spurious 60 → 60 % |
| EBMKF_ta4 | MPI_SSP126 | -0.015 | -0.019 | -0.0044 | **yes** | p05 -0.0281 → -0.0366; spurious 59 → 61 % |
| FaIR_comb_unB | MPI_SSP126 | -0.007 | -0.007 | +0.0000 | no | p05 -0.0136 → -0.0136; spurious 54 → 54 % |
| FaIR_nonat_unB | MPI_SSP126 | -0.008 | -0.008 | -0.0000 | no | p05 -0.0152 → -0.0152; spurious 56 → 56 % |
| GAM_AR1 | MPI_SSP126 | -0.023 | -0.023 | +0.0000 | no | p05 -0.0477 → -0.0477; spurious 62 → 62 % |
| GWI_anthro_CGWL | MPI_SSP126 | -0.003 | -0.003 | -0.0000 | no | p05 -0.0051 → -0.0051; spurious 42 → 42 % |
| GWI_tot_CGWL | MPI_SSP126 | -0.003 | -0.003 | +0.0000 | no | p05 -0.0054 → -0.0054; spurious 44 → 44 % |
| Kal_flexLin | MPI_SSP126 | -0.020 | -0.020 | -0.0000 | no | p05 -0.0402 → -0.0403; spurious 61 → 61 % |
| PIVW_anthro | MPI_SSP126 | -0.007 | -0.007 | -0.0000 | no | p05 -0.0128 → -0.0127; spurious 52 → 52 % |
| PIVW_real | MPI_SSP126 | -0.011 | -0.012 | -0.0005 | **yes** | p05 -0.0221 → -0.0232; spurious 55 → 56 % |
| lowess1dt36wnc | MPI_SSP126 | -0.015 | -0.015 | +0.0000 | no | p05 -0.0291 → -0.0291; spurious 57 → 57 % |
| CGWL10y_sfUKCP | MPI_SSP245 | -0.021 | -0.021 | -0.0000 | no | p05 -0.0381 → -0.0381; spurious 97 → 97 % |
| EBMKF_ta4 | MPI_SSP245 | -0.010 | -0.014 | -0.0041 | **yes** | p05 -0.0141 → -0.0232; spurious 98 → 98 % |
| FaIR_comb_unB | MPI_SSP245 | -0.004 | -0.004 | +0.0000 | no | p05 -0.0011 → -0.0011; spurious 99 → 99 % |
| FaIR_nonat_unB | MPI_SSP245 | -0.005 | -0.005 | -0.0000 | no | p05 -0.0026 → -0.0025; spurious 98 → 98 % |
| GAM_AR1 | MPI_SSP245 | -0.017 | -0.017 | +0.0000 | no | p05 -0.0285 → -0.0285; spurious 98 → 98 % |
| GWI_anthro_CGWL | MPI_SSP245 | -0.003 | -0.003 | +0.0000 | no | p05 +0.0071 → +0.0071; spurious 100 → 100 % |
| GWI_tot_CGWL | MPI_SSP245 | -0.004 | -0.004 | +0.0000 | no | p05 +0.0062 → +0.0062; spurious 100 → 100 % |
| Kal_flexLin | MPI_SSP245 | -0.016 | -0.016 | -0.0000 | no | p05 -0.0261 → -0.0261; spurious 98 → 98 % |
| PIVW_anthro | MPI_SSP245 | -0.005 | -0.005 | -0.0000 | no | p05 -0.0004 → -0.0004; spurious 99 → 99 % |
| PIVW_real | MPI_SSP245 | -0.008 | -0.008 | -0.0005 | no | p05 -0.0084 → -0.0096; spurious 97 → 97 % |
| lowess1dt36wnc | MPI_SSP245 | -0.010 | -0.010 | +0.0000 | no | p05 -0.0144 → -0.0144; spurious 97 → 97 % |
| CGWL10y_sfUKCP | MPI_SSP370 | -0.020 | -0.020 | +0.0000 | no | p05 -0.0246 → -0.0246; spurious 100 → 100 % |
| EBMKF_ta4 | MPI_SSP370 | -0.006 | -0.010 | -0.0040 | **yes** | p05 +0.0005 → -0.0084; spurious 100 → 100 % |
| FaIR_comb_unB | MPI_SSP370 | -0.005 | -0.005 | -0.0000 | no | p05 +0.0094 → +0.0094; spurious 100 → 100 % |
| FaIR_nonat_unB | MPI_SSP370 | -0.006 | -0.006 | -0.0001 | no | p05 +0.0084 → +0.0085; spurious 100 → 100 % |
| GAM_AR1 | MPI_SSP370 | -0.016 | -0.016 | +0.0000 | no | p05 -0.0171 → -0.0171; spurious 100 → 100 % |
| GWI_anthro_CGWL | MPI_SSP370 | -0.003 | -0.003 | +0.0000 | no | p05 +0.0122 → +0.0122; spurious 100 → 100 % |
| GWI_tot_CGWL | MPI_SSP370 | -0.004 | -0.004 | -0.0000 | no | p05 +0.0109 → +0.0109; spurious 100 → 100 % |
| Kal_flexLin | MPI_SSP370 | -0.012 | -0.012 | -0.0000 | no | p05 -0.0110 → -0.0110; spurious 100 → 100 % |
| PIVW_anthro | MPI_SSP370 | -0.004 | -0.004 | -0.0000 | no | p05 +0.0112 → +0.0112; spurious 100 → 100 % |
| PIVW_real | MPI_SSP370 | -0.006 | -0.007 | -0.0009 | **yes** | p05 +0.0047 → +0.0036; spurious 100 → 100 % |
| lowess1dt36wnc | MPI_SSP370 | -0.008 | -0.008 | +0.0000 | no | p05 -0.0018 → -0.0018; spurious 100 → 100 % |
| CGWL10y_sfUKCP | NorESM_RCP45_VolcConst | -0.017 | -0.017 | -0.0000 | no | p05 -0.0264 → -0.0264; spurious 98 → 98 % |
| EBMKF_ta4 | NorESM_RCP45_VolcConst | -0.007 | -0.010 | -0.0027 | **yes** | p05 -0.0073 → -0.0138; spurious 97 → 98 % |
| FaIR_comb_unB | NorESM_RCP45_VolcConst | -0.005 | -0.005 | -0.0000 | no | p05 -0.0020 → -0.0020; spurious 95 → 95 % |
| FaIR_nonat_unB | NorESM_RCP45_VolcConst | -0.006 | -0.006 | +0.0001 | no | p05 -0.0027 → -0.0029; spurious 95 → 96 % |
| GAM_AR1 | NorESM_RCP45_VolcConst | -0.013 | -0.013 | +0.0000 | no | p05 -0.0193 → -0.0193; spurious 98 → 98 % |
| GWI_anthro_CGWL | NorESM_RCP45_VolcConst | -0.004 | -0.004 | +0.0000 | no | p05 +0.0073 → +0.0073; spurious 100 → 100 % |
| GWI_tot_CGWL | NorESM_RCP45_VolcConst | -0.002 | -0.002 | +0.0000 | no | p05 +0.0071 → +0.0071; spurious 100 → 100 % |
| Kal_flexLin | NorESM_RCP45_VolcConst | -0.011 | -0.011 | +0.0000 | no | p05 -0.0157 → -0.0158; spurious 98 → 98 % |
| PIVW_anthro | NorESM_RCP45_VolcConst | -0.005 | -0.005 | +0.0001 | no | p05 -0.0004 → -0.0004; spurious 94 → 95 % |
| PIVW_real | NorESM_RCP45_VolcConst | -0.006 | -0.006 | -0.0003 | no | p05 -0.0042 → -0.0047; spurious 96 → 96 % |
| lowess1dt36wnc | NorESM_RCP45_VolcConst | -0.008 | -0.008 | +0.0000 | no | p05 -0.0087 → -0.0087; spurious 97 → 97 % |

### (d) cent20y monotonicity

| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |
|---|---|---|---|---|---|
| heads/tails member-years 1966–2015 decreasing (%) | 0.5 | 0.5 | +0.00 | no | from step2025_and_checks_1018.log; min step −0.0026 both |
| HadCRUT5.1 cent20y decreasing years 1966–2015 | 0 | 0 | +0.0 | no |  |
