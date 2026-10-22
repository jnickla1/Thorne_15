# Why the current-warming estimate fell slightly from the 2024 to the 2025 update

*Draft supplement section, **21 Oct 2026 refresh (1850–1900 datum)**. Numbers come from:*
- *2024 update: the HadCRUT5.0-based state at git `93793643` (20 Dec 2025, data through 2024). Unchanged.*
- *2025 update: the outputs regenerated on 21 Oct 2026 after the switch to the 1850–1900 datum and the anchored-method datum correction (commits `216d3fde`, `acf2df22`, `2d95aaec`). These are the heads/tails ensemble (02:18–03:35), `hist_results.pickle` and `current_methods_statistics_261021True.csv` (02:16), and the weights (`Results2`, 02:09). Files, mtimes and md5s are in `inputs_1021_manifest.csv`.*
- *The rates and dashboard (`rates_histens/`, `climate-dashboard/data/*_corrected_headlines.csv`) were regenerated at 11:37–11:39 on 21 Oct. They agree with the levels computed here from the pickles to 1e-12 °C.*
- *The future ensembles on the 1850–1900 datum are still rerunning on Oscar. The future columns of the (d) table use the 18 Oct future runs. A constant datum offset does not change year-on-year steps, so these columns are expected to stand.*
- *Table S4 (CSCM) and Figure 17 do not feed these numbers. See `changes_1019_vs_1021.md` for every number that changed.*

*Scripts and output tables are in `thorne_2024v2025_investigation/` (`*_1021.py`, `*_1021.csv`, `baseline_hypothesis_1021.py`).*

## Headline

The 2024 values are from the ESSD preprint abstract (essd-2025-825). The 2025 values are the 21 Oct dashboard export on the 1850–1900 datum: 1.3646 [1.228, 1.530] (realised) and 1.3341 [1.206, 1.493] (anthropogenic). The realised value rounds to 1.36, matching the manuscript. The anthropogenic value rounds to **1.33**, but the manuscript says 1.34.

| | 2024 (preprint) | 2025 (21 Oct outputs) | Change | Expected at ~0.3 °C/decade |
|---|---|---|---|---|
| Realised | 1.40 [1.23–1.58] °C, "unlikely" >1.5 °C | 1.36 [1.23–1.53] °C, "very unlikely" | −0.04 | ≈ +0.03 |
| Anthropogenic | 1.34 [1.18–1.50] °C, "very unlikely" | 1.33 [1.21–1.49] °C, "extremely unlikely" | −0.01 | ≈ +0.03 |

**Where the numbers come from.**
- **2024 values:** the PIVW values in the 2024 rate CSV (`rates_histens/`, git `876b4b8c`/`93793643`), 1.406 and 1.340 °C. The script took each method's last *retrospective* value (`LT_trend[2][-1]`) and used the current output only when no retrospective value existed.
- **2025 values:** the dashboard export of *current* outputs, 1.365 [1.228, 1.530] and 1.334 [1.206, 1.493] °C.
- **Exceedance probabilities:** the mean member P(>1.5 °C) is always computed from current outputs. It fell from 7.5 % to 6.2 % (realised) and from 7.1 % to 1.9 % (anthropogenic).

**Explanations tested.**
- (a) Input data
- (b) Method weights
- (c) Internal method changes
- (d) Methods decreasing while tracking the 20-yr mean
- (e) The change of pre-industrial datum from 1850–1899 to 1850–1900

The shortfall against expectation comes from a reporting correction, (a), and a smaller FaIR contribution under (c). The retuned EBMKF_ta4 is a small (c) effect in the opposite direction. Neither (b) nor (d) contributes materially. (e) lowers both years equally, by 0.001 °C (realised) and 0.003 °C (anthropogenic). It is a level correction, not a cause of the small 2024→2025 step.

## Decomposition (PIVW, heads/tails ensemble, 1850–1900 baseline)

| Step | Realised level (°C) | Change | Anthropogenic level (°C) | Change |
|---|---|---|---|---|
| 2024 update as reported (retrospective endpoints) | 1.406 | | 1.340 | |
| (i) Same run, current endpoints | ≈1.361 | ≈ −0.045 | ≈1.355 | ≈ +0.015 |
| (ii) 2025 inputs and method versions, still year 2024 (old datum) | 1.337 | ≈ −0.023 | 1.319 | ≈ −0.036 |
| (ii-e) Anchored-method datum correction, 1850–1899 → 1850–1900 | 1.336 | −0.001 | 1.317 | −0.003 |
| (iii) One more year, 2024 → 2025 | 1.365 | +0.029 | 1.334 | +0.017 |

### (i) Reporting correction (12 Oct 2026)

FaIR_comb_unB's retrospective branch had a mixing-weight error (weights summing to 7/6). This put its retrospective endpoint about 0.185 °C above its current value. That excess was measured twice:
- 1.552 vs 1.368 in the 2025 ensemble
- 1.581 vs 1.393 in the 2024 observational run

**Effect on the realised headline.** At FaIR's 24.6 % weight, the excess inflated the 2024 realised headline by about 0.045 °C. CGWL's retrospective endpoint was about 0.012 °C *below* its current value, which offsets a little of this.

Two independent checks agree with the result:
- The single HadCRUT5.0 observational run gives a current realised PIVW of 1.360 °C in 2024.
- The 2025 code under the old convention would itself have reported 1.408 °C (git `51de5173`).

**Effect on the anthropogenic headline.** The sign is reversed. The old retrospective FaIR_nonat endpoint was 0.022–0.029 °C *below* its current value. The anthropogenic headline was therefore slightly understated.

### (ii) Change of vintage at fixed year 2024

The 2024-vintage per-member pickles are not in git. Instead, we solved for the uniform level shift that reproduces the saved 2024-vintage P(>1.5 °C) from the 2025-vintage members:
- The Gaussian reconstruction reproduces the saved 2025-vintage probabilities exactly.
- The 2024 ensemble's wider spread (SD ratio 1.10) is applied before solving.

**Pure input effect.** The three temperature-only methods (lowess, GAM, Kal_flexLin) have unchanged code. Their shift is +0.026 to +0.035 °C, which measures the input effect alone. This matches the inputs directly:
- In the new heads/tails ensemble, 2024 itself is 0.036 °C cooler and 2015–2024 is 0.030 °C cooler.
- Warming from 1850–1900 to 1981–2010 is 0.017 °C less.
- The between-member SD is 0.084 vs 0.093, which also lowers exceedance probabilities.

**Effect on the PIVW.** The shift for the realised PIVW itself is 0.023 °C, including the 0.001 °C datum correction (e). This is consistent with the direct estimate in (i) to 0.002 °C: 1.336 + 0.023 = 1.359, against 1.361. It is smaller than the input effect alone (+0.026 to +0.035) because the retuned EBMKF_ta4 now sits slightly *above* its 2024-vintage level (shift −0.006).

**HadCRUT version.** HadCRUT5.0.2 → 5.1 is negligible. It changes the 20-yr means by about −0.005 °C. In the single observational run, on a common 1850–1899 datum, it changes each non-FaIR method's 2024 value by −0.006 to +0.008 °C. The 21 Oct run on 1850–1900 lowers every method in that run by a further 0.0025 °C (CGWL10y_sfUKCP by 0.001–0.002 °C); see (e). The exception is the retuned EBMKF_ta4: +0.015 to +0.035 °C in 2015–2023, a method change rather than a data change. The change of vintage is therefore driven by the new ensemble members, which are drawn from the updated set of datasets.

### (iii) The extra year

Every realised method rose from 2024 to 2025 in **every** ensemble member. The realised PIVW rose by +0.029 °C (minimum +0.021). This happened even though annual temperature fell by about 0.09 °C after the 2024 peak.

Compared with post-1960 non-volcanic years that followed an annual drop of ≥0.10 °C, the non-FaIR methods' 2025 steps fall at the 97th to >99th percentile (EBMKF_ta4: 99.7–99.8). So the realised one-year step matches the expected trend. The realised shortfall of about −0.07 °C against expectation is the sum of (i) and (ii); (e) adds −0.001 °C.

## (c) Method changes: FaIR

FaIR moved beyond the pure input effect.

**Spread-adjusted shifts in the ensemble:**
- FaIR_comb_unB: +0.050 °C (+0.047 before the datum correction)
- FaIR_nonat_unB: +0.067 °C (+0.064)
- For comparison: about +0.03 °C for the temperature-only methods and +0.023 °C for GWI_anthro_CGWL (+0.021)

**Same-year differences in the single observational run:**
- FaIR_nonat: 0.038 °C lower (0.040 including the 1850–1900 datum shift of the 21 Oct run)
- FaIR_comb: 0.013 °C lower (0.016 including the datum shift)
- EBMKF_ta4 (retuned 18–19 Oct): 0.015–0.035 °C *higher* (2015–2023; the 2024 run's EBMKF_ta4 ends in 2023)
- All other methods: within ±0.008 °C

FaIR is also the nearly flat method in 2025. Its steps of +0.011 to +0.013 °C sit at the 30th–35th percentile of its own post-cooling-year steps, and below its median of about 0.02 °C.

FaIR_nonat carries 69 % of the anthropogenic weight. This explains most of the anthropogenic shortfall:
- the larger vintage step, −0.036 vs −0.023 for realised (−0.038 vs −0.025 including the datum correction)
- the one-year step of +0.017, about 0.013 °C short of trend

Shifts for the other methods, for completeness:
- CGWL10y_sfUKCP: +0.056 °C. Its analogue selection is more sensitive to the inputs.
- EBMKF_ta4: −0.006 °C (was +0.017 °C before the 18–19 Oct retuning). Its 2024 level in the new ensemble is 1.332, against 1.323 before, and its 2025 step is +0.039.
- GWI_tot_CGWL: +0.027 °C (+0.024 before the datum correction).

## (b) Weights

The PIVW weights barely changed between vintages: the final weights equal the paper weights to within 0.1 pp. Recomputing the 2025 PIVW level with different weight sets:

| Weights | 2025 PIVW level (°C) |
|---|---|
| Final weights (1/err_var100, `Results2` 21 Oct 02:09) | 1.365 |
| Paper weights (7.5/7.4/7.9/24.6/16.4/25.7/10.5) | 1.365 |
| 1/RMS² from the HadCRUT5.0 (Dec 2025) statistics | 1.364 |
| 1/RMS² from the HadCRUT5.1 statistics (21 Oct, 1850–1900) | 1.364 |
| Equal weights | 1.373 |

The final weights are 7.46/7.44/7.86/24.57/16.46/25.68/10.54 % (realised) and 69.0/31.0 % (anthropogenic). EBMKF_ta4's weight has varied between 10.9 % and 19.0 % across the October recalibrations (14.1 % on 18 Oct). All non-equal weight sets agree to 0.0005 °C.

## (d) How often do methods decrease outside volcanic periods?

The benchmark itself does not decrease:
- The HadCRUT5.1 centred 20-yr mean rises in every year 1966–2015.
- In the heads/tails ensemble, only 0.5 % of member-years since 1966 decrease, by at most 0.003 °C.

**Method.** We computed year-on-year changes Δ(t) of each method's *current* estimate:
- **Volcanic exclusion:** years with AOD > 0.015, plus the 3 following years, are excluded. A +5-year sensitivity gives similar or lower rates. This flags Krakatoa, Santa María, Katmai, Agung, Fuego, El Chichón and Pinatubo; Hunga Tonga is not flagged.
- **Futures:** only years after the switch from historical forcing (MPI ≥2016, NorESM VolcConst ≥2007). The future pickles are the 18 Oct evening runs (20:20–20:52) on the old datum; the 1850–1900 rerun is still on Oscar.

Percentage of years with Δ < 0 (in brackets: Δ < −0.01 °C):

| Method | HadCRUT5.1 ≥1960 | heads/tails ≥1960 | MPI SSP1-2.6 | MPI SSP2-4.5 | MPI SSP3-7.0 | NorESM VolcConst |
|---|---|---|---|---|---|---|
| lowess1dt36wnc | 8 (3) | 8 (2) | 42 (25) | 20 (8) | 6 (2) | 15 (4) |
| GAM_AR1 | 21 (16) | 18 (10) | 43 (31) | 27 (16) | 14 (8) | 23 (12) |
| Kal_flexLin | 16 (3) | 13 (5) | 44 (30) | 26 (16) | 12 (5) | 22 (10) |
| FaIR_comb_unB | 8 (0) | 7 (0.1) | 42 (10) | 6 (1) | 1 (0.2) | 7 (1) |
| EBMKF_ta4 | 8 (0) | 13 (1) | 42 (29) | 25 (14) | 10 (4) | 19 (8) |
| GWI_tot_CGWL | 16 (3) | 9 (1) | 33 (0.1) | 1 (0) | 0.3 (0) | 0.2 (0) |
| CGWL10y_sfUKCP | 26 (21) | 23 (15) | 41 (29) | 28 (19) | 17 (11) | 23 (14) |
| **PIVW (realised)** | 11 (0) | 9 (1) | 41 (20) | 15 (5) | 3 (0.8) | 11 (2) |
| FaIR_nonat_unB | 11 (0) | 7 (0.1) | 42 (13) | 8 (1) | 1 (0.3) | 8 (2) |
| GWI_anthro_CGWL | 11 (3) | 7 (1) | 32 (0.1) | 0.5 (0) | 0.1 (0) | 0.2 (0) |
| **PIVW (anthro)** | 3 (0) | 4 (0) | 40 (9) | 6 (1) | 1 (0.1) | 6 (1) |

**Size of decreases.** Non-volcanic decreases are small:
- Forcing-based methods and PIVW: the mean decrease is −0.002 to −0.008 °C, and the historical 5th-percentile step is above −0.01 °C.
- GAM, Kal_flexLin and CGWL: the mean decrease is −0.014 to −0.023 °C.

**Spurious vs tracking.** In the warming futures (SSP2-4.5, SSP3-7.0, NorESM), 95–100 % of decreases are spurious: they happen while the benchmark is rising. In SSP1-2.6, about 40 % follow genuine stabilisation or cooling.

The full distributions are in `decrease_summary_*.csv`: percentiles, minimum, and results conditional on a preceding ≥0.1 °C annual drop.

**Conclusion for (d).** Such decreases occur, but not in 2025. No realised method fell from 2024 to 2025 in any ensemble member. The only anthropogenic exception is a single FaIR_nonat member at −0.0004 °C. So (d) does not explain the 2024→2025 change.

**Reconciling with the email to Thorne.** The "flat or slightly decreasing" methods noted there are a cross-vintage effect, from the preprint to now:
- FaIR_comb fell (1.566 → 1.365) because of the bug fix.
- FaIR_nonat (1.347 → 1.335) and GWI_tot_CGWL (1.320 → 1.330) were nearly flat because of the input and FaIR changes above.
- EBMKF_ta4 now rises (1.350 → 1.371), after the 18–19 Oct retuning.

## (e) Change of pre-industrial datum, 1850–1899 → 1850–1900

On 21 Oct the datum moved from 1850–1899 to 1850–1900 inclusive (`216d3fde`). The 18 methods anchored to fixed constants (FaIR, GWI, eROF, KCC_human) were then lowered by c = 0.00254 °C (`acf2df22`, `datum_shift.py`). This c is mean(1850–1900) − mean(1850–1899) of HadCRUT5.1. Script: `baseline_hypothesis_1021.py`; tables: `baseline_hypothesis_1021.csv` and `baseline_hypothesis_hist_1021.csv`.

**Inputs.** The heads/tails members were *already* on 1850–1900 in both vintages, because `metaensemblefactory.py` anomalises them to [1850, 1900]. Their 1850–1899 means are −0.0031 °C (2024 vintage) and −0.0029 °C (2025 vintage). Only the single HadCRUT run and the anchored methods' constants were on 1850–1899.

**Ensemble effect, split by state (2025 ensemble means):**

| Change | Effect |
|---|---|
| Seeding of random draws (20 Oct, `0671536e`) | ≤ 1e-5 °C for every method |
| Datum switch in the code, without the correction | 0 for every method; +0.0001 °C on the PIVW from the weight refit |
| Anchored-method correction | −0.00254 °C for FaIR_comb, GWI_tot, FaIR_nonat and GWI_anthro; 0 for the others |
| Net PIVW change | −0.0013 °C (realised, ≈50 % anchored weight) and −0.0025 °C (anthropogenic, 100 % anchored) |

The correction applies equally in 2024 and 2025, so **every 2024→2025 step is unchanged** to 1e-5 °C. In the single HadCRUT run, every method moves by −0.00254 °C in every year, except CGWL10y_sfUKCP (−0.0009 to −0.0022 °C).

**What this means for the question.** The preprint's anchored methods carried the same +0.0025 °C datum offset as the 19 Oct 2025 vintage. The correction is therefore a small level correction applied to the 2025 update only. It adds −0.001 °C (realised) and −0.003 °C (anthropogenic) to the gap against the preprint. It does not affect the one-year step, the exceedance probabilities beyond about 0.2 pp, or the explanation of the shortfall. It does move the rounded headlines: realised 1.37 → **1.36** and anthropogenic 1.34 → **1.33**.

## Caveats

- **The 2024-vintage current levels are inferred.** The values ≈1.361 (realised) and ≈1.355 (anthropogenic) come from the 2024 rate-CSV method endpoints minus the retrospective/current differences measured in the 2025 vintage. The realised value is confirmed independently by the spread-adjusted probability inversion (1.359). The anthropogenic inversion gives 1.370, so treat that step as −0.036 to −0.051 °C (FaIR is sample-based, so the Gaussian inversion is weaker there). A backdated rerun of `93793643` would give exact values if wanted.
- **Confirm the preprint headline source.** Check that the preprint's 1.40 [1.23–1.58] came from this rate CSV (1.406). Its interval is asymmetric, so it may have used a percentile method like the dashboard's.
- **Vintage.** The 2025 historical and heads/tails numbers are from the 21 Oct outputs on the 1850–1900 datum. The (d) future columns are from the 18 Oct runs (old datum) while the 1850–1900 future rerun finishes on Oscar. The datum offset is constant in time, so those step statistics are not expected to change; refresh them once the rerun lands.
- **Rounding.** The realised headline 1.3646 rounds to 1.36, matching the manuscript. The anthropogenic headline 1.3341 now rounds to **1.33**; the 17 Oct manuscript text says 1.34 [1.21-1.50]. Check that the manuscript and this note agree before sending.
