# Why the current-warming estimate fell slightly (realised) or stayed flat (anthropogenic) from the 2024 to the 2025 update

*Draft supplement section, **19 Oct 2026 (final) refresh**. Numbers come from:*
- *2024 update: the HadCRUT5.0-based state at git `93793643` (20 Dec 2025, data through 2024). Unchanged.*
- *2025 update: the current outputs, all regenerated 18–19 Oct 2026. These are the heads/tails ensemble (18 Oct 18:47–22:58), the futures (18 Oct 20:20–20:52), `hist_results.pickle` and the var_adjusted statistics (19 Oct 01:23–01:32), the weights (`Results2` var_scale 18 Oct 22:37; comb scales 19 Oct 13:49) and the dashboard (19 Oct 16:19). Files, mtimes and md5s are in `inputs_1019_manifest.csv`.*
- *Relative to the 17 Oct note, only EBMKF_ta4 changed materially (retuned), and the weights are back to the paper weights to within 0.1 pp. See `changes_1017_vs_1018.md` and `changes_1018_vs_1019.md`.*
- *Table S4 (CSCM crossing comparisons) and Figure 17 (Kalman-filter attributable warming) are still being finalised. Neither feeds the PIVW numbers here: CSCM is a separate combination, and Figure 17 uses the attribution output, not EBMKF_ta4's current estimate. If the Figure 17 work changes EBMKF_ta4 itself, rerun the `*_1019.py` scripts.*

*Scripts and output tables are in `thorne_2024v2025_investigation/` (`*_1019.py`, `*_1019.csv`).*

## Headline

The 2024 values are from the ESSD preprint abstract (essd-2025-825). The 2025 values are from the 19 Oct dashboard export. The realised value, 1.3659, rounds to 1.37; the manuscript may still say 1.36.

| | 2024 (preprint) | 2025 (19 Oct outputs) | Change | Expected at ~0.3 °C/decade |
|---|---|---|---|---|
| Realised | 1.40 [1.23–1.58] °C, "unlikely" >1.5 °C | 1.37 [1.23–1.53] °C, "very unlikely" | −0.03 | ≈ +0.03 |
| Anthropogenic | 1.34 [1.18–1.50] °C, "very unlikely" | 1.34 [1.21–1.50] °C, "extremely unlikely" | 0.00 | ≈ +0.03 |

**Where the numbers come from.**
- **2024 values:** the PIVW values in the 2024 rate CSV (`rates_histens/`, git `876b4b8c`/`93793643`), 1.406 and 1.340 °C. The script took each method's last *retrospective* value (`LT_trend[2][-1]`) and used the current output only when no retrospective value existed.
- **2025 values:** the dashboard export of *current* outputs, 1.366 [1.230, 1.532] and 1.337 [1.209, 1.496] °C.
- **Exceedance probabilities:** the mean member P(>1.5 °C) is always computed from current outputs. It fell from 7.5 % to 6.4 % (realised) and from 7.1 % to 2.1 % (anthropogenic).

**Explanations tested.**
- (a) Input data
- (b) Method weights
- (c) Internal method changes
- (d) Methods decreasing while tracking the 20-yr mean

The shortfall against expectation comes from a reporting correction, (a), and a smaller FaIR contribution under (c). The retuned EBMKF_ta4 is a small (c) effect in the opposite direction. Neither (b) nor (d) contributes materially.

## Decomposition (PIVW, heads/tails ensemble, 1850–1900 baseline)

| Step | Realised level (°C) | Change | Anthropogenic level (°C) | Change |
|---|---|---|---|---|
| 2024 update as reported (retrospective endpoints) | 1.406 | | 1.340 | |
| (i) Same run, current endpoints | ≈1.361 | ≈ −0.045 | ≈1.355 | ≈ +0.015 |
| (ii) 2025 inputs and method versions, still year 2024 | 1.337 | ≈ −0.023 | 1.319 | ≈ −0.036 |
| (iii) One more year, 2024 → 2025 | 1.366 | +0.029 | 1.337 | +0.017 |

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

**Effect on the PIVW.** The shift for the realised PIVW itself is 0.022 °C. This is consistent with the direct estimate in (i) to 0.002 °C: 1.337 + 0.022 = 1.359, against 1.361. It is smaller than the input effect alone (+0.026 to +0.035) because the retuned EBMKF_ta4 now sits slightly *above* its 2024-vintage level (shift −0.006).

**HadCRUT version.** HadCRUT5.0.2 → 5.1 is negligible. It changes the 20-yr means by about −0.005 °C. In the single observational run, it changes each non-FaIR method's 2024 value by −0.006 to +0.008 °C. The exception is the retuned EBMKF_ta4: +0.015 to +0.035 °C in 2015–2023, a method change rather than a data change. The change of vintage is therefore driven by the new ensemble members, which are drawn from the updated set of datasets.

### (iii) The extra year

Every realised method rose from 2024 to 2025 in **every** ensemble member. The realised PIVW rose by +0.029 °C (minimum +0.021). This happened even though annual temperature fell by about 0.09 °C after the 2024 peak.

Compared with post-1960 non-volcanic years that followed an annual drop of ≥0.10 °C, the non-FaIR methods' 2025 steps fall at the 97th to >99th percentile (EBMKF_ta4: 99.7–99.8). So the realised one-year step matches the expected trend. The realised shortfall of about −0.07 °C against expectation is the sum of (i) and (ii).

## (c) Method changes: FaIR

FaIR moved beyond the pure input effect.

**Spread-adjusted shifts in the ensemble:**
- FaIR_comb_unB: +0.047 °C
- FaIR_nonat_unB: +0.064 °C
- For comparison: about +0.03 °C for the temperature-only methods and +0.021 °C for GWI_anthro_CGWL

**Same-year differences in the single observational run:**
- FaIR_nonat: 0.038 °C lower
- FaIR_comb: 0.013 °C lower
- EBMKF_ta4 (retuned 18–19 Oct): 0.015–0.035 °C *higher* (2015–2023; the 2024 run's EBMKF_ta4 ends in 2023)
- All other methods: within ±0.008 °C

FaIR is also the nearly flat method in 2025. Its steps of +0.011 to +0.013 °C sit at the 30th–35th percentile of its own post-cooling-year steps, and below its median of about 0.02 °C.

FaIR_nonat carries 69 % of the anthropogenic weight. This explains most of the anthropogenic shortfall:
- the larger vintage step, −0.036 vs −0.023 for realised
- the one-year step of +0.017, about 0.013 °C short of trend

Shifts for the other methods, for completeness:
- CGWL10y_sfUKCP: +0.056 °C. Its analogue selection is more sensitive to the inputs.
- EBMKF_ta4: −0.006 °C (was +0.017 °C before the 18–19 Oct retuning). Its 2024 level in the new ensemble is 1.332, against 1.323 before, and its 2025 step is +0.039.
- GWI_tot_CGWL: +0.024 °C.

## (b) Weights

The PIVW weights barely changed between vintages: the final weights equal the paper weights to within 0.1 pp. Recomputing the 2025 PIVW level with different weight sets:

| Weights | 2025 PIVW level (°C) |
|---|---|
| Final weights (1/err_var100, `Results2` var_scale 18 Oct 22:37) | 1.366 |
| Paper weights (7.5/7.4/7.9/24.6/16.4/25.7/10.5) | 1.366 |
| 1/RMS² from the HadCRUT5.0 (Dec 2025) statistics | 1.366 |
| 1/RMS² from the HadCRUT5.1 statistics (19 Oct var_adjusted) | 1.366 |
| Equal weights | 1.374 |

The final weights are 7.46/7.44/7.86/24.58/16.47/25.69/10.50 % (realised) and 69.0/31.0 % (anthropogenic). EBMKF_ta4's weight has varied between 10.9 % and 19.0 % across the October recalibrations (14.1 % on 18 Oct). All non-equal weight sets agree to 0.0005 °C.

## (d) How often do methods decrease outside volcanic periods?

The benchmark itself does not decrease:
- The HadCRUT5.1 centred 20-yr mean rises in every year 1966–2015.
- In the heads/tails ensemble, only 0.5 % of member-years since 1966 decrease, by at most 0.003 °C.

**Method.** We computed year-on-year changes Δ(t) of each method's *current* estimate:
- **Volcanic exclusion:** years with AOD > 0.015, plus the 3 following years, are excluded. A +5-year sensitivity gives similar or lower rates. This flags Krakatoa, Santa María, Katmai, Agung, Fuego, El Chichón and Pinatubo; Hunga Tonga is not flagged.
- **Futures:** only years after the switch from historical forcing (MPI ≥2016, NorESM VolcConst ≥2007). The future pickles are the 18 Oct evening runs (20:20–20:52), including the SSP245 rerun.

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
| **PIVW (realised)** | 13 (0) | 9 (1) | 41 (20) | 15 (5) | 3 (0.8) | 11 (2) |
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
- FaIR_comb fell (1.566 → 1.368) because of the bug fix.
- FaIR_nonat (1.347 → 1.338) and GWI_tot_CGWL (1.320 → 1.332) were nearly flat because of the input and FaIR changes above.
- EBMKF_ta4 now rises (1.350 → 1.371), after the 18–19 Oct retuning.

## Suggested reply to Peter Thorne (one paragraph)

> Thanks, Peter. You were right that the new datasets pulled things down. This year's heads/tails ensemble, drawn from the expanded set of datasets, is about 0.03 °C cooler than last year's over the last decade, with a slightly narrower spread. On its own, that lowers the methods' end-of-2024 estimates by about 0.03 °C, roughly one year's warming. The preprint's realised 1.40 °C was also inflated by about 0.04–0.05 °C. That headline read each method's retrospective endpoint, and the retrospective branch of one FaIR combination had a weighting error, which we have since corrected. The 1.5 °C exceedance probabilities and the current estimates were unaffected. Like for like, realised warming at the end of 2024 was about 1.36 °C on last year's data and 1.34 °C on this year's. It then rose by 0.03 °C, in line with the trend, to 1.37 °C (1.366) in 2025. Every method increased in every ensemble member, even though 2025 was about 0.1 °C cooler than 2024. The anthropogenic estimate is flat at 1.34 °C, mainly because the updated FaIR runs put human-induced warming about 0.03–0.04 °C lower than before and rose only about 0.01 °C in 2025. The changes in method weighting and from HadCRUT5.0 to 5.1 have negligible effect (≤0.005 °C). Year-on-year decreases in the methods do occur outside volcanic periods, in about 4–20 % of years and typically by less than 0.01 °C, but none occurred in 2025.

## Caveats

- **The 2024-vintage current levels are inferred.** The values ≈1.361 (realised) and ≈1.355 (anthropogenic) come from the 2024 rate-CSV method endpoints minus the retrospective/current differences measured in the 2025 vintage. The realised value is confirmed independently by the spread-adjusted probability inversion (1.359). The anthropogenic inversion gives 1.370, so treat that step as −0.036 to −0.051 °C (FaIR is sample-based, so the Gaussian inversion is weaker there). A backdated rerun of `93793643` would give exact values if wanted.
- **Confirm the preprint headline source.** Check that the preprint's 1.40 [1.23–1.58] came from this rate CSV (1.406). Its interval is asymmetric, so it may have used a percentile method like the dashboard's.
- **Vintage.** All 2025 numbers are from the 18–19 Oct outputs, including the retuned EBMKF_ta4 (≈16 % weight) and the 18 Oct evening future ensembles. Table S4 (CSCM) and Figure 17 were still in progress at the time but do not feed these numbers.
- **Rounding.** The realised headline 1.3659 rounds to 1.37, not 1.36. Check that the manuscript and this note agree before sending.
