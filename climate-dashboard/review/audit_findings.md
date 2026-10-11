# Current and retrospective endpoint audit

> **Historical audit, before the in-place corrections.** The findings below describe the original code. See [the correction log](../CHANGELOG.md) for applied changes and current-output checks.


The FaIR retrospective formula remains unchanged, as requested. Final dashboard figures are pending this review. All original scripts, saved results, and input files remain unchanged.

## Why PIVW gives two 2025 warming levels

For the 100-member realized ensemble with the preindustrial calibration:

| Method | Current 2025 mean (°C) | Retrospective 2025 mean (°C) | PIVW weight | Contribution to the PIVW difference (°C) |
|---|---:|---:|---:|---:|
| FaIR_comb_unB | 1.367696 | 1.551741 | 24.608% | +0.045290 |
| CGWL10y_sfUKCP | 1.397890 | 1.386392 | 10.509% | −0.001208 |
| Other five methods | Agree, or no retrospective output exists | | 64.883% | 0.000000 |
| PIVW combination | **1.364154** | **1.408236** | 100% | **+0.044082** |

GWI has no retrospective curve, so both versions use its current output. The other four selected methods agree at the 2025 endpoint.

The original `histens_rates.py` selects retrospective histories whenever they exist, then uses their final values for both the headline warming level and the rate calculation. This explains why its saved PIVW warming level is 1.408236°C. The warming level should be read independently from the current outputs if the retrospective history is intended only for fitting rates. The separate working copy `../source/histens_rates_fixed.py` makes that distinction, retaining the retrospective histories for the 2011–2025 rate fits.

## FaIR has two separate discrepancies

The saved FaIR endpoints are reproducible: rerunning the original method for all 100 members exactly reproduced both saved current and retrospective endpoint means. These results are not stale pickles.

**Mixing formula.** In `FaIR_comb_unB_review_copy.py`, line 163 uses the current branch's coefficients 2.45/4, 1/4, and 0.55/4 (0.6125, 0.25, 0.1375), with a four-year lag. They sum to 1. Line 155 uses retrospective coefficients 2/3, 1/3, and 1/6, with a five-year lag. They sum to 7/6. Bias corrections also enter the two expressions differently. Therefore, even identical source endpoint samples do not make these two branches equal. The scientific intent of this difference needs review before altering the formula.

**Confirmed row alignment from the supplied cluster code.** The upstream export uses `ff.define_time(1750, 2027, 1)` and saves temperatures directly on the FaIR `timebounds` axis. Consequently, row 275 is the boundary at the start of 2025, and row 276 is the boundary at the start of 2026. The slicer intentionally uses `year - 1750 + 1`, so its row labeled 2025 is boundary 2026. The retrospective method instead crops from row 100 and reads output index 175 for 2025, selecting boundary 2025.

For all 100 members, and all three source components (all, nonat, and anthro), the finite samples in the current row labeled 2025 exactly equal retrospective row 276, in the same order. All current arrays contain the expected 96 rows for 1930–2025. This rules out a missing-row shift in these supplied combined arrays.

FaIR distinguishes instantaneous temperature boundaries from annual or mid-year averages; averaging the adjacent boundaries would be a different convention. See the official [FaIR time-axis documentation](https://docs.fairmodel.net/en/latest/intro.html#time). The local `core_code_copy_2.2.2/constraining/01_constrain-gsat-rmse-only.py`, lines 124–150, explicitly chooses the boundary following each observation and uses `temp_in[1:n_obs+1]`. Its comments explain that convention. The corresponding 2.1.3 copy makes the same selection. The subsequently supplied Oscar observation-cutoff and RMSE selection blocks match the local 2.2.2 copy exactly. The current slicer's `+1` is therefore consistent with the supplied calibration, and the retrospective extraction is the branch that needs alignment. Changing both branches to annual averages would require reviewing the calibration and baseline definitions as well.

**Effect of aligning rows while retaining the retrospective formula.** A read-only calculation shifts both contemporaneous retrospective inputs and their five-year-lag input forward one row, holding coefficients, biases, and the saved endpoint's random perturbation fixed. It gives:

| FaIR_comb_unB version | Mean across 100 members (°C) |
|---|---:|
| Saved current endpoint | 1.367696 |
| Saved retrospective endpoint | 1.551741 |
| Retrospective endpoint with aligned rows, same formula | 1.579964 |

The index correction adds 0.028223°C to the retrospective mean. Thus it increases the discrepancy from 0.184046°C to 0.212269°C rather than resolving it. The two mixing expressions remain different after alignment. In addition, the current branch's lagged input comes from an earlier-cutoff fit, whereas the retrospective lagged input comes from the full 2025 fit; aligning the endpoint alone does not make these lagged samples identical.

The supplied export also subtracts a simple mean of boundaries 1850–1900, while the local RMSE constraint subtracts a half-weighted-boundary baseline spanning 1850–1901. This baseline difference is shared by the saved current and retrospective source arrays, so it does not explain their relative row mismatch, but should be reconciled before adopting an annual-average convention.

No FaIR coefficient, lag, slice, or cached result has been changed. The retained formula remains an outstanding scientific review item.

**Manual offsets checked in both branches.** All four 2025 `_unB` methods apply their constant bias corrections to both branches. The standalone methods apply equal output shifts. The combined method's net current correction is +0.0026625°C, while its net retrospective correction is −0.0086667°C, because the corrections enter its mixing expressions at different stages. Their difference of −0.0113292°C reduces rather than creates the large retrospective-minus-current discrepancy. The satellite input rebase is also applied to both branches, but its output effect is multiplied by 7/6 in the retrospective combination. Full branch tracing, ensemble verification, and early-period/uncertainty adjustments are documented in `fair_offset_review.md`.

## CGWL also has different endpoint definitions

The current CGWL implementation combines a trailing 10-year observational average with a 10-year forecast average, giving each half the weight. Its retrospective taper uses a 21-year window, combining the available observations and remaining forecast years. At the endpoint this includes 11 observational years and 10 forecast years. These definitions can yield different final means even though they have the same observational endpoint. This smaller difference accounts for −0.001208°C in PIVW.

## The 1.0°C rate inconsistency

The original rate script calculates a 1.0°C endpoint-difference rate and then calculates an OLS slope and slope SE. It updates the SE but fails to assign the OLS slope to `rate10a`. Thus the reported rate and its uncertainty describe different estimators. The working copy adds `rate10a = slope10 * 10.0`, matching the existing 0.5°C OLS convention.

## Seven quartic crossing estimates

The centered 20-year observational benchmark ends in 2015. Seven of the 100 preindustrial members have not crossed 1.0°C by then. The requested quartic procedure uses the last 25 available centered 10-year means, covering 1996–2020. For each member it selects the first real, upward quartic crossing after 2015.

| Zero-based member | Estimated crossing year |
|---:|---:|
| 0 | 2015.817926 |
| 5 | 2016.160830 |
| 29 | 2015.508399 |
| 64 | 2016.624042 |
| 69 | 2017.442898 |
| 72 | 2016.532381 |
| 84 | 2017.331668 |

The combined crossing-year mean is **2011.218093**, using 93 observed crossings and 7 quartic estimates. The empirical 2.5th–97.5th percentile member range is **2004.024272–2016.580504**. These are ensemble ranges, without added uncertainty in the quartic fitting procedure.

For the seven estimated crossing-time rates, the centered 10-year series supplies the available values in the ±5-year OLS window. The quartic supplies only window values after the last available centered 10-year mean (2020). The combined rate is **0.236426°C/decade** with a combined SE of **0.062447°C/decade**. This rate uncertainty combines local OLS fit error and member dispersion, and does not add quartic parameter uncertainty. Both crossing and rate estimates retain explicit estimated-status flags in the member CSV.

The quartic roots lie within the fitted centered 10-year period, although they estimate crossings beyond the last available centered 20-year year. No 1.5°C definitive crossing has been extrapolated.

## Files and reproducibility

- `method_endpoints_2025.csv`: every member's current and retrospective method endpoints.
- `endpoint_difference_contributions.csv`: weighted attribution of the PIVW discrepancy.
- `fair_source_endpoints_2025.csv`: source-array endpoint means and sample counts.
- `fair_rerun_vs_saved.csv`: full 100-member FaIR endpoint rerun comparison.
- `fair_alignment_counterfactual.csv`: per-member effect of shifting the retrospective source rows while holding its formula and endpoint noise fixed.
- `fair_offset_review.md`, `fair_offset_audit_members.csv`, and `fair_offset_audit_summary.csv`: checks of every local `_unB` bias correction and satellite rebase in both branches.
- `quartic_crossing_estimates.csv` and `quartic_crossing_review.png`: the seven quartic estimates.
- `../source/histens_rates_dashboard.py`: independent analysis pipeline, including the centered benchmark, corrected OLS rates, selected-method comparisons, and quartic estimates.
- `../data/rpi_real_member_crossings_rates.csv`: member results with observed/estimated status.

Run the independent pipeline in cleanPy:

```sh
python /absolute/path/to/histens_rates_dashboard.py all --source-root /absolute/path/to/Thorne_15_codefigurestats
```

The copied legacy script reads the original input directory and writes only under the new output directory. It retains the legacy nearest-year crossing convention, so its threshold rows are a regression comparison rather than the final quartic benchmark. The independent dashboard pipeline provides the quartic benchmark.
