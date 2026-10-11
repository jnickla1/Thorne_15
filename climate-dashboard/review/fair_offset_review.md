# FaIR local offset review

> **Historical audit, before the in-place corrections.** The findings below describe the original code. See [the correction log](../CHANGELOG.md) for applied changes and current-output checks.


## Conclusion

The 2025 `_unB` constant bias corrections are applied to both current and retrospective outputs. The three standalone methods apply equal corrections to their two branches. `FaIR_comb_unB` applies corrections to both branches, but their net effects differ because bias subtraction occurs at different stages of the two mixing expressions. This difference reduces, rather than creates, the large retrospective-minus-current warming-level discrepancy.

The newly supplied Oscar observation-cutoff and RMSE boundary-selection blocks match the local 2.2.2 constraint copy exactly. For cutoff year 2025, `n_obs = 176`, the observations are `gmst[0:176]`, and the model selection is `temp_in[1:177]`. With the model array starting at boundary 1850, its final selected boundary is 2026. This confirms that the current slicer's following-boundary convention matches the calibration. The retrospective extraction remains one row earlier.

No original method, source array, cached result, coefficient, lag, or bias has been modified. The retrospective formula remains retained for review.

## Constant bias corrections in the 2025 heads/tails branches

The table gives the change in the output relative to the same method with `curbias` set to zero, keeping all other processing and random perturbations fixed. Positive values increase warming. These corrections are identical for the preindustrial and satellite-calibrated experiments.

| Method | Current bias effect (°C) | Retrospective bias effect (°C) |
|---|---:|---:|
| `FaIR_nonat_unB` | +0.004300 | +0.004300 |
| `FaIR_all_unB` | −0.020000 | −0.020000 |
| `FaIR_anthro_unB` | +0.026200 | +0.026200 |
| `FaIR_comb_unB` | +0.0026625 | −0.0086667 |

`FaIR_nonat_unB` is the FaIR method selected in the dashboard's anthropogenic combination. It is distinct from `FaIR_anthro_unB`, which was also checked here.

### How the combined method handles its corrections

The combined method uses a +0.017°C correction to its non-natural input and a −0.020°C correction to its all-forcing input. Its current expression, for the recent years under review, is:

```text
current = 0.6125 × corrected_nonat_current
        + 0.25   × corrected_all_current
        + 0.1375 × corrected_all_current_lag4

net current bias = 0.6125 × 0.017 − (0.25 + 0.1375) × 0.020
                 = +0.0026625°C
```

The retained retrospective expression is:

```text
retrospective = (2/3) × (raw_nonat_retrospective + 0.017)
              + (1/3) × raw_all_retrospective
              + (1/6) × raw_all_retrospective_lag5
              − 0.020

net retrospective bias = (2/3) × 0.017 − 0.020
                       = −0.0086667°C
```

Thus the all-forcing correction is applied to both current all-forcing components before mixing, while the retrospective expression subtracts 0.020°C once after mixing. It is not missing from the retrospective output, but is applied differently.

The bias contributions to the retrospective-minus-current difference total **−0.0113292°C**. The saved preindustrial difference is **+0.1840457°C**. Without these constant corrections, it would be **+0.1953749°C**. After aligning the retrospective boundary rows while retaining its formula and endpoint noise, the difference is **+0.2122686°C**; removing the biases from that comparison would give **+0.2235978°C**.

These constant shifts do not change the fitted 2011–2025 rate within an individual member. They change levels and potentially threshold crossing times. The mixing coefficients, lags, and row alignment can affect both levels and rates.

## Satellite calibration shift

Every `_unB` method subtracts the same per-member `offset = -mean(temperature[0:50])` from both raw current and retrospective arrays. The mean offset across this ensemble is **0.694354°C**.

For the standalone methods, the final mean shifts by the same amount in both branches. For the combined method, its current coefficients sum to 1, but the retained retrospective coefficients sum to 7/6. Consequently, the common input shift has different output effects:

| Method | Current rebase effect (°C) | Retrospective rebase effect (°C) |
|---|---:|---:|
| Three standalone methods | −0.694354 | −0.694354 |
| `FaIR_comb_unB` | −0.694354 | −0.810080 |

The excess retrospective shift is **−0.115726°C**, exactly one sixth of the mean offset. This is a consequence of the unnormalized retrospective mixing expression. It explains why its current/retrospective discrepancy changes with the calibration baseline, even though the same input offset is applied to both branches.

The code uses the first **50** observations for this local rebase, so this is specifically the implemented 1850–1899 input mean. It is distinct from the 51 annual observations in an inclusive 1850–1900 observational baseline and from the half-weighted boundary baseline in the Oscar constraint. The distinction is recorded here without changing those conventions.

## Other local adjustments

- The uncertainty contraction acts on the current branch only: factors 0.2 for non-natural, 0.4 for all-forcing, and 0.23 for anthropogenic samples. Where its center is the actual ensemble mean, this changes spread without changing the mean, including at 2025. Retrospective spread is not contracted by these factors.
- The current branches replace their provisional centers for 1930–1964 with `0.28 - offset`; the retrospective branches have no corresponding replacement. Because contraction follows the replacement, the final current mean is a weighted blend of this provisional center and the original source mean. This early-period adjustment does not affect the 2011–2025 endpoint and rate windows.
- The time-varying `newcorrecton` alternatives are inactive: the executed branch is guarded by `if (True)`. Those alternatives are not contributing to the saved historical-ensemble results.
- Both branches receive small random perturbations after these steps. The audit holds those perturbations fixed when isolating a correction's effect.

## Verification and saved tables

The review reran all four `_unB` methods for 100 members under both calibration options, then repeated the calculations with the constant biases removed in memory and, for the satellite option, with its rebase removed in memory. The original source files were read only. For all methods present in the saved result pickles, both rerun endpoint means match the saved endpoints to the numerical precision checked (maximum difference below 1e-12°C; observed maximum zero).

- `fair_offset_audit_members.csv`: 800 method/member/calibration rows with isolated bias and rebase effects.
- `fair_offset_audit_summary.csv`: ensemble means of those quantities.
- `fair_alignment_counterfactual.csv`: the separately checked retrospective row-alignment effect.

The remaining major issue is the difference between the current and retrospective mixing definitions, particularly the retrospective coefficient sum of 7/6. No change to that formula has been adopted.
