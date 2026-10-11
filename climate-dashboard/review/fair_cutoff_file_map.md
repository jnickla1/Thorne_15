# Where the FaIR cutoff inputs are stored

`cut2021` names an Oscar trajectory constrained using observations through 2021. It has now been transferred as a separate full trajectory in the local `retrospective_4/` folder. The original reslicer selected one row from each cutoff trajectory and stacked those rows into the local `combined_...` array. `cut2025` full trajectories were also transferred into the local `retrospective/` folder.

The original trajectories are created before reslicing. The distinction is which trajectory supplied a value: the 2021 cutoff fit or the 2025 cutoff fit. The current mixing expression combines values from different cutoff fits. The latest user instruction keeps that expression and makes the retrospective lag use the 2021 cutoff trajectory.

## Local files used by the combined method

Reader:

`/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats/Methods/43_Forcing_Based/1_ERF_FaIR/FaIR_comb_unB_method.py`

For the first heads/tails member (`model_run=0`, which maps to `r1`):

| Input | Absolute local filename | Relevant zero-based row |
|---|---|---|
| Current all-forcing history | `/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats/Methods/43_Forcing_Based/1_ERF_FaIR/resliced_headstails25/combined_all-headstails25r1_all.npy` | Row 91 is labeled 2021; row 95 is labeled 2025 |
| Contemporary retrospective all-forcing trajectory | `/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats/Methods/43_Forcing_Based/1_ERF_FaIR/resliced_headstails25/retrospective/all-headstails25r1cut2025_temp_all.npy` | Row 276 supplies the 2025-labeled value |
| Earlier-cutoff retrospective lag trajectory, transfer verified | `/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats/Methods/43_Forcing_Based/1_ERF_FaIR/resliced_headstails25/retrospective_4/all-headstails25r1cut2021_temp_all.npy` | Row 272 supplies the 2021-labeled lag for 2025 |

The current array has shape `(96, 841)` and labels 1930–2025. The retrospective array has shape `(278, 841)` and temperature boundaries 1750–2027. The reslicer uses the following boundary: label 2021 selects boundary 2022, and label 2025 selects boundary 2026. Comments in the existing readers that describe the current array as starting in 1750, or the retrospective array as starting in 1930, are reversed; the inspected array shapes and supplied export code establish the mapping above.

The historical central-observation branch uses the same filenames with `r0`. Heads/tails members use `r1` through `r100`.

## Original Oscar filenames

From the supplied reslicer, the first-member source files are:

```
/users/jnickla1/data/jnickla1/fair-calibrate/output/all-headstails25/all-headstails25r1cut2021_temp_all.npy
/users/jnickla1/data/jnickla1/fair-calibrate/output/all-headstails25/all-headstails25r1cut2025_temp_all.npy
```

Before the latest correction, the current 2025 mixture took its lagged input from local combined row 91, originally extracted from the first file's row 272. The retrospective 2025 mixture took its lagged input from the second file's row 272. These are the same temperature boundary but different observation cutoffs. Both can be generated before the reslicer is run and still differ.

For example, the inspected r1 raw all-forcing sample means are 1.250629675°C in combined row 91 and 1.222834725°C in cut2025 row 272, before local bias, uncertainty, calibration, and mixing corrections. The difference exists in the stored inputs.

## Selected correction: change the retrospective lag only

The user's latest instruction supersedes the proposal to change the current estimator or export matched rows from every cutoff. The 2025 retrospective mixture reads a second complete all-forcing matrix constrained through 2021 and use its slice `[97:273]` as the four-year-lagged history. Contemporary non-natural and all-forcing inputs still come from the cut2025 fits. The coefficients and biases remain the corrected values, and full retrospective spread is retained.

The satellite-calibration offset is applied once to the additional matrix. Its columns are padded or cropped to the final retrospective matrix's width, retaining source column order. This matches the available prefix used by the current mixture and preserves the number of subsequent random draws. Older and future input branches retain their existing retrospective lag selection.

For the final 2025 estimate, the current and retrospective lag now both originate from the cut2021 fit's 2021-labeled value. The earlier values in the new retrospective lag remain retrospective estimates from that 2021 fit; they should not be described as estimates available in each earlier year. Current spread contraction and separate random perturbations remain as defined, so exact equality of current and retrospective sample distributions is not being imposed.

## One transfer command

Run this on the Mac. It copies only the full **all-forcing cut2021** trajectories for runs r0–r100, using their existing names:

```sh
scp 'jnickla1@ssh.ccv.brown.edu:/users/jnickla1/data/jnickla1/fair-calibrate/output/all-headstails25/all-headstails25r*cut2021_temp_all.npy' '/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats/Methods/43_Forcing_Based/1_ERF_FaIR/resliced_headstails25/retrospective_4/'
```

No nonat or anthro cut2021 matrices are needed. No new FaIR simulations or reslicing are needed if the original Oscar files remain available. No additional method modules are introduced.

The prior `export_fair_same_cut_lag4.py` exporter and its NPZ are not needed for this selected approach. They belong to the superseded current-branch proposal.

Transfer and validation are complete. All 101 raw trajectories' row 272 values match the corresponding combined arrays' 2021 rows exactly. The 400 corrected retrospective cache copies, all four dashboard datasets, and the figure library have been refreshed. Current samples and derived current curves remain unchanged. Existing current caches need not be regenerated for this correction. See [the change log](../CHANGELOG.md) for validation and data status.
