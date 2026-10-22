# FaIR runs on Oscar: pipeline pieces outside the constraint-set scripts

The constrained FaIR results used in the paper (`FaIR_*_unB_method.py`, `FaIR_comb_unB_method.py`) and the TOA ensemble read by `2_Kalman/EBMKF_ta4_method.py` were produced on Brown's Oscar cluster. The code ran in a checkout of John Nicklas's fork of fair-calibrate:

- Fork: https://github.com/jnickla1/fair-calibrate, commit `7c96064` on top of `2c6ad33`. The fork was cloned from upstream https://github.com/chrisroadmap/fair-calibrate at tag `v1.4.5` (`120ccfa`).
- Constraint set: `input/fair-2.2.2/v1.4.3/all-headstails25/`, run with FaIR 2.2.2. The folder is named v1.4.3, but `04_dump-calibration.py` reads the v1.4.1 calibrated parameters and species configs from `data2/calibration/v1.4.1/`.

## Where each piece is in this repository

| Oscar location (relative to the fair-calibrate root) | Here |
|---|---|
| `input/fair-2.2.2/v1.4.3/all-headstails25/{calibration,sampling,constraining}/` | `../core_code_copy_2.2.2/`. All 38 numbered and TOA scripts are byte-identical to Oscar, checked by md5 on 2026-10-22. |
| `.env222` (environment for the run) | `env/.env222` |
| `run_mult_cutoffs_sampALLsteps2` (SLURM: calibration + sampling) | `drivers/run_mult_cutoffs_sampALLsteps2` |
| `run_mult_const_rep_headstails2` (SLURM array: constraining per run × cutoff) | `drivers/run_mult_const_rep_headstails2` |
| `data/forcing/sst_pseudo_25_hadcrut51.csv` (observational ensemble used for the GMST constraint; column 1 = central HadCRUT5.1, i.e. run r0; columns 2–101 = heads/tails members r1–r100) | `data/sst_pseudo_25_hadcrut51.csv` |
| `reslice/reslice_current_HT3.py` | `reslice/reslice_current_HT3.py` |
| `reslice/reslice_current_HT3_lag4.py` | `reslice/reslice_current_HT3_lag4.py` |

## Workflow

The drivers are submitted from the fair-calibrate root with `sbatch`. They read `.env222` from there, use the `fair2` conda environment (FaIR 2.2.2), and need the R `EBM` package for calibration step 02.

1. **Calibration and prior ensemble.** `drivers/run_mult_cutoffs_sampALLsteps2` runs the numbered scripts in `calibration/` and `sampling/` in order. No cutoff is set at this stage.
   - As archived, its loop runs `sampling` only; the `calibration` entry is commented out because calibration had already been run with the same script.
2. **Constraining.** `drivers/run_mult_const_rep_headstails2` is a SLURM array over `ENSEMBLE_RUN` 0–100. Each task uses GNU `parallel` across `CUTOFF_YEAR` 1930–2025.
   - For each cutoff it first runs `constraining/01_constrain-gsat-rmse-only.py`, `02_run-1pct.py` and `03_reweight-rmse-posterior-multiple-constraints.py`, and then `04_dump-calibration.py`.
   - `04` writes `output/all-headstails25/all-headstails25r{run}cut{year}_temp_{all,nonat,anthro}.npy`. Each array has shape (278, 841), on timebounds 1750–2027.
3. **Reslicing.** `reslice/reslice_current_HT3.py` takes, from each cutoff-year trajectory, the row for the following timebound (`year - 1750 + 1`). It stacks them into `combined_all-headstails25r{run}_{all,nonat,anthro}.npy`, with shape (96, 841) and annual labels 1930–2025.
   - The results are the files in `../resliced_headstails25/`.
   - The full `cut2025` and `cut2021` trajectories were copied unchanged into `../resliced_headstails25/retrospective/` and `retrospective_4/`.
   - The script only accepts trajectories modified in Sep–Nov 2026, so stale outputs are skipped.
4. **TOA ensemble for EBMKF_ta4.**
   - `constraining/16_TOA-forcing_25.py` draws one constrained parameter set per member (seed 726371 + member). It drives FaIR with total ERF and writes per-member TOA imbalance and surface temperature for 1850–2025.
   - `17_concat_toa_25.py` stacks those into `all-headstails_current_TOA_ensemble_25.npy`, which is read from `2_Kalman/heads_tails_forcing/`.

`reslice_current_HT3_lag4.py` (Oct 2026) exports matched current/four-year-lag rows from each cutoff. It is archived because it was written during the cut2021 review, but the paper results do not use its output. The adopted correction instead reads the full `cut2021` trajectories in `retrospective_4/` (see `climate-dashboard/review/fair_cutoff_file_map.md`).

Run outputs are not archived here; they total several GB on Oscar. The exported arrays that the method readers consume are in the sibling `resliced_*` and `retrospective*` folders.
