# fut_LOO: future leave-one-out combinations (Supplement Table S4)

Run every command from the **repository root**. The scripts import `custom_stack_utils` and read `Results2/` relative to it.

## What's here

| path | purpose |
|---|---|
| `create_combination_method_LLOffcompute.py` | The leave-one-out kernel. For one selection mode, scenario and held-out member, it fits PIVW and CSCM on the other members and scores them on the held-out one. |
| `make_oscar_worker.py` | Regenerates `LLOffcompute_fix_worker.py` from the kernel, and `fut2_fix_worker.py` from `fut_2_script.py`. |
| `fut_2_script.py` | The individual reference methods (mode 2: EBMKF_ta4, GWI_tot_CGWL), scored over 2000–2089. |
| `fut2_fix_worker.py` | Generated; don't edit. One member per call, writes `Results3_fix/nmes2run_*.npz`. |
| `LLOffcompute_fix_worker.py` | Generated file; don't edit it. It is the kernel with routing-only changes: 261018 selection statistics, `Results3_fix/` output, and saved `n_valid_kl`. |
| `run_LOO_fix.sh` | Oscar array: 1,350 tasks (modes 7, 16, 32, 8 × 270 members, then mode 2 × 270), 24 GB and 6 h each. Finished cases are skipped. |
| `check_and_aggregate.py` | Coverage and finite checks, optional cross-check against the 2026-10-19 Mac run, Table S4 and the comparison with Dec 2025. |
| `selection.json` | The actual method lists for each mode: 7, 8, 18 and 33 methods for modes 7, 8, 16 and 32. |
| `newsamples_patch/` | The one-time 2026-10-19 fix of `Results2/*_newsamples.npy`, already applied on the Mac and on Oscar, with the input-hash lists used to confirm both sides matched. |
| `results/Results3_fix_oscar/` | The 1,080 per-case results from Oscar array 7637718. |
| `results/mac_fallback/` | The Mac computation of mode 32 / NorESM Volc / member 59, used as a cross-check. |
| `output/` | Table S4 (`table-S4.md`, `-numeric.csv`, `-formatted.csv`, `-for-docs.xlsx/.csv`), `final_combined_*`, both KL normalizations and `comparison_old_vs_new.csv`. |
| `legacy/` | Pre-October entry points (`run_LOO*.sh`, `create_comb_method_LLOffcompute_percase.py`, `create_combination_method_LLOffcombine.py`, `fut_2_script.py`, the mode-2 individual-method comparison) and the Dec 2025 `final_combined_*` outputs. |
| `OSCAR_STEPS_2026-10-19.md` | The step-by-step log of the October 2026 Oscar run. Its file paths predate this folder. |

## Regenerating Table S4

1. **Inputs.** `Results2/{scenario}{member}_*` must hold the 270 future calibration sets (MPI SSP126/245/370 × 50, NorESM Volc/VolcConst × 60). They need `newsamples` produced by the fixed `calibrate_method_uncertainty.py`, which uses `*best_alter_scale`.
   - Sets produced before 2026-10-19 need `newsamples_patch/patch_newsamples_scale.py`, applied once only; the manifest prevents double application.
   - The `averaged_runs*.csv` files and `261018_HadCRUT5.1_var_unadjusted/current_methods_statistics_261018True.csv` choose the methods for each mode.
2. **If the kernel changed,** run `python fut_LOO/make_oscar_worker.py`.
3. **On Oscar:** `mkdir -p Results3_fix/logs && sbatch fut_LOO/run_LOO_fix.sh`. Most cases take under 1 h. Some 33-method cases need 2–4 h, because the SLSQP stacking fit runs to its 100-iteration limit; the weights settle well before that.
4. **Back on the Mac:** `rsync -av --include='*.npz' --exclude='*' <oscar>:<repo>/Results3_fix/ fut_LOO/results/Results3_fix_oscar/`, then `python fut_LOO/check_and_aggregate.py fut_LOO/results/Results3_fix_oscar fut_LOO/output`.

## Pasting Table S4 into Google Docs

1. Open `output/table-S4-for-docs.xlsx` in Excel or Google Sheets.
2. Select the filled range and copy it.
3. Paste into the Google Doc. If you paste from Sheets, choose "Paste unlinked"; it comes in as a native table.

The CSV version is UTF-8 and opens in Excel.

## Method notes

- **PIVW weights** use inverse error variances over 2000–2089 for every ensemble. Before 2026-10-19, MPI used 2000–2090.
- **PIVW predictive SD:** a single constant, chosen by KL over 2000–2089.
- **Held-out scoring** (KL and RMSE, full and from 2025) uses 2000–2089 for every ensemble; MPI's 251-year benchmark also has 2090, which is excluded since 2026-10-20. Crossing metrics use the full series.
- **CSCM:** stacking weights from Gauss–Hermite log-likelihoods, then a compression factor `nsharp` from a bounded scalar search on [1, 10].
- **KL** in the main table uses the legacy 90/65-year normalization. `combination-summary-per-valid-year.csv` normalizes by the actual number of finite benchmark years instead.
- **Reproducibility:** Mac and Oscar results agree to floating-point precision (CSCM KL within 0.34 % per case). The per-case seeds are fixed in the kernel.
- **Provenance:** the seven 33-method MPI cases listed in `results/mac_reruns_2089eval/` were computed on the Mac, because their Oscar tasks hit the 6 h limit. They are identical to Oscar runs to floating-point precision.
- **Individual reference blocks:** `check_and_aggregate.py` appends EBMKF_ta4 and GWI_tot_CGWL blocks to `table-S4-for-docs` when mode-2 results are found (optional third argument). The current blocks come from `results/mode2_individual_2026-10-18/`, which scored MPI over 2000–2090; a note row in the xlsx says so until mode 2 is rerun.
