# 2024 → 2025 headline investigation (read-only analysis)

This folder answers Peter Thorne's question of why the 2025 update's numbers are not higher than 2024's. Nothing outside this folder was modified.

## Documents

- `supplement_note_2024v2025.md`: the findings, a decomposition table, the (d) decrease statistics, and a one-paragraph reply.
- `typos_paper_and_supplement.md`: typos found in the main paper and supplement drafts.

## Inputs

- **2024 update (HadCRUT5.0, data to 2024):** git `93793643`. These were extracted into `snap24/` (heads/tails statistics and probability arrays) and `hist_results_93793643.pickle`.
- **2025 update (paper state):**
  - git `f36ba84b`, extracted into `snap25/` and `hist_results_f36ba84b.pickle`
  - heads/tails pickles and `Results2` from the 10 Oct 2026 backup, `.../outputs/kalman-rerun-2026-10-13/before/`. This backup is byte-identical to the `f36ba84b` heads/tails CSVs.
  - future pickles from `f36ba84b`, extracted to the session scratchpad.

## Scripts

Each script writes a table:

| Script | Covers | Output |
|---|---|---|
| `compare_inputs.py` | (a) HadCRUT5.0.2 vs 5.1, and the 2024 vs 2025 heads/tails ensembles | printed |
| `compare_probs.py` | per-method P(>1.5 °C) by vintage | `p15_by_method_vintage.csv` |
| `compare_hist_run.py` | same-year comparison of the single observational run | `hadcrut_run_2024v2025.csv` |
| `equiv_shift.py` | 2024-vintage level shift implied by the saved probabilities | `equivalent_level_shift.csv` |
| `equiv_shift_spread.py` | spread-adjusted level shift, plus the direct estimate (used in the note) | `equivalent_level_shift_spread_adjusted.csv` |
| `compare_weights.py` | (b) weights across saved calibrations | `weights_by_vintage.csv` |
| `step2025_and_checks.py` | 2025 steps, cent20y monotonicity, weight sensitivity | `step_2024to2025_headstails.csv` |
| `decrease_analysis.py` | (d) non-volcanic year-on-year decreases | `decrease_summary_*.csv`, `steps_nonvolc_*.csv.gz` |

`invert_probs.py` was a failed approach, kept for the record: P saturates at the 1.0 °C threshold.

## 18 Oct 2026 refresh

The 2025 vintage was swapped for the outputs regenerated on 18 Oct 2026: the retuned EBMKF_ta4, the heads/tails ensemble, `hist_results.pickle`, the rates/dashboard, the futures and the `Results2` weights. The 17 Oct files above are kept unchanged for comparison.

- `changes_1017_vs_1018.md` (+ `.csv`): every number in the note, 17 Oct vs 18 Oct, with a flag for whether the rounded text changes. `build_changes_1018.py` writes the whole file (header from `changes_header_1018.md`).
- `supplement_note_2024v2025_1018.md`: the refreshed note and Thorne reply.
- `cfg_1018.py`: shared paths, and weights recomputed as 1/`err_var100` from `Results2/historical_names_var_scale.csv`.
- `freeze_inputs_1018.py` → `inputs_1018_manifest.csv`: path, mtime and md5 of every input. `snap25_1018/` and `hist_results_1018.pickle` are frozen copies.
- `diff_old_new_1018.py` → `diff_old_new_1018.txt`; `step2025_and_checks_1018.log` holds the cent20y check: which methods changed. Only EBMKF_ta4 changed materially.
- `*_1018.py` → `*_1018.csv`: 18 Oct versions of `compare_probs`, `compare_hist_run`, `equiv_shift`, `equiv_shift_spread`, `compare_weights`, `step2025_and_checks` (plus `weight_sensitivity_1018.csv` and `step_2024to2025_members_1018.csv`) and `decrease_analysis` (log in `decrease_analysis_1018.log`).
- `step_context_post1960.py`: rebuilds `step2025_in_context*.csv` and `decrease_summary_hist_post1960_n*.csv`. These were produced inline on 17 Oct. Run with `''` to verify against the 17 Oct files (exact match), or with `_1018` for the new vintage.

Run order: `freeze_inputs_1018`, `compare_probs_1018`, `compare_hist_run_1018`, `equiv_shift_1018`, `equiv_shift_spread_1018`, `compare_weights_1018`, `step2025_and_checks_1018`, `decrease_analysis_1018`, `step_context_post1960 _1018`, `build_changes_1018`. Run them all from this folder.

## 19 Oct 2026 final refresh

This refresh uses the outputs regenerated on 18–19 Oct. Table S4 (CSCM) and Figure 17 were still being worked on, but neither feeds these numbers.

- `supplement_note_2024v2025_1019.md`: the final note and Thorne reply.
- `changes_1018_vs_1019.md` / `.csv`: 18 Oct vs 19 Oct, written by `build_changes_1019.py` (header in `changes_header_1019.md`).
- `cfg_1019.py` reads the weights and comb scales live from `../Results2/` (the most current files). `results2_1019/` keeps a copy of the 19 Oct versions, and the config prints a NOTE if the live files differ from it. The other inputs are frozen in `snap25_1019/` and `hist_results_1019.pickle`; the manifest is `inputs_1019_manifest.csv`.
- `cfg_1018.py` now hard-codes the 18 Oct weights, because `Results2` was overwritten on 18 Oct at 22:37.
- `diff_1018_vs_1019.py` → `diff_1018_vs_1019.txt`. Only EBMKF_ta4 changed (+0.0015 °C), and the weights returned to the paper weights.
- Pipeline: the same scripts with the `_1019` suffix. The log is in `pipeline_1019.log`; to rebuild the context tables, run `step_context_post1960.py _1019`.
