from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt

#only does realized warming for the moment

from hist_heads_tails_evaluation_script import evalmins, evalmaxs

n_methods, n_fine, n_thr = 5, 36*12, 3

from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

import pandas as pd

THRESH_COLS = ["xyear0.5", "xyear1.0", "xyear1.5"]

def discover_histens_files(
    experiment_type: str,
    results_dir: str | Path = "Results"
) -> Tuple[List[Path], List[Path]]:
    """
    Return lists of CSV and NPY result files for a given historical ensemble experiment.
    Prints counts for visibility.
    """
    results_dir = Path(results_dir)
    csv_files = sorted(results_dir.glob(f"headstails_statistics_{experiment_type}*.csv"))
    npy_files = sorted(results_dir.glob(f"headstails_method_fineprobs_{experiment_type}*.npy"))
    print(f"[{experiment_type}] found {len(csv_files)} CSV(s) and {len(npy_files)} NPY(s) in {results_dir}")
    return csv_files, npy_files

def read_method_names_from_csvs(csv_files: List[Path]) -> List[str]:
    """
    Read method names once from the first CSV (first col named 'method_name'; 
    falls back to the first column if not present). Warn if later CSVs differ.
    """
    if not csv_files:
        return []
    df0 = pd.read_csv(csv_files[0])
    col0 = "method_name" if "method_name" in df0.columns else df0.columns[0]
    names0 = df0[col0].astype(str).tolist()

    mismatches = 0
    for p in csv_files[1:]:
        df = pd.read_csv(p)
        coli = "method_name" if "method_name" in df.columns else df.columns[0]
        if list(df[coli].astype(str)) != names0:
            mismatches += 1
    if mismatches:
        print(f"Warning: method name/order mismatch in {mismatches} CSV(s); using the first CSV’s order.")
    return names0

def collect_crossing_years(
    csv_files: List[Path],
    method_names: List[str]
) -> Dict[str, Dict[str, List[float]]]:
    """
    Gather all crossing years from ALL CSVs into lists per method.
    Returns a dict like:
      {
        'xyear0.5': {'Method A': [years...], 'Method B': [...], ...},
        'xyear1.0': {...},
        'xyear1.5': {...}
      }
    """
    out: Dict[str, Dict[str, List[float]]] = {
        t: {m: [] for m in method_names} for t in THRESH_COLS
    }

    for p in csv_files:
        df = pd.read_csv(p)
        method_col = "method_name" if "method_name" in df.columns else df.columns[0]

        # Prefer named columns; if any missing, fall back to “last three columns”
        thresh_cols_present = [c for c in THRESH_COLS if c in df.columns]
        if len(thresh_cols_present) != 3:
            thresh_cols_present = list(df.columns[-3:])  # assume they are the last 3

        for _, row in df.iterrows():
            m = str(row[method_col])
            if m not in method_names:
                # Unknown method in this CSV; skip to keep ordering consistent
                continue
            for c in thresh_cols_present:
                v = row[c]
                if pd.notna(v):
                    try:
                        out[c][m].append(float(v))
                    except Exception:
                        pass
    return out

def _cdf_to_pdf(cdf: np.ndarray, axis: int = 1) -> np.ndarray:
    """
    Discrete PDF from CDF via forward difference along `axis`.
    Keeps the same length via 'prepend' and clips small negatives to 0.
    """
    diffs = np.diff(cdf, axis=axis, prepend=np.take(cdf, [0], axis=axis))
    return np.clip(diffs, 0, None)

def load_and_average_pdfs(npy_files: List[Path]) -> np.ndarray | None:
    """
    Load each npy (shape expected: [n_methods, n_fine, 3] CDF),
    convert to PDF, pad fine axis to the maximum length, and average across files.
    Returns mean_pdf with shape [n_methods, n_fine, 3] (nanmean over ensemble).
    """
    if not npy_files:
        return None

    pdfs = []
    n_methods_ref = None
    max_fine = 0

    # First pass: read and convert CDF → PDF, track shapes
    for p in npy_files:
        a = np.load(p, allow_pickle=False)
        if a.ndim != 3 or a.shape[2] != 3:
            print(f"Skipping {p.name}: unexpected shape {a.shape} (expected [n_methods, n_fine, 3])")
            continue
        pdf = _cdf_to_pdf(a, axis=1)  # diff along fine-time axis
        pdfs.append(pdf)
        n_methods_ref = a.shape[0] if n_methods_ref is None else min(n_methods_ref, a.shape[0])
        max_fine = max(max_fine, a.shape[1])

    if not pdfs:
        return None

    # Second pass: align shapes (truncate methods if needed, pad fine axis with NaN)
    aligned = []
    for pdf in pdfs:
        pdf = pdf[:n_methods_ref, :, :]  # truncate methods if needed
        pad_len = max_fine - pdf.shape[1]
        if pad_len > 0:
            pad = np.full((pdf.shape[0], pad_len, pdf.shape[2]), np.nan, dtype=pdf.dtype)
            pdf = np.concatenate([pdf, pad], axis=1)
        aligned.append(pdf)

    stack = np.stack(aligned, axis=0)  # [n_files, n_methods, n_fine, 3]
    mean_pdf = np.nanmean(stack, axis=0)  # [n_methods, n_fine, 3]
    return mean_pdf

def combine_histens_data(
    experiment_type: str,
    results_dir: str | Path = "Results"
) -> Tuple[List[str], np.ndarray | None, Dict[str, Dict[str, List[float]]]]:
    """
    High-level convenience wrapper:
      - prints counts of discovered files,
      - returns (method_names, mean_pdf, crossing_years).
    """
    csvs, npys = discover_histens_files(experiment_type, results_dir)
    method_names = read_method_names_from_csvs(csvs)
    crossing_years = collect_crossing_years(csvs, method_names)
    mean_pdf = load_and_average_pdfs(npys)

    if method_names:
        print(f"[{experiment_type}] methods: {len(method_names)}")
    if mean_pdf is not None:
        print(f"[{experiment_type}] mean_pdf shape: {mean_pdf.shape}  (methods × fine × thresholds)")
    else:
        print(f"[{experiment_type}] mean_pdf: None (no valid npy files)")

    # quick sanity on the crossing-year lists
    for t in THRESH_COLS:
        total_entries = sum(len(lst) for lst in crossing_years[t].values())
        print(f"[{experiment_type}] collected {total_entries} entries for {t}")
    return method_names, mean_pdf, crossing_years


method_names_rpi, mean_pdf_rpi, crossing_rpi = combine_histens_data("histens_rpi_real", results_dir="Results")
method_names_sat, mean_pdf_sat, crossing_sat = combine_histens_data("histens_satcal_real", results_dir="Results")

method_names = method_names_rpi
mean_cdfs=[ mean_pdf_rpi,  mean_pdf_sat]

# Plot (reuse earlier PANEL_INFO)
PANEL_INFO = [
    {"row_key":0, "threshold_idx":0, "title":"1850–1900 baseline — 0.5°C crossing (~1985)", "histogram":True},
    {"row_key":0, "threshold_idx":1, "title":"1850–1900 baseline — 1.0°C crossing (~2010)", "histogram":True},
    {"row_key":0, "threshold_idx":2, "title":"1850–1900 baseline — 1.5°C crossing ", "histogram":False},
    {"row_key":1, "threshold_idx":0, "title":"1981–2010 baseline — −0.25°C crossing (~1985)", "histogram":True},
    {"row_key":1, "threshold_idx":1, "title":"1981–2010 baseline — +0.25°C crossing (~2010)", "histogram":True},
    {"row_key":1, "threshold_idx":2, "title":"1981–2010 baseline — +0.75°C crossing ", "histogram":False},
]

fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex="col", sharey=False)
inum = 12
for j, info in enumerate(PANEL_INFO):
    row = j // 3; col = j % 3
    ax = axes[row, col]
    row_key = info["row_key"]; t_idx = info["threshold_idx"]
    ax.set_title(info["title"], fontsize=10)
    # PDFs (stacked)
    pdf = _cdf_to_pdf(mean_cdfs[row_key])
    n_methods = pdf.shape[0]; n_fine = pdf.shape[1]
    years = np.arange(evalmins[t_idx], evalmins[t_idx]+ (n_fine)*1/inum ,1/inum)
    y = np.nan_to_num(pdf[:,:,t_idx], nan=0.0, posinf=0.0, neginf=0.0)
    ax.stackplot(years, *y, labels=method_names[row_key], alpha=0.35, linewidth=0.0)
    for m in range(n_methods):
        ax.plot(years, y[m], linewidth=1.0)
    # Histograms when applicable
    if info["histogram"]:
        all_years = []
        for arr in hist_data[row_key][t_idx].values():
            all_years.extend(arr)
        if len(all_years) > 0:
            bins = np.arange(1960, 2051, 1)
            ax.hist(all_years, bins=bins, density=True, alpha=0.25, edgecolor="none")
    if row == 1:
        ax.set_xlabel("Year")
    if col == 0:
        ax.set_ylabel("Probability density")
    ax.grid(True, alpha=0.2, linewidth=0.5)

legend_labels = method_names["preind"]
handles = [plt.Line2D([0],[0]) for _ in legend_labels]
fig.legend(legend_labels, loc="upper center", ncol=min(len(legend_labels), 6), frameon=False)

fig.tight_layout(rect=[0,0,1,0.93])
outpath = "draft_heads_tails_probplot.png"
fig.savefig(outpath, dpi=600)
print(f"Saved figure: {outpath}")
