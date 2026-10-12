from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple



from evaluate_heads_tails_ensemble import evalmins, evalmaxs, relthresh, sel_methods_list_real, sel_methods_list_anthro

goals = ["real", "anthro"]

from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

import pandas as pd

THRESH_COLS = ["xyear0.5", "xyear1.0", "xyear1.5"]

from scipy.ndimage import gaussian_filter1d

def refine_pdf_repeat_gaussian(pdf_year, start_year, L=12, sigma_years=1.0):
    # pdf_year: 1D array over whole years, sum ≈ 1
    pdf_year = np.clip(np.asarray(pdf_year, float), 0, None)

        # ---- NEW: trim trailing zeros on the right ----
    nz = np.nonzero(pdf_year > 0)[0]
    if nz.size > 0:
        last = nz[-1]
        pdf_year = pdf_year[:last + 1]
    # If all zeros, leave as-is

    # Upsample: split each year into L bins; divide by L so mass per year is preserved
    pdf_fine = np.repeat(pdf_year / L, L)

    # Smooth on the fine grid; sigma must be in *fine bins*
    sigma_bins = sigma_years * L
    pdf_fine = gaussian_filter1d(pdf_fine, sigma=sigma_bins, mode="nearest")
    pdf_fine = np.clip(pdf_fine, 0, None)
    pdf_fine =  pdf_fine * L  # keep integral = L

    # Fine-grid time axis at bin centers
    dt = 1.0 / L
    
    t_fine = np.arange(start_year, start_year + len(pdf_year), dt) + 0.5 * dt

    pdf_fine = pdf_fine[:-int(L/2)]
    t_fine = t_fine[:-int(L/2)]

    return t_fine, pdf_fine


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
                        out[c][m].append(2050) #adding 2050s here to represent Nans so density=True still works!
                        pass
                else:
                    out[c][m].append(2050)
    return out

def make_time_axis(evalmin, n_yrs, inum=1):
    return evalmin + np.arange(n_yrs)/inum

def cdf_to_pdf(cdf: np.ndarray, axis: int = 1) -> np.ndarray:
    """
    Make the CDF monotone non-decreasing along `axis`, then take a discrete
    derivative to get a PDF. Length along `axis` is preserved via `prepend`.
    """
    cdf = np.asarray(cdf, dtype=float)
    # keep CDF in [0, 1]
    cdf = np.clip(cdf, 0.0, 1.0)

    # Treat NaNs as -inf so they don't raise the cummax; replace leading -inf with 0
    cdf_masked = np.where(np.isnan(cdf), -np.inf, cdf)
    cdf_mono = np.maximum.accumulate(cdf_masked, axis=axis)
    cdf_mono = np.where(np.isneginf(cdf_mono), 0.0, cdf_mono)
    cdf_mono = np.clip(cdf_mono, 0.0, 1.0)

    # Discrete derivative (same length) and non-negativity
    pdf = np.diff(cdf_mono, axis=axis, prepend=np.take(cdf_mono, [0], axis=axis))
    pdf = np.clip(pdf, 0.0, None)
    return pdf

def load_and_average_pdfs(npy_files: List[Path], delete_slices = []) -> np.ndarray | None:
    """
    Load each npy (shape expected: [n_methods, n_yrs, 3] CDF),
    convert to PDF, pad fine axis to the maximum length, and average across files.
    Returns mean_pdf with shape [n_methods, n_yrs, 3] (nanmean over ensemble).
    """
    if not npy_files:
        return None

    pdfs = []
    n_methods_ref = None
    max_fine = 0

    # First pass: read and convert CDF → PDF, track shapes
    arr0 = np.stack([np.load(p, allow_pickle=False) for p in npy_files], axis=0)
    arr = np.delete(arr0, delete_slices, axis=1)
    nan_n = np.isnan(arr).sum(axis=0)

    print("partially missing cells by threshold:",
          [np.sum((nan_n[:,:,j] > 0) & (nan_n[:,:,j] < arr.shape[0]))
           for j in range(3)])

    for j, name in enumerate(THRESH_COLS):
        print(name, "all-NaN indices (method, year-index):",
              np.argwhere(nan_n[:,:,j] == arr.shape[0]).tolist())
    
    allnan = np.all(np.isnan(arr), axis=0)
    print("all-NaN cells by threshold:", allnan.sum(axis=(0,1)))
    print("overall NaN fraction by threshold:", np.isnan(arr).mean(axis=(0,1,2)))
    
    mean_cdf = np.nanmean(arr, axis=0)
    
    if mean_cdf.ndim != 3 or mean_cdf.shape[2] != 3:
        print(f"Error unexpected shape {a.shape} (expected [n_methods, n_yrs, 3])")
        exit()

    pdf = cdf_to_pdf(mean_cdf, axis=1)  # diff along fine-time axis
    
    return pdf,mean_cdf



    
def combine_histens_data(goal:str,
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

    if goal=="real":
        to_remove = ["GWI_tot_SR15"]
    else:
        to_remove = ["cent20y","GWI_anthro_SR15","GWI_anthro"]
    indices_to_delete = []
    for index, item in enumerate(method_names):
        if item in to_remove:
            indices_to_delete.append(index)

    indices_to_delete.sort(reverse=True)
    for index in indices_to_delete:
        del method_names[index]
    
    crossing_years = collect_crossing_years(csvs, method_names)
    mean_pdf,mean_cdf = load_and_average_pdfs(npys, delete_slices=indices_to_delete)

    if method_names:
        print(f"[{experiment_type}] methods: {len(method_names)}")
    if mean_pdf is not None:
        print(f"[{experiment_type}] mean_pdf shape: {mean_pdf.shape}  (methods × fine × thresholds)")
    else:
        print(f"[{experiment_type}] mean_pdf: None (no valid npy files)")

    # quick sanity on the crossing-year lists
    for t in THRESH_COLS:
        allvals=crossing_years[t].values()
        allvalsarr = np.array([lst for lst in allvals])
        total_entries = sum(sum(lst<2050) for lst in allvalsarr[:-2])
        print(f"[{experiment_type}] collected {total_entries} entries for {t}")
        print(f"{method_names[-1]} collected {sum(allvalsarr[-1]<2050)} entries for {t}")
    return method_names, mean_pdf,mean_cdf, crossing_years, indices_to_delete



def load_goal(goal: str):
    method_names_rpi, mean_pdf_rpi, mean_cdf_rpi, crossing_rpi, del_idx = combine_histens_data(goal,
        "histens_rpi_" + goal, results_dir="Results"
    )
    method_names_sat, mean_pdf_sat, mean_cdf_sat, crossing_sat, del_idx2 = combine_histens_data(goal,
        "histens_satcal_" + goal, results_dir="Results"
    )

    # sanity: deletion indices should match; if not, keep the stricter intersection
    # (keeping simple here and trusting your existing pipeline)
    method_names = method_names_rpi
    mean_pdfs = [mean_pdf_rpi, mean_pdf_sat]   # row_key: 0=rpi, 1=satcal
    mean_cdfs = [mean_cdf_rpi, mean_cdf_sat]
    hist_data = [crossing_rpi, crossing_sat]

    return dict(
        goal=goal,
        method_names=method_names,
        mean_pdfs=mean_pdfs,
        mean_cdfs=mean_cdfs,
        hist_data=hist_data,
        del_idx=del_idx
    )

data_by_goal = {g: load_goal(g) for g in goals}

# Plot (reuse earlier PANEL_INFO)
##PANEL_INFO = [
##    {"row_key":0, "threshold_idx":0, "title":"1850–1900 baseline: 0.5°C crossing (~1985)", "histogram":True},
##    {"row_key":0, "threshold_idx":1, "title":"1850–1900 baseline: 1.0°C crossing (~2010)", "histogram":True},
##    {"row_key":0, "threshold_idx":2, "title":"1850–1900 baseline: 1.5°C crossing ", "histogram":True},
##    {"row_key":1, "threshold_idx":0, "title":f"1981–2010 baseline: {0.5+relthresh :.1f}°C crossing (~1985)", "histogram":True},
##    {"row_key":1, "threshold_idx":1, "title":f"1981–2010 baseline: {1.0+relthresh :.1f}°C crossing (~2010)", "histogram":True},
##    {"row_key":1, "threshold_idx":2, "title":f"1981–2010 baseline: {1.5+relthresh :.1f}°C crossing ", "histogram":True},
##]

csv_threshold_names = ['xyear0.5','xyear1.0','xyear1.5']
fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex="col", sharey=False)

from evaluate_historical_observations import gen_color

def build_method_colors_for_goal(goal: str, method_names: list[str], del_idx: list[int]) -> dict[str, str]:
    if goal == "real":
        sel_methods_list_sel = sel_methods_list_real
    else:
        sel_methods_list_sel = sel_methods_list_anthro

    parts = [p.split("/") for p in sel_methods_list_sel]
    last_parts = [seg[-1].replace("_method.py", "") for seg in parts]

    # apply same deletions as you do to method_names
    for index in sorted(del_idx, reverse=True):
        del last_parts[index]
        del parts[index]

    if goal == "anthro":
        last_parts.append('FaIR_anthro_unB_retro')

    last_parts.append('equal_weight_comb')
    last_parts.append('PIVW_comb')

    assert last_parts == method_names, "Method name mismatch!"

    combined_keys = [("/".join(seg[1:3]) if len(seg) >= 4 else seg[1]) for seg in parts]
    sel_methods_colors = [gen_color(k) for k in combined_keys]

    # NOTE: for the appended methods, pick explicit colors:
    # equal_weight_comb / PIVW_comb you can set as desired.
    # We'll only need PIVW here; equal_weight_comb is not plotted as a line in your current code.
    color_map = {m: c for m, c in zip(method_names, sel_methods_colors)}
    # If appended methods weren't in sel_methods_colors list, set fallbacks:
    color_map.setdefault("equal_weight_comb", "#777777")
    color_map.setdefault("PIVW_comb", "#444444")
    if goal == "anthro":
        color_map.setdefault("FaIR_anthro_unB_retro", "#555555")

    return color_map


color_maps = {
    g: build_method_colors_for_goal(g, data_by_goal[g]["method_names"], data_by_goal[g]["del_idx"])
    for g in goals
}


import matplotlib.ticker as ticker

csv_threshold_names = ['xyear0.5', 'xyear1.0', 'xyear1.5']
evalmins2 = evalmins.copy()
evalmins2[1] = 1995
ylims = [0.4, 0.4, 0.4]
inum = 1
bins = np.arange(1960, 2051, 1)

from plot_fut_results import sel_methods, sel_methodsA

dfw = pd.read_csv("Results2/historical_names_var_scale.csv")

for label, methods in [("Realized", sel_methods), ("Anthropogenic", sel_methodsA)]:
    x = dfw[dfw["method_name"].isin(methods)].copy()
    x["weight"] = 1 / x["err_var100"]
    x["percent"] = 100 * x["weight"] / x["weight"].sum()

    print(f"\nPIVW weights — {label}:")
    for name, pct in zip(x["method_name"], x["percent"]):
        print(f"  {name}: {pct:.2f}%")

for goal, methods in [("real", sel_methods), ("anthro", sel_methodsA)]:
    d = data_by_goal[goal]["hist_data"][0]["xyear1.5"]  # 1850–1900 baseline
    vals = np.concatenate([np.asarray(d[m]) for m in methods])
    print(goal, np.sum(vals < 2050), "/", len(vals))

    

def methods_for_main_stack(method_names: list[str]) -> list[int]:
    # mimic your existing logic
    return [i for i, m in enumerate(method_names) if not (m.endswith("comb") or m.startswith("FaIR_anthro"))]

def plot_baseline_row(ax_row, goal: str, baseline_row_key: int, suprow_label: str):
    """
    ax_row: list/array of 3 axes
    goal: "real" or "anthro"
    baseline_row_key: 0=rpi (1850–1900), 1=satcal (1981–2010)
    suprow_label: used only for any conditional annotations if needed
    """
    d = data_by_goal[goal]
    method_names = d["method_names"]
    pdf = d["mean_pdfs"][baseline_row_key]
    cdf = d["mean_cdfs"][baseline_row_key]
    hist_data = d["hist_data"][baseline_row_key]
    cmap = color_maps[goal]

    keep_idx = methods_for_main_stack(method_names)

    for col in range(3):
        ax = ax_row[col]
        ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.set_xlim([evalmins2[col], evalmaxs[col]])
        ax.set_ylim([0, ylims[col]])
        ax.grid(True, alpha=0.2, linewidth=0.5)

        t_idx = col
        n_methods = pdf.shape[0]
        years = make_time_axis(evalmins[t_idx], pdf.shape[1], inum=inum)

        # ---- lines (main set) ----
        for m in range(n_methods):
            name = method_names[m]
            tfine, pdffine = refine_pdf_repeat_gaussian(pdf[m, :, t_idx], evalmins[col], L=12, sigma_years=1.0)

            if name == "PIVW_comb":
                # special: thick black line
                ax.plot(tfine, pdffine, linewidth=3.0, color="k", label=name, zorder=5)
            elif not (name.endswith("comb") or name.startswith("FaIR_anthro")):
                ax.plot(tfine, pdffine, linewidth=1.0, color=cmap[name], label=name, zorder=4)

        # ---- stacked histogram (main set only; no legend labels here) ----
        this_hist = hist_data[csv_threshold_names[t_idx]]
        data_list = [np.asarray(this_hist.get(method_names[i], []), dtype=float) for i in keep_idx]
        colors = [cmap[method_names[i]] for i in keep_idx]

        if len(data_list) > 0:
            ax.hist(
                data_list,
                bins=bins,
                density=True,
                stacked=True,
                color=colors,
                alpha=0.5,
                edgecolor="none",
                zorder=1
            )

        # ---- PIVW histogram overlay (gap-filled dark gray, so colors show behind) ----
        if "PIVW_comb" in method_names:
            pivw_years = np.asarray(this_hist.get("PIVW_comb", []), dtype=float)
            if pivw_years.size > 0:
                ax.hist(
                    pivw_years,
                    bins=bins,
                    density=True,
                    histtype="bar",
                    color="#444444",
                    alpha=0.6,
                    rwidth=0.4,      # <-- gaps
                    edgecolor="none",
                    zorder=3
                )


def make_figure_for_baseline(baseline_row_key: int, suptitle: str, outpath: str):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5), sharex="col", sharey=False)

    # top = real, bottom = anthro
    plot_baseline_row(axes[0, :], goal="real",   baseline_row_key=baseline_row_key, suprow_label="real")
    plot_baseline_row(axes[1, :], goal="anthro", baseline_row_key=baseline_row_key, suprow_label="anthro")

    # per-panel titles (now baseline-specific)
    #
    if baseline_row_key==0:
        titles = [
            "0.5°C crossing (~1985)",
            "1.0°C crossing (~2010)",
            "1.5°C crossing",
        ]
    else:
        titles = [
            f"{0.5+relthresh :.1f}°C crossing (~1985)",
            f"{1.0+relthresh :.1f}°C crossing (~2010)",
            f"{1.5+relthresh :.1f}°C crossing",
        ]

    for c in range(3):
        axes[0, c].set_title("Realized warming: "+titles[c], fontsize=10)
        axes[1, c].set_title("Anthropogenic warming: " + titles[c], fontsize=10)

    for r in range(2):
        for c in range(3):
            if r == 1:
                axes[r, c].set_xlabel("Year")
            if c == 0:
                axes[r, c].set_ylabel("Probability density / year")

    # your "No Detection" annotation (anthro only; 1981–2010 baseline only)
    if baseline_row_key == 1:
        axes[1, 2].text(2015, 0.2, "No Detection", fontsize=18, ha="center", va="center")

    fig.suptitle(suptitle, y=0.98)

    # legend: create composite (hist patch + line) entries
    # We'll use REAL’s main methods for ordering, and include PIVW at the end if present.
    real_names = data_by_goal["real"]["method_names"]
    keep_real = methods_for_main_stack(real_names)
    legend_names = [real_names[i] for i in keep_real if real_names[i] != "PIVW_comb"]
    if "PIVW_comb" in real_names:
        legend_names.append("PIVW_comb")

    handles = []
    for name in legend_names:
        if name == "PIVW_comb":
            patch = mpatches.Patch(facecolor="#444444", alpha=0.6, edgecolor="none")
            line = Line2D([0], [0], color="k", lw=3.0)
        else:
            col = color_maps["real"].get(name, "#333333")
            patch = mpatches.Patch(facecolor=col, alpha=0.5, edgecolor="none")
            line = Line2D([0], [0], color=col, lw=1.5)
        handles.append((patch, line))

    # --- TOP LEGEND (nudged down from suptitle) ---
    top_legend = fig.legend(
        handles,
        legend_names,
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc="upper center",
        bbox_to_anchor=(0.5, 0.93),   # <-- push below suptitle
        ncol=min(len(legend_names), 6),
        frameon=False
    )

    anthro_names = data_by_goal["anthro"]["method_names"]
    keep_anthro = methods_for_main_stack(anthro_names)
    legend_names_anthro = [anthro_names[i] for i in keep_anthro if anthro_names[i] != "PIVW_comb"]
    if "PIVW_comb" in anthro_names:
        legend_names_anthro.append("PIVW_comb")

    handles = []
    for name in legend_names_anthro:
        if name == "PIVW_comb":
            patch = mpatches.Patch(facecolor="#444444", alpha=0.6, edgecolor="none")
            line = Line2D([0], [0], color="k", lw=3.0)
        else:
            col = color_maps["anthro"].get(name, "#333333")
            patch = mpatches.Patch(facecolor=col, alpha=0.5, edgecolor="none")
            line = Line2D([0], [0], color=col, lw=1.5)
        handles.append((patch, line))


    # --- BOTTOM LEGEND (below bottom row) ---
    bottom_legend = fig.legend(
        handles,
        legend_names_anthro,
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),   # <-- below panels
        ncol=min(len(legend_names), 6),
        frameon=False
    )
    for ax in axes.flat:
        ax.tick_params(axis="x", labelbottom=True)
    fig.tight_layout(rect=[0, 0.06, 1, 0.88])
    fig.savefig(outpath, dpi=600)
    print(f"Saved figure: {outpath}")
    return fig

make_figure_for_baseline(
    baseline_row_key=0,
    suptitle="1850–1900 Baseline",
    outpath="baseline_1850-1900_heads_tails_probplot.png"
)

make_figure_for_baseline(
    baseline_row_key=1,
    suptitle="1981–2010 Baseline",
    outpath="baseline_1981-2010_heads_tails_probplot.png"
)

plt.show()

#outpath = goal+"_heads_tails_probplot.png"
#fig.savefig(outpath, dpi=600)
#print(f"Saved figure: {outpath}")
