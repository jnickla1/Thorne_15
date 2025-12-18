from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import sys
#only does realized warming for the moment

from hist_heads_tails_evaluation_script import evalmins, evalmaxs, relthresh, sel_methods_list_real, sel_methods_list_anthro

goal=sys.argv[1]

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
    mean_cdf = np.nanmean(arr, axis=0)
    
    if mean_cdf.ndim != 3 or mean_cdf.shape[2] != 3:
        print(f"Error unexpected shape {a.shape} (expected [n_methods, n_yrs, 3])")
        exit()

    pdf = cdf_to_pdf(mean_cdf, axis=1)  # diff along fine-time axis
    
    return pdf,mean_cdf



    
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
        total_entries = sum(sum(np.array(lst)<2050) for lst in crossing_years[t].values())
        print(f"[{experiment_type}] collected {total_entries} entries for {t}")
    return method_names, mean_pdf,mean_cdf, crossing_years, indices_to_delete


method_names_rpi, mean_pdf_rpi,mean_cdf_rpi, crossing_rpi, del_idx = combine_histens_data("histens_rpi_"+goal, results_dir="Results")
method_names_sat, mean_pdf_sat,mean_cdf_sat, crossing_sat, del_idx = combine_histens_data("histens_satcal_"+goal, results_dir="Results")
method_names = method_names_rpi
mean_pdfs= [mean_pdf_rpi,  mean_pdf_sat]
mean_cdfs= [mean_cdf_rpi,  mean_cdf_sat]
hist_data = [crossing_rpi,crossing_sat]
# Plot (reuse earlier PANEL_INFO)
PANEL_INFO = [
    {"row_key":0, "threshold_idx":0, "title":"1850–1900 baseline: 0.5°C crossing (~1985)", "histogram":True},
    {"row_key":0, "threshold_idx":1, "title":"1850–1900 baseline: 1.0°C crossing (~2010)", "histogram":True},
    {"row_key":0, "threshold_idx":2, "title":"1850–1900 baseline: 1.5°C crossing ", "histogram":True},
    {"row_key":1, "threshold_idx":0, "title":f"1981–2010 baseline: {0.5+relthresh :.1f}°C crossing (~1985)", "histogram":True},
    {"row_key":1, "threshold_idx":1, "title":f"1981–2010 baseline: {1.0+relthresh :.1f}°C crossing (~2010)", "histogram":True},
    {"row_key":1, "threshold_idx":2, "title":f"1981–2010 baseline: {1.5+relthresh :.1f}°C crossing ", "histogram":True},
]

csv_threshold_names = ['xyear0.5','xyear1.0','xyear1.5']
fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex="col", sharey=False)

if goal=="real":
    sel_methods_list_sel = sel_methods_list_real
else:
    sel_methods_list_sel = sel_methods_list_anthro
    

parts = [p.split("/") for p in sel_methods_list_sel]

# verify last segment matches method_names_rpi (strip "_method.py")
last_parts = [seg[-1].replace("_method.py", "") for seg in parts]

for index in del_idx:
    del last_parts[index]
    del parts[index]

if goal=="anthro": last_parts.append('FaIR_anthro_unB_retro')

last_parts.append('equal_weight_comb'); last_parts.append('PIVW_comb')
print(method_names_rpi)
assert last_parts == method_names_rpi, "Method name mismatch!"

# build color keys: use "dir2/dir3" if present, otherwise just "dir2"
combined_keys = [("/".join(seg[1:3]) if len(seg) >= 4 else seg[1]) for seg in parts]

# combine 2nd and 3rd segments and get colors
from hist_evaluation_script import gen_color

sel_methods_colors = [gen_color(k) for k in combined_keys]

#if goal=="anthro": sel_methods_colors.append('#341539')

#if goal=="real":
#    sel_methods_colors.append('#F527E7')
#    sel_methods_colors.append('#27F546')


ylims=[0.4,0.4,0.4]
evalmins2=evalmins.copy()
evalmins2[1]=1995
import matplotlib.ticker as ticker
inum = 1 #adding 12 subdivisions didn't do anything since we linearly interpolated then took a diff

keep = [i for i, method in enumerate(method_names) if not (method.endswith("comb") or method.startswith("FaIR_anthro"))]

for j, info in enumerate(PANEL_INFO):
    row = j // 3; col = j % 3
    ax = axes[row, col]
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.set_xlim([evalmins2[col], evalmaxs[col]])
    ax.set_ylim([0,ylims[col]])
    row_key = info["row_key"]
    t_idx = info["threshold_idx"]
    ax.set_title(info["title"], fontsize=10)
    # PDFs (stacked)
    pdf = mean_pdfs[row_key]
    cdf = mean_cdfs[row_key]
    n_methods = pdf.shape[0]; n_yrs = pdf.shape[1]
    years = make_time_axis(evalmins[t_idx], n_yrs, inum=inum)
    #y = np.nan_to_num(pdf[:,:,t_idx], nan=0.0, posinf=0.0, neginf=0.0)
    y = pdf[:,:,t_idx]
    yc=cdf[:,:,t_idx]
        # ---- write panel CSV ----
    calibration = "satcal" if row else "rpi"
    outname = f"rates_histens/{calibration}_{goal}_thresh{(col*5):02d}.csv"

    # Build header: "method_name, years[0], years[1], ..."
    header = ["method_name"] + years.tolist()

    # Build rows: each is [method_name[m], yc[m,0], yc[m,1], ...]
    rows = [
        [method_names[m], *yc[m]]
        for m in range(n_methods)
    ]

    df_out = pd.DataFrame(rows, columns=header)
    df_out.to_csv(outname, index=False)
    #ax.stackplot(years, *y, labels=method_names[row_key], alpha=0.35, linewidth=0.0)
    for m in range(n_methods):
        #pdf_year_smooth = gaussian_filter1d(y[m], sigma=1.0, mode="nearest")
        tfine,pdffine = refine_pdf_repeat_gaussian(y[m], evalmins[col], L=12, sigma_years=1.0)
        #ax.plot(years,y[m] , linewidth=1.0, color= sel_methods_colors[m], label=method_names[m])
        if col==2:
            print(method_names[m])
            ind_latest=np.where(years==2024)[0]
            print(f"{(yc[m][ind_latest][0]):.2%} 1.5C crossing chance" )
        if(not (method_names[m].endswith("comb") or method_names[m].startswith("FaIR_anthro"))):
            ax.plot(tfine,pdffine , linewidth=1.0, color= sel_methods_colors[m], label=method_names[m])
    # Histograms when applicable
    if info["histogram"]:
        all_years = []
        this_hist_data= hist_data[row_key][csv_threshold_names[t_idx]]
        bins = np.arange(1960, 2051, 1)
        #for m in range(n_methods):
        #   ax.hist(this_hist_data[method_names[m]], bins=bins, density=True, color = sel_methods_colors[m], alpha=0.25, edgecolor="none")
         # Build data in the same order as method_names (so colors line up)
        data_list = [np.asarray(this_hist_data.get(m, []), dtype=float) for m in method_names]

        # Drop empty series to avoid blank legend entries (optional)
        #keep = [i for i, arr in enumerate(data_list) if arr.size > 0]
        
        if keep:
            data_list = [data_list[i] for i in keep]
            colors    = [sel_methods_colors[i] for i in keep]
            labels    = [method_names[i] for i in keep]

            ax.hist(data_list,bins=bins,density=True,stacked=True,
                color=colors,alpha=0.5,edgecolor="none",label=labels           # if you want a legend from the hist
            )
    if row == 1:
        ax.set_xlabel("Year")
    if col == 0:
        ax.set_ylabel("Probability density / year")

    if row==1 and col == 2 and goal=="anthro":
        ax.text(2015,0.2, 'No Detection', fontsize=18,ha='center',va='center')
        
    ax.grid(True, alpha=0.2, linewidth=0.5)

legend_labels = [method_names[i] for i in keep]
handles = [plt.Line2D([0],[0]) for _ in legend_labels]
fig.legend(legend_labels, loc="upper center", ncol=min(len(legend_labels), 6), frameon=False)

fig.tight_layout(rect=[0,0,1,0.93])
plt.show()
outpath = goal+"_heads_tails_probplot.png"
fig.savefig(outpath, dpi=600)
print(f"Saved figure: {outpath}")
