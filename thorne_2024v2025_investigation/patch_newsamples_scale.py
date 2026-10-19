"""Exact post-hoc fix for the sampled-method newsamples scale bug (see cscm_newsamples_scale_bug.md).

Old files hold c + (r - c)/s; the corrected calibrate_method_uncertainty.py:380 gives c + (r - c)*s.
So new = c + (old - c)*s**2, applied only to rows of sampled (dict-returning) methods, with s and the
row order taken from the same run's _names_var_scale.csv.  Idempotent via a manifest.

usage: python patch_newsamples_scale.py historical [fut_..._Volc0 ...]
"""
import sys, os, json, numpy as np, pandas as pd

SAMPLED = {"CGWL10y_for_halfU", "CGWL10y_forec", "CGWL10y_pUKCP", "CGWL10y_sUKCP", "CGWL10y_sfUKCP",
           "FaIR_all", "FaIR_all_unB", "FaIR_anthro", "FaIR_anthro_unB", "FaIR_comb_unB", "FaIR_nonat",
           "FaIR_nonat_unB", "GWI_anthro", "GWI_anthro_AR6", "GWI_anthro_CGWL", "GWI_anthro_SR15",
           "GWI_anthro_orig", "GWI_tot", "GWI_tot_AR6", "GWI_tot_CGWL", "GWI_tot_SR15", "GWI_tot_orig",
           "KCC_all", "KCC_human", "eROF_anthro", "eROF_tot"}  # methods whose run_method returns a dict with 'resample'
MANIFEST = "Results2/newsamples_scale_patch_manifest.json"

done = json.load(open(MANIFEST)) if os.path.exists(MANIFEST) else {}
for name in sys.argv[1:]:
    if name in done:
        print(f"skip {name} (already patched)"); continue
    df = pd.read_csv(f"Results2/{name}_names_var_scale.csv", index_col=0)
    c = np.load(f"Results2/{name}_central_est.npy")
    path = f"Results2/{name}_newsamples.npy"
    S = np.load(path)
    assert S.shape[0] == len(df) == c.shape[0], (name, S.shape, len(df), c.shape)
    rows = [i for i, m in enumerate(df["method_name"]) if m in SAMPLED]
    for i in rows:
        s = df["best_alter_scale"].iloc[i]
        S[i] = c[i][:, None] + (S[i] - c[i][:, None]) * s**2
    tmp = path[:-4] + ".tmp.npy"
    np.save(tmp, S); os.replace(tmp, path)
    done[name] = {"rows": rows, "methods": [df["method_name"].iloc[i] for i in rows]}
    json.dump(done, open(MANIFEST, "w"), indent=1)
    print(f"patched {name}: {len(rows)} sampled rows")
