"""rms_z of sampled-method newsamples vs the run's own standard (2000-2089): ~1 if correctly scaled,
~s^2 if unpatched, ~s^4 if patched twice. usage: python check_sample_calibration.py RUNNAME [RUNNAME ...]"""
import sys, numpy as np, pandas as pd
SAMPLED = {"FaIR_comb_unB", "GWI_tot_CGWL", "CGWL10y_sfUKCP", "GWI_anthro_CGWL", "FaIR_nonat_unB", "CGWL10y_pUKCP"}
for name in sys.argv[1:]:
    df = pd.read_csv(f"Results2/{name}_names_var_scale.csv", index_col=0)
    S = np.load(f"Results2/{name}_newsamples.npy", mmap_mode="r"); C = np.load(f"Results2/{name}_central_est.npy")
    st = np.load(f"Results2/{name}_standard.npy"); sl = slice(150, 240); out = []
    for i, m in enumerate(df.method_name):
        if m not in SAMPLED: continue
        sd = np.nanstd(np.asarray(S[i][sl]), axis=1); ok = np.isfinite(sd) & (sd > 0) & np.isfinite(st[sl])
        z = (st[sl] - C[i][sl])[ok] / sd[ok]
        out.append(f"{m}:s={df.best_alter_scale.iloc[i]:.2f},rmsz={np.sqrt(np.mean(z**2)):.2f}")
    print(name, " ".join(out))
