"""Regenerate only the eROF rows of Results2/historical_newsamples.npy after the 2026-10-19
ppf_resample k=1 fix, using the stored best_alter_scale and the corrected (c + s*(r-c)) convention."""
import os, numpy as np, pandas as pd, importlib
data = pd.read_csv("./Common_Data/HadCRUT.5.1.csv")
temps_obs = data["Anomaly"].to_numpy(); temps_obs = temps_obs - np.mean(temps_obs[0:50])
uncert = (data["Lower"].to_numpy(), data["Upper"].to_numpy()); years = data["Time"].to_numpy()
df = pd.read_csv("Results2/historical_names_var_scale.csv", index_col=0)
C = np.load("Results2/historical_central_est.npy"); path = "Results2/historical_newsamples.npy"; S = np.load(path)
np.random.seed(20261019)
for m in ["eROF_anthro", "eROF_tot"]:
    i = list(df.method_name).index(m); s = df.best_alter_scale.iloc[i]
    res = importlib.import_module(f"Methods.43_Forcing_Based.3_Human_Induced.{m}_method").run_method(years, temps_obs, uncert, '', 'historical')
    c = res['mean'](years, 0)
    assert np.allclose(c, C[i], equal_nan=True), m          # same central estimate as stored
    r = res['resample'](years, S.shape[2], 0)
    S[i] = c[:, None] + (r - c[:, None]) * s
    print(m, "row", i, "scale", round(s, 3))
tmp = path[:-4] + ".tmp.npy"; np.save(tmp, S); os.replace(tmp, path)
