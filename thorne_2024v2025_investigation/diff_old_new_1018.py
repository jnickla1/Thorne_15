"""Which methods changed between the 17 Oct investigation's 2025 vintage (Oct-10 backup / f36ba84b) and the 18 Oct outputs?
Also: SSP126 vs SSP245/370 agreement over the shared historical period, and hist_results vs the var_adjusted export."""
import pickle, numpy as np, pandas as pd
BK = '/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-rerun-2026-10-13/before/Results/'
NEW = '../Results/'
REAL = ['lowess1dt36wnc', 'GAM_AR1', 'Kal_flexLin', 'FaIR_comb_unB', 'EBMKF_ta4', 'GWI_tot_CGWL', 'CGWL10y_sfUKCP', 'cent20y']; ANTH = ['FaIR_nonat_unB', 'GWI_anthro_CGWL']
rows = []
for exp, ms in [('rpi_real', REAL), ('rpi_anthro', ANTH)]:
    for k in range(100):
        a = pickle.load(open(BK + f'resultshistens_{exp}{k}.pickle', 'rb')); b = pickle.load(open(NEW + f'resultshistens_{exp}{k}.pickle', 'rb'))
        for m in ms:
            for j, lab in [(0, 'cur'), (1, 'se'), (2, 'retro')]:
                x = np.asarray(a[m]['LT_trend'][j], float); y = np.asarray(b[m]['LT_trend'][j], float)
                n = min(len(x), len(y)); dd = np.abs(y[:n] - x[:n]); dd = dd[np.isfinite(dd)]
                rows.append(dict(exp=exp, member=k, method=m, comp=lab, maxabs=dd.max() if len(dd) else np.nan,
                                 d2024=y[174] - x[174] if lab != 'retro' else np.nan, d2025=y[175] - x[175] if lab != 'retro' else np.nan))
d = pd.DataFrame(rows)
print('heads/tails pickles, new minus Oct-10 backup (max over members):')
print(d.groupby(['exp', 'method', 'comp'])[['maxabs', 'd2024', 'd2025']].agg(['max', 'mean']).round(5).to_string())
h0 = pickle.load(open('hist_results_f36ba84b.pickle', 'rb')); h1 = pickle.load(open('hist_results_1018.pickle', 'rb'))
print('\nhist_results f36ba84b -> 1018: methods with any change in current (max abs, 2024, 2025)')
for m in sorted(set(h0) & set(h1)):
    try:
        x = np.asarray(h0[m]['LT_trend'][0], float); y = np.asarray(h1[m]['LT_trend'][0], float)
    except Exception: continue
    n = min(len(x), len(y)); dd = np.nanmax(np.abs(y[:n] - x[:n])) if np.isfinite(y[:n] - x[:n]).any() else 0
    if dd > 1e-9: print(f'  {m:28s} {dd:.4f}  d2024={y[174]-x[174]:+.4f}  d2025={y[175]-x[175]:+.4f}')
print('only in old:', sorted(set(h0) - set(h1)), ' only in new:', sorted(set(h1) - set(h0)))
# SSP126 consistency with SSP245/370 over 1850-2014
for k in [0, 7, 23, 49]:
    P = {s: pickle.load(open(NEW + f'resultsfut_ESM1-2-LR_{s}_constVolc{k}.pickle', 'rb')) for s in ['SSP126', 'SSP245', 'SSP370']}
    out = []
    for m in REAL[:-1] + ANTH:
        a = np.asarray(P['SSP126'][m]['LT_trend'][0], float)[:165]
        out.append(f"{m}:{max(np.nanmax(np.abs(a - np.asarray(P[s][m]['LT_trend'][0], float)[:165])) for s in ['SSP245', 'SSP370']):.1e}")
    print('MPI member', k, 'SSP126 vs 245/370 max |diff| 1850-2014:', ' '.join(out))
# hist_results vs var_adjusted export
ex = np.load('../261018_HadCRUT5.1_var_adjusted/hist_results_export_2x176x125.npy'); man = pd.read_csv('../261018_HadCRUT5.1_var_adjusted/hist_results_manifest.csv')
print('\nexport', ex.shape, man.columns.tolist()[:6]); print(man.head(3).to_string())
