"""What changed between the 18 Oct and 19 Oct 2025-vintage outputs?
Heads/tails: the 18 Oct pickles were overwritten in Results/, so compare against the per-member 2024/2025 levels saved in
step_2024to2025_members_1018.csv. Single run: hist_results_1018.pickle vs hist_results_1019.pickle. Also weights."""
import pickle, numpy as np, pandas as pd
from cfg_1019 import RES, REAL, ANTH, W
import cfg_1018
old = pd.read_csv('step_2024to2025_members_1018.csv').set_index(['method', 'member'])
rows = []
for exp, ms in [('rpi_real', REAL), ('rpi_anthro', ANTH)]:
    for k in range(100):
        d = pickle.load(open(RES + f'resultshistens_{exp}{k}.pickle', 'rb'))
        for m in ms:
            a = np.asarray(d[m]['LT_trend'][0], float); o = old.loc[(m, k)]
            rows.append(dict(method=m, member=k, d2024=a[174] - o.lvl2024, d2025=a[175] - o.lvl2025))
d = pd.DataFrame(rows)
print('heads/tails current, 19 Oct minus 18 Oct (mean, max|.|):')
print(d.groupby('method')[['d2024', 'd2025']].agg(['mean', lambda x: x.abs().max()]).round(5).to_string())
h0 = pickle.load(open('hist_results_1018.pickle', 'rb')); h1 = pickle.load(open('hist_results_1019.pickle', 'rb'))
print('\nsingle run: methods whose current series changed (max abs, d2024, d2025)')
for m in sorted(set(h0) & set(h1)):
    try: x = np.asarray(h0[m]['LT_trend'][0], float); y = np.asarray(h1[m]['LT_trend'][0], float)
    except Exception: continue
    n = min(len(x), len(y)); df = y[:n] - x[:n]
    if np.isfinite(df).any() and np.nanmax(np.abs(df)) > 1e-9:
        print(f'  {m:28s} {np.nanmax(np.abs(df)):.4f}  d2024={df[174] if n>174 else np.nan:+.4f}  d2025={df[175] if n>175 else np.nan:+.4f}')
print('\nweights 18 Oct:', np.round(cfg_1018.W['real'], 2), np.round(cfg_1018.W['anthro'], 2))
print('weights 19 Oct:', np.round(W['real'], 2), np.round(W['anthro'], 2))
