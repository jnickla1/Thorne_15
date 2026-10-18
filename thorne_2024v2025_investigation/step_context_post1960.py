"""Rebuilds step2025_in_context*.csv and decrease_summary_hist_post1960_n{3,5}*.csv from a steps_nonvolc file.
(The 17 Oct versions were produced inline; this script reproduces them exactly from the 17 Oct steps files — run with
argument '' to verify, '_1018' for the 18 Oct vintage.)
step2025_in_context: for each method, the per-member 2025 step is percentiled against all heads/tails non-volcanic
steps 1960-2024 (pct_all) and those following an annual drop <= thr (pct), averaged over members."""
import sys, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from decrease_analysis import summarise
suf = sys.argv[1] if len(sys.argv) > 1 else '_1018'; tag = suf if suf else '_verify1017'
d3 = pd.read_csv(f'steps_nonvolc_n3{suf}.csv.gz')
h = d3[d3.source == 'heads/tails']; rows = []
for m, g in h.groupby('method'):
    s25 = g[g.year == 2025].d_cur.to_numpy(); p = g[(g.year >= 1960) & (g.year < 2025)]
    pall = (p.d_cur.to_numpy()[:, None] < s25).mean() * 100
    for thr in [-0.08, -0.1]:
        c = p[p.d_raw <= thr].d_cur.to_numpy()
        rows.append(dict(method=m, thr=thr, n=len(c), med=np.median(c), pct=(c[:, None] < s25).mean() * 100, pct_all=pall, step25=s25.mean()))
pd.DataFrame(rows).to_csv(f'step2025_in_context{tag}.csv', index=False)
for nv in [3, 5]:
    d = pd.read_csv(f'steps_nonvolc_n{nv}{suf}.csv.gz')
    summarise(d[d.source.isin(['HadCRUT5.1', 'heads/tails']) & (d.year >= 1960)]).to_csv(f'decrease_summary_hist_post1960_n{nv}{tag}.csv')
if not suf:
    for a, b in [('step2025_in_context.csv', 'step2025_in_context_verify1017.csv')] + [(f'decrease_summary_hist_post1960_n{n}.csv', f'decrease_summary_hist_post1960_n{n}_verify1017.csv') for n in [3, 5]]:
        x = pd.read_csv(a); y = pd.read_csv(b); num = x.select_dtypes('number').columns
        print(a, 'reproduced:', x.shape == y.shape and np.allclose(x[num].to_numpy(float), y[num].to_numpy(float), equal_nan=True, atol=1e-12))
