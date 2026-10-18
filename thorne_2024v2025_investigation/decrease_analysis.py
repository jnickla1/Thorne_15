"""(d) How often, and by how much, does each selected method's CURRENT estimate decrease
from one year to the next outside volcanic periods?

Delta(t) = current(t) - current(t-1), where current(t) uses data only through year t.
A decrease is classed as "tracking" if the benchmark cent20y (retrospective centred
20-yr mean of the same series) also fell over t-1 -> t, otherwise "spurious".

Sources (all read-only snapshots; nothing in the main tree is modified):
  hist   : HadCRUT5.1 observational run, git f36ba84b hist_results.pickle
  histens: 100-member heads/tails ensemble, 2025 paper state (Oct-10 backup == f36ba84b CSVs)
  MPI SSP126/245/370 constVolc (50 each) and NorESM RCP45 VolcConst (60), git f36ba84b
Futures: only steps after the switch from historical forcing (MPI: Delta from 2016,
NorESM: Delta from 2007). Volcanic mask (historical/heads-tails only): AOD (9c column D,
x1e-3) > 0.015 -> that year plus the next N years excluded (N=3 main, N=5 sensitivity).
"""
import pickle, numpy as np, pandas as pd, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
BK = Path('/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-rerun-2026-10-13/before')
FUT = Path('/private/tmp/claude-501/-Users-JohnMatthew-Downloads-Thorne-15-codefigurestats/252b6cf2-e2ee-4c12-bbc1-45fc95f574e4/scratchpad/fut')

REAL = ['lowess1dt36wnc', 'GAM_AR1', 'Kal_flexLin', 'FaIR_comb_unB', 'EBMKF_ta4', 'GWI_tot_CGWL', 'CGWL10y_sfUKCP']
ANTH = ['FaIR_nonat_unB', 'GWI_anthro_CGWL']
W_REAL = np.array([7.5, 7.4, 7.9, 24.6, 16.4, 25.7, 10.5]); W_ANTH = np.array([69.0, 31.0])  # paper weights

aod = pd.read_csv(ROOT / 'Common_Data/toyKFmodelData9c.csv', header=None).set_index(0)[3]
erupt = aod.index[aod > 15].to_numpy()

def volc_mask(years, n):
    bad = np.zeros(len(years), bool)
    for e in erupt:
        bad |= (years >= e) & (years <= e + n)
    return bad

def cur(d, m):
    t = d[m]['LT_trend']
    return np.asarray(t[0], float)

def add_pivw(d, methods, w):
    st = np.array([cur(d, m) for m in methods]); ww = w[:, None] * np.isfinite(st)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.nansum(st * w[:, None], 0) / ww.sum(0)

def series(d, kind):
    """Return dict name -> current series, plus benchmark cent20y and raw annual."""
    s = {}
    meths = REAL + ANTH if kind == 'hist' else (REAL if kind != 'anthro' else ANTH)
    for m in meths:
        if m in d: s[m] = cur(d, m)
    if all(m in d for m in REAL): s['PIVW_real'] = add_pivw(d, REAL, W_REAL)
    if all(m in d for m in ANTH): s['PIVW_anthro'] = add_pivw(d, ANTH, W_ANTH)
    bench = np.asarray(d['cent20y']['LT_trend'][2], float)
    return s, bench

rows = []
def collect(source, member, s, bench, years, keep, raw=None):
    db = np.diff(bench)
    draw = np.diff(raw) if raw is not None else np.full(len(years) - 1, np.nan)
    for name, c in s.items():
        dc = np.diff(c)
        for i in np.where(keep[1:] & np.isfinite(dc))[0]:
            rows.append((source, member, name, int(years[i + 1]), dc[i], db[i], draw[i]))

def run(nvolc):
    rows.clear()
    # 1. HadCRUT5.1 single observational run
    h = pickle.load(open(OUT / 'hist_results_f36ba84b.pickle', 'rb'))
    yrs = np.arange(1850, 1850 + len(cur(h, 'lowess1dt36wnc')))
    raw = pd.read_csv(ROOT / 'Common_Data/HadCRUT.5.1.csv').set_index('Time')['Anomaly'].reindex(yrs).to_numpy()
    s, b = series(h, 'hist')
    collect('HadCRUT5.1', 0, s, b, yrs, ~volc_mask(yrs, nvolc), raw)
    # 2. heads/tails ensemble (realised and anthropogenic pickles), pre-industrial baseline
    ens = np.loadtxt(ROOT / 'Hist_heads_tails/sst_pseudo.csv', delimiter=',')
    for exp, kind in [('rpi_real', 'real'), ('rpi_anthro', 'anthro')]:
        for m in range(100):
            d = pickle.load(open(BK / f'Results/resultshistens_{exp}{m}.pickle', 'rb'))
            s, b = series(d, kind)
            collect('heads/tails', m, s, b, yrs, ~volc_mask(yrs, nvolc), ens[:, m + 1])
    # 3. futures, post-switch only, no volcanic mask needed (constant background volcanism)
    for e, n, y0 in [('ESM1-2-LR_SSP126_constVolc', 50, 2015), ('ESM1-2-LR_SSP245_constVolc', 50, 2015),
                     ('ESM1-2-LR_SSP370_constVolc', 50, 2015), ('NorESM_RCP45_VolcConst', 60, 2006)]:
        for m in range(n):
            d = pickle.load(open(FUT / f'{e}{m}.pickle', 'rb'))
            L = len(cur(d, 'lowess1dt36wnc')); fy = np.arange(1850, 1850 + L)
            s = {k: cur(d, k) for k in REAL + ANTH if k in d}
            s['PIVW_real'] = add_pivw(d, REAL, W_REAL); s['PIVW_anthro'] = add_pivw(d, ANTH, W_ANTH)
            b = np.asarray(d['cent20y']['LT_trend'][2], float)
            collect(e.replace('ESM1-2-LR_', 'MPI_').replace('_constVolc', ''), m, s, b, fy, fy >= y0)
    return pd.DataFrame(rows, columns=['source', 'member', 'method', 'year', 'd_cur', 'd_bench', 'd_raw'])

def summarise(df):
    def f(g):
        d = g.d_cur; neg = d < 0; hb = np.isfinite(g.d_bench)
        return pd.Series(dict(
            n_steps=len(d), frac_neg=neg.mean(), frac_lt_m01=(d < -0.01).mean(), frac_lt_m02=(d < -0.02).mean(),
            mean_step=d.mean(), p05=d.quantile(.05), p01=d.quantile(.01), min_step=d.min(),
            mean_neg_size=d[neg].mean() if neg.any() else np.nan,
            frac_neg_spurious=((neg & hb & (g.d_bench > 0)).sum() / max((neg & hb).sum(), 1)) if hb.any() else np.nan,
            frac_neg_given_rawdrop10=(d[g.d_raw <= -0.10] < 0).mean() if (g.d_raw <= -0.10).any() else np.nan,
            mean_step_given_rawdrop10=d[g.d_raw <= -0.10].mean() if (g.d_raw <= -0.10).any() else np.nan,
            n_rawdrop10=int((g.d_raw <= -0.10).sum())))
    return df.groupby(['source', 'method']).apply(f)

if __name__ == '__main__':
    print('eruption-flag years (AOD>0.015):', list(erupt))
    for nv in [3, 5]:
        df = run(nv)
        df.to_csv(OUT / f'steps_nonvolc_n{nv}.csv.gz', index=False)
        sm = summarise(df); sm.to_csv(OUT / f'decrease_summary_nonvolc_n{nv}.csv')
        pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
        print(f'\n=== volcanic exclusion: eruption year + {nv} ===')
        print(sm[['n_steps', 'frac_neg', 'frac_lt_m01', 'frac_lt_m02', 'mean_step', 'p05', 'min_step',
                  'frac_neg_spurious', 'frac_neg_given_rawdrop10', 'mean_step_given_rawdrop10', 'n_rawdrop10']].round(3))
