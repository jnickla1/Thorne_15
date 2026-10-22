"""(e) Baseline hypothesis: does the switch of the pre-industrial datum from 1850-1899 to 1850-1900 (commits 216d3fde,
acf2df22, 2d95aaec; 21 Oct 2026) explain any of the 2024->2025 change?
1. Datum of the inputs: heads/tails members (both vintages) vs HadCRUT5.1.
2. 2025-vintage heads/tails ensemble-mean levels per method through four states:
   19 Oct (c5d9051e rates) -> seeded draws (0671536e rates) -> 1850-1900 code, no anchored correction (rates_histens 01:12)
   -> with anchored-method datum correction (current pickles, 21 Oct 02:18-03:35).
3. Same for the 2024 step and for the single HadCRUT5.1 run (hist_results_1019 vs _1021)."""
import pickle, numpy as np, pandas as pd
from cfg_1021 import RES, REAL, ANTH, W
out = []
# 1. input datum
h = pd.read_csv('../Common_Data/HadCRUT.5.1.csv').set_index('Time')['Anomaly']
print('HadCRUT5.1: mean(1850-1900) - mean(1850-1899) = %+.5f (datum_shift.HISTORICAL = 0.00254)' % (h.loc[1850:1900].mean() - h.loc[1850:1899].mean()))
for f, lab in [('../Hist_heads_tails/to_2024/sst_pseudo_rpi.csv', '2024 vintage'), ('../Hist_heads_tails/sst_pseudo.csv', '2025 vintage')]:
    e = np.loadtxt(f, delimiter=','); y = e[:, 0]; x = e[:, 1:]
    a = x[(y >= 1850) & (y <= 1899)].mean(0); b = x[(y >= 1850) & (y <= 1900)].mean(0)
    print(f'heads/tails {lab}: member mean over 1850-1900 = {b.mean():+.5f} (already on 1850-1900); over 1850-1899 = {a.mean():+.5f} '
          f'(range {a.min():+.4f}..{a.max():+.4f})')
# 2. ensemble-mean 2025 levels through the states
def rates(tag, exp):
    d = pd.read_csv(f'rates_states/rates_{exp}_{tag}.csv'); d = d[d.level == '2025'].set_index('method'); return d.iloc[:, 1].astype(float)
cur = {}
for exp, ms, g in [('rpi_real', REAL, 'real'), ('rpi_anthro', ANTH, 'anthro')]:
    L = {m: [] for m in ms}; L24 = {m: [] for m in ms}
    for k in range(100):
        d = pickle.load(open(RES + f'resultshistens_{exp}{k}.pickle', 'rb'))
        for m in ms:
            a = np.asarray(d[m]['LT_trend'][0], float); L[m].append(a[175]); L24[m].append(a[174])
    s = pd.Series({m: np.mean(L[m]) for m in ms}); s24 = pd.Series({m: np.mean(L24[m]) for m in ms})
    s['PIVW_comb'] = (s[ms] * W[g]).sum() / W[g].sum(); s24['PIVW_comb'] = (s24[ms] * W[g]).sum() / W[g].sum()
    st = pd.DataFrame({'19Oct': rates('c5d9051e_1019', exp), 'seeds_20Oct': rates('0671536e_seeds', exp),
                       'datum1900_nocorr_21Oct0112': rates('2d95aaec_0112', exp)}).loc[ms + ['PIVW_comb']]
    st['datum1900_corrected_21Oct'] = s
    st['seed_effect'] = st.seeds_20Oct - st['19Oct']; st['datum_switch_effect'] = st.datum1900_nocorr_21Oct0112 - st.seeds_20Oct
    st['anchored_correction_effect'] = st.datum1900_corrected_21Oct - st.datum1900_nocorr_21Oct0112
    st['total_19to21'] = st.datum1900_corrected_21Oct - st['19Oct']
    old = pd.read_csv('step_2024to2025_headstails_1019.csv', index_col=0).rename(index={'PIVW_' + g: 'PIVW_comb'})
    st['step2025_1019'] = old.loc[ms + ['PIVW_comb'], 'mean_step']; st['step2025_1021'] = s - s24
    print(f'\n{exp}: 2025 ensemble-mean level by state (PIVW rows use each state\'s own weights where from rates files)')
    print(st.round(5).to_string()); st.insert(0, 'exp', exp); out.append(st)
pd.concat(out).to_csv('baseline_hypothesis_1021.csv')
# 3. single run
h0 = pickle.load(open('hist_results_1019.pickle', 'rb')); h1 = pickle.load(open('hist_results_1021.pickle', 'rb'))
rows = []
for m in REAL + ANTH:
    a = np.asarray(h0[m]['LT_trend'][0], float); b = np.asarray(h1[m]['LT_trend'][0], float)
    rows.append(dict(method=m, d2024=b[174] - a[174], d2025=b[175] - a[175], d1950=b[100] - a[100], step_change=(b[175] - b[174]) - (a[175] - a[174])))
r = pd.DataFrame(rows).set_index('method'); print('\nsingle HadCRUT5.1 run, 21 Oct minus 19 Oct:'); print(r.round(5).to_string()); r.to_csv('baseline_hypothesis_hist_1021.csv')
