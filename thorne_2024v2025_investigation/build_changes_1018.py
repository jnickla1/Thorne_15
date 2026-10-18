"""Builds changes_1017_vs_1018.md: every number used in supplement_note_2024v2025.md, 17 Oct value vs 18 Oct value.
Old values come from the 17 Oct tables in this folder (or the 17 Oct note text where no table exists); new values from the *_1018 tables."""
import pickle, numpy as np, pandas as pd
from cfg_1018 import W, W_PAPER, REAL, ANTH
rows = []
def add(sec, q, old, new, dp, unit='°C', note=''):
    fmt = lambda v: '—' if v is None or (isinstance(v, float) and np.isnan(v)) else (f'{v:.{dp}f}' if isinstance(v, (int, float, np.floating)) else str(v))
    if isinstance(old, (int, float, np.floating)) and isinstance(new, (int, float, np.floating)):
        ch = '' if np.isnan(new - old) else f'{new - old:+.{dp + 1}f}'; txt = 'no' if (np.isnan(old) and np.isnan(new)) else ('yes' if round(old + 1e-12, dp) != round(new + 1e-12, dp) else 'no')
    else: ch = ''; txt = 'yes' if str(old) != str(new) else 'no'
    rows.append((sec, q, fmt(old), fmt(new), ch, txt, note))
rd = lambda f, **k: pd.read_csv(f, **k)
# ---------- headline
add('Headline', 'Realised PIVW 2025 level (rates/dashboard)', 1.3642, 1.3655464, 3, note='17 Oct text 1.364; rounds to 1.37 (was 1.36) at 2 dp')
add('Headline', 'Realised PIVW 2025, 2 dp as in table', 1.36, 1.37, 2, note='1.36555 → 1.37: conflicts with manuscript 1.36')
add('Headline', 'Realised 2025 interval lower', 1.230, 1.2286, 3, note='[1.23–1.53] at 2 dp unchanged')
add('Headline', 'Realised 2025 interval upper', 1.528, 1.5317, 3)
add('Headline', 'Realised change vs preprint (2 dp)', -0.04, -0.03, 2, note='1.37 − 1.40; unrounded 1.3655 − 1.406 = −0.041')
add('Headline', 'Anthropogenic PIVW 2025 level', 1.337, 1.3366395, 3, note='bit-identical inputs except FaIR_nonat (<3e-4)')
add('Headline', 'Anthropogenic 2025 interval', '[1.209, 1.496]', '[1.209, 1.496]', 3)
p15o = rd('p15_by_method_vintage.csv', header=[0, 1], index_col=[0, 1]); p15n = rd('p15_by_method_vintage_1018.csv', header=[0, 1], index_col=[0, 1])
for e, lab in [('rpi_real', 'Realised'), ('rpi_anthro', 'Anthropogenic')]:
    add('Headline', f'{lab} mean-member P(>1.5) 2025 (%)', 100 * p15o.loc[(e, 'PIVW_comb'), ('2025', '2025')], 100 * p15n.loc[(e, 'PIVW_comb'), ('2025', '2025')], 1, '%',
        'likelihood wording unchanged (very unlikely <10 %, extremely unlikely <5 %)')
    add('Headline', f'{lab} P(>1.5) 2024, 2024 vintage (%)', 100 * p15o.loc[(e, 'PIVW_comb'), ('2024', '2024')], 100 * p15n.loc[(e, 'PIVW_comb'), ('2024', '2024')], 1, '%')
    add('Headline', f'{lab} P(>1.5) 2024, new vintage (%)', 100 * p15o.loc[(e, 'PIVW_comb'), ('2025', '2024')], 100 * p15n.loc[(e, 'PIVW_comb'), ('2025', '2024')], 1, '%')
# ---------- weights
wo = dict(zip(REAL + ANTH, list(W_PAPER['real']) + list(W_PAPER['anthro']))); wn = dict(zip(REAL + ANTH, list(W['real']) + list(W['anthro'])))
for m in REAL + ANTH: add('Weights (PIVW %)', m, wo[m], wn[m], 1, '%', 'new = 1/err_var100 from Results2/historical_names_var_scale.csv (17 Oct 23:58)')
# ---------- decomposition
so = rd('step_2024to2025_headstails.csv', index_col=0); sn = rd('step_2024to2025_headstails_1018.csv', index_col=0)
i_old = {'real': 1.3606, 'anthro': 1.3550}; i_new = {'real': 1.3609, 'anthro': 1.3551}   # equiv_shift_spread(_1018).py direct estimate
rep = {'real': 1.406, 'anthro': 1.340}
for g, lab in [('real', 'Realised'), ('anthro', 'Anthropogenic')]:
    k = 'PIVW_' + g
    add('Decomposition', f'{lab}: (i) 2024 vintage, current endpoints', i_old[g], i_new[g], 3, note='18 Oct weights; paper weights give 1.3606 (realised)')
    add('Decomposition', f'{lab}: (i) step vs reported {rep[g]}', i_old[g] - rep[g], i_new[g] - rep[g], 3)
    add('Decomposition', f'{lab}: (ii) level, new vintage at 2024', so.loc[k, 'mean2024'], sn.loc[k, 'mean2024'], 3)
    add('Decomposition', f'{lab}: (ii) step', so.loc[k, 'mean2024'] - i_old[g], sn.loc[k, 'mean2024'] - i_new[g], 3)
    add('Decomposition', f'{lab}: (iii) 2025 level', so.loc[k, 'mean2025'], sn.loc[k, 'mean2025'], 3)
    add('Decomposition', f'{lab}: (iii) one-year step (mean)', so.loc[k, 'mean_step'], sn.loc[k, 'mean_step'], 3)
    add('Decomposition', f'{lab}: (iii) min member step', so.loc[k, 'min_step'], sn.loc[k, 'min_step'], 3)
    add('Decomposition', f'{lab}: (i)+(ii) shortfall', so.loc[k, 'mean2024'] - rep[g], sn.loc[k, 'mean2024'] - rep[g], 2)
for m in REAL + ANTH:
    add('Per-method 2025 step (heads/tails mean)', m, so.loc[m, 'mean_step'], sn.loc[m, 'mean_step'], 3, note=f"min member {sn.loc[m, 'min_step']:+.4f}; frac<0 {sn.loc[m, 'frac_members_neg']:.2f}")
for m in REAL + ANTH:
    add('Per-method 2024 level, new vintage', m, so.loc[m, 'mean2024'], sn.loc[m, 'mean2024'], 3)
co = rd('step2025_in_context.csv'); cn = rd('step2025_in_context_1018.csv')
for m in REAL + ANTH + ['PIVW_real', 'PIVW_anthro']:
    for thr in [-0.08, -0.1]:
        a = co[(co.method == m) & (co.thr == thr)].iloc[0]; b = cn[(cn.method == m) & (cn.thr == thr)].iloc[0]
        add('2025 step percentile (after annual drop ≤ thr)', f'{m}, thr {thr}', a.pct, b.pct, 0, 'pct', f'pct_all {a.pct_all:.0f} → {b.pct_all:.0f}; median step {a.med:+.4f} → {b.med:+.4f}')
# ---------- equivalent shifts
eo = rd('equivalent_level_shift_spread_adjusted.csv'); en = rd('equivalent_level_shift_spread_adjusted_1018.csv')
for _, r in en.iterrows():
    o = eo[(eo.exp == r.exp) & (eo.method == r.method)].iloc[0]
    add('Spread-adjusted equivalent shift', f'{r.exp} {r.method}', o.shift_spread_adjusted, r.shift_spread_adjusted, 3, note=f'unadjusted {o.shift_unadjusted:+.4f} → {r.shift_unadjusted:+.4f}')
add('Spread-adjusted equivalent shift', 'Realised PIVW implied 2024-vintage level (caveat)', 1.3356 + 0.02604, en[(en.exp == 'rpi_real') & (en.method == 'PIVW_comb')].implied_v24_level.iloc[0], 3)
add('Spread-adjusted equivalent shift', 'Anthro PIVW implied 2024-vintage level (caveat)', 1.3194 + 0.05110, en[(en.exp == 'rpi_anthro') & (en.method == 'PIVW_comb')].implied_v24_level.iloc[0], 3)
# ---------- HadCRUT single run
ho = rd('hadcrut_run_2024v2025.csv', index_col=0); hn = rd('hadcrut_run_2024v2025_1018.csv', index_col=0); ho.index = ho.index.str.replace('PIVW_anth$', 'PIVW_anthro', regex=True)
for m in REAL + ANTH + ['PIVW_real', 'PIVW_anthro']:
    for c in ['vintage_eff_2024', 'vintage_eff_2015', 'year_step_25']:
        if m in ho.index: add('HadCRUT single run (5.1 run − 5.0.2 run, same year)', f'{m} {c}', ho.loc[m, c], hn.loc[m, c], 3)
h24 = pickle.load(open('hist_results_93793643.pickle', 'rb')); e24 = np.asarray(h24['EBMKF_ta4']['LT_trend'][0], float)
for f, lab in [('hist_results_f36ba84b.pickle', 'old'), ('hist_results_1018.pickle', 'new')]:
    a = np.asarray(pickle.load(open(f, 'rb'))['EBMKF_ta4']['LT_trend'][0], float); globals()['eb_' + lab] = (a[173] - e24[173], a[170] - e24[170])
add('HadCRUT single run (5.1 run − 5.0.2 run, same year)', 'EBMKF_ta4 vintage_eff_2023 (2024 run ends 2023)', eb_old[0], eb_new[0], 3, note='breaks "all others within ±0.008"')
add('HadCRUT single run (5.1 run − 5.0.2 run, same year)', 'EBMKF_ta4 vintage_eff_2020', eb_old[1], eb_new[1], 3)
add('HadCRUT single run (5.1 run − 5.0.2 run, same year)', 'PIVW_real level 2024, 2025-vintage run', ho.loc['PIVW_real', 'v25_cur2024'], hn.loc['PIVW_real', 'v25_cur2024'], 3)
add('HadCRUT single run (5.1 run − 5.0.2 run, same year)', 'PIVW_real level 2024, 2024 run (text 1.360)', ho.loc['PIVW_real', 'v24_cur2024'], hn.loc['PIVW_real', 'v24_cur2024'], 3)
# ---------- weight sensitivity
ws = rd('weight_sensitivity_1018.csv').set_index('weights')
for lab_o, v_o, lab_n in [('Paper weights', 1.3642, '18 Oct var_scale (1/err_var100)'), ('Paper weights (17 Oct) on new levels', 1.3642, '17 Oct paper weights'),
                          ('1/RMS² Dec 2025 HadCRUT5.0', 1.364, '1/RMS^2 Dec2025 HadCRUT5.0'), ('1/RMS² HadCRUT5.1 (Oct 10 False → Oct 18 var_adjusted True)', 1.364, '1/RMS^2 Oct18 HadCRUT5.1 var_adjusted (True)'),
                          ('1/RMS² HadCRUT5.1 Oct 10 False', 1.364, '1/RMS^2 Oct10 HadCRUT5.1 (False)'), ('Equal weights', 1.373, 'equal')]:
    add('(b) Weight sensitivity: 2025 PIVW level', lab_o, v_o, ws.loc[lab_n, 'PIVW2025'], 3, note=f"EBMKF_ta4 weight {ws.loc[lab_n, 'EBMKF_ta4']:.1f} %")
add('(b) Weight sensitivity: 2025 PIVW level', 'Range of 2025 PIVW across non-equal weight sets', 0.000, np.ptp(ws.loc[[i for i in ws.index if i != 'equal'], 'PIVW2025'].to_numpy()), 3, note='"≤0.005 °C" claim holds')
# ---------- (d)
def load(s):
    a = rd(f'decrease_summary_hist_post1960_n3{s}.csv'); b = rd(f'decrease_summary_nonvolc_n3{s}.csv'); b = b[b.source.str.startswith(('MPI', 'NorESM'))]
    a['source'] = a.source.replace({'HadCRUT5.1': 'HadCRUT5.1 ≥1960', 'heads/tails': 'heads/tails ≥1960'}); return pd.concat([a, b]).set_index(['source', 'method'])
do, dn = load(''), load('_1018')
for (src, m), r in dn.iterrows():
    o = do.loc[(src, m)]
    add('(d) % of non-volcanic years with Δ<0', f'{m} | {src}', 100 * o.frac_neg, 100 * r.frac_neg, 0, '%')
    add('(d) % of non-volcanic years with Δ<−0.01', f'{m} | {src}', 100 * o.frac_lt_m01, 100 * r.frac_lt_m01, 1 if r.frac_lt_m01 < 0.01 else 0, '%')
    add('(d) mean size of decreases', f'{m} | {src}', o.mean_neg_size, r.mean_neg_size, 3, note=f'p05 {o.p05:+.4f} → {r.p05:+.4f}; spurious {100*o.frac_neg_spurious:.0f} → {100*r.frac_neg_spurious:.0f} %')
# cent20y check (benchmark is cent20y of the inputs; inputs unchanged)
add('(d) cent20y monotonicity', 'heads/tails member-years 1966–2015 decreasing (%)', 0.52, 0.52, 1, '%', 'from step2025_and_checks_1018.log; min step −0.0026 both')
add('(d) cent20y monotonicity', 'HadCRUT5.1 cent20y decreasing years 1966–2015', 0, 0, 0, 'count')
d = pd.DataFrame(rows, columns=['section', 'quantity', 'old (17 Oct)', 'new (18 Oct)', 'change', 'rounded text changes?', 'note'])
d.to_csv('changes_1017_vs_1018.csv', index=False)
with open('changes_1017_vs_1018.md', 'w') as f:
    f.write(open('changes_header_1018.md').read())
    for sec, g in d.groupby('section', sort=False):
        f.write(f'\n### {sec}\n\n| quantity | old (17 Oct) | new (18 Oct) | change | rounded text changes? | note |\n|---|---|---|---|---|---|\n')
        for _, r in g.iterrows(): f.write(f"| {r.quantity} | {r['old (17 Oct)']} | {r['new (18 Oct)']} | {r.change} | {'**yes**' if r['rounded text changes?']=='yes' else 'no'} | {r.note} |\n")
print(len(d), 'rows;', (d['rounded text changes?'] == 'yes').sum(), 'flagged'); print(d[d['rounded text changes?'] == 'yes'][['section', 'quantity', 'old (17 Oct)', 'new (18 Oct)']].to_string())
