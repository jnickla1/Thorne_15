"""Check the Oscar leave-one-out outputs and build Table S4 (Mac side).

usage (from the repo root):
  python fut_LOO/check_and_aggregate.py fut_LOO/results/Results3_fix_oscar fut_LOO/output [MODE2_DIR [FAIR_DIR]]
MODE2_DIR holds the individual-reference results (nmes2run_*.npz: EBMKF_ta4, GWI_tot_CGWL); default = the first
directory. FAIR_DIR holds nmesFaIRrun_*.npz (FaIR_comb_unB, from fut_LOO/fut_fair_worker.py); default = MODE2_DIR.
Each complete set adds one block per method to table-S4-for-docs and to table-S4-individual-references.csv.

1. Coverage: every (mode in 7,16,32,8) x 270 member file exists and all fields are finite.
2. Cross-check: cases already computed by the stopped local pipeline run must agree.
3. Aggregation: same formulas as combination_pipeline.aggregate (legacy 90/65-year KL normalization and
   per-valid-year normalization); writes final_combined_{mode}methods_wide_ff.csv, table-S4-numeric.csv,
   table-S4-formatted.csv and table-S4.md into OUTDIR.
4. Paste-ready copies (table-S4-for-docs.xlsx / .csv, readable labels) and a comparison with the Dec 2025
   outputs in fut_LOO/legacy/dec2025_outputs (comparison_old_vs_new.csv).
"""
import sys, json, glob
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import gaussian_kde

SCENARIOS = ['fut_ESM1-2-LR_SSP126_constVolc', 'fut_ESM1-2-LR_SSP245_constVolc', 'fut_ESM1-2-LR_SSP370_constVolc',
             'fut_NorESM_RCP45_Volc', 'fut_NorESM_RCP45_VolcConst']
COUNTS = [50, 50, 50, 60, 60]; MODES = [7, 16, 32, 8]
FIELDS = ['firstcross15_diff', 'ncrosses', 'rmse', 'kl', 'rmse75', 'kl75']
HERE = Path(__file__).resolve().parent
# Optional cross-check source: the stopped 2026-10-19 Mac pipeline run (skipped if absent).
LOCAL = Path('/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/'
             'future-refresh-newsamples-fix-2026-10-19/combination-run')
LABELS = {'fut_ESM1-2-LR_SSP126_constVolc': 'MPI-ESM1-2-LR SSP1-2.6', 'fut_ESM1-2-LR_SSP245_constVolc': 'MPI-ESM1-2-LR SSP2-4.5',
          'fut_ESM1-2-LR_SSP370_constVolc': 'MPI-ESM1-2-LR SSP3-7.0', 'fut_NorESM_RCP45_Volc': 'NorESM RCP4.5 (volcanic)',
          'fut_NorESM_RCP45_VolcConst': 'NorESM RCP4.5 (constant volcanic)'}

src, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
mode2_src = Path(sys.argv[3]) if len(sys.argv) > 3 else src
fair_src = Path(sys.argv[4]) if len(sys.argv) > 4 else mode2_src
selections = json.loads((HERE/'selection.json').read_text())

# 1. coverage
missing, bad = [], []
for mode in MODES:
    for s, n in enumerate(COUNTS):
        for m in range(n):
            p = src/f'nmes{mode}run_exp{s}_r{m}.npz'
            if not p.exists(): missing.append(p.name); continue
            with np.load(p) as z:
                if any(z[k].shape != (2,) or not np.isfinite(z[k]).all() for k in FIELDS): bad.append(p.name)
print(f'coverage: {1080-len(missing)}/1080 present, {len(bad)} with non-finite fields')
if missing: print('  first missing:', missing[:10])
if bad: print('  bad:', bad[:10])

# 2. cross-check against the local pipeline cases
# CSCM (index 1) agrees to cross-platform floating-point precision; PIVW (index 0) also shifts because the local run used the
# stored 2000-2090 MPI err_var100, while this run uses 2000-2089 for all ensembles.
worst = {k: [0.0, 0.0] for k in FIELDS}; n = 0
for p in sorted(glob.glob(str(LOCAL/'run/Results3/nmes*run_exp*_r*.npz'))) if LOCAL.exists() else []:
    q = src/Path(p).name
    if not q.exists() or Path(p).name.startswith('nmes2run'): continue
    with np.load(p) as a, np.load(q) as b:
        for k in FIELDS:
            for j in range(2): worst[k][j] = max(worst[k][j], float(np.abs(a[k][j] - b[k][j])))
    n += 1
print(f'cross-check vs local pipeline: {n} overlapping cases')
print('  CSCM max abs diff (Mac vs Oscar: float-level, expect <0.5% in KL):', {k: f'{v[1]:.2e}' for k, v in worst.items()})
print('  PIVW max abs diff (expect small):', {k: f'{v[0]:.2e}' for k, v in worst.items()})
if missing or bad:
    sys.exit('Not aggregating: incomplete or invalid outputs.')

# 3. aggregate (mirrors combination_pipeline.aggregate)
rows, actual = [], []
for mode in MODES:
    for s, n in enumerate(COUNTS):
        arrays = []
        for i in range(n):
            with np.load(src/f'nmes{mode}run_exp{s}_r{i}.npz') as z: arrays.append({k: z[k].copy() for k in z.files})
        for method in range(2):
            row = {'experiment': SCENARIOS[s], 'approach': ['inverse_variance', 'sharpened_blended'][method],
                   'mode': mode, 'method_count': selections[str(mode)]['count']}
            kde = gaussian_kde(np.array([a['firstcross15_diff'][method] for a in arrays]))
            for radius in [1, 2, 5]: row[f'firstcross{radius}'] = float(kde.integrate_box_1d(-radius, radius))
            row['ncrosses'] = float(np.mean([a['ncrosses'][method] for a in arrays]))
            for key in ['rmse', 'rmse75']: row[key] = float(np.sqrt(np.mean([a[key][method]**2 for a in arrays])))
            for key, denom in [('kl', 90), ('kl75', 65)]: row[key] = float(np.mean([a[key][method] for a in arrays])/denom)
            rows.append(row)
            other = dict(row)
            for key, count in [('kl', 'n_valid_kl'), ('kl75', 'n_valid_kl75')]:
                other[key] = float(np.mean([a[key][method]/int(a[count]) for a in arrays]))
            actual.append(other)
frame = pd.DataFrame(rows); frame.to_csv(out/'combination-summary-legacy-normalization.csv', index=False)
pd.DataFrame(actual).to_csv(out/'combination-summary-per-valid-year.csv', index=False)
for mode in MODES: frame[frame['mode'] == mode].to_csv(out/f'final_combined_{mode}methods_wide_ff.csv', index=False)
table = frame[frame['mode'].isin([7, 16, 32])].copy()
table['combination'] = table.approach.map({'inverse_variance': 'PIVW', 'sharpened_blended': 'CSCM'})
cols = ['combination', 'method_count', 'mode', 'experiment', 'firstcross1', 'firstcross2', 'firstcross5', 'kl', 'kl75', 'ncrosses', 'rmse', 'rmse75']
table[cols].to_csv(out/'table-S4-numeric.csv', index=False)
fmt = table[cols].copy()
for c in ['firstcross1', 'firstcross2', 'firstcross5']: fmt[c] = fmt[c].map(lambda x: f'{100*x:.1f}%')
for c in ['kl', 'kl75', 'ncrosses']: fmt[c] = fmt[c].map(lambda x: f'{x:.3f}')
for c in ['rmse', 'rmse75']: fmt[c] = fmt[c].map(lambda x: f'{x:.4f}')
fmt.to_csv(out/'table-S4-formatted.csv', index=False)
heads = ['Combination', 'Actual method count', 'Selection mode', 'Ensemble', 'Within 1 yr', 'Within 2 yr', 'Within 5 yr',
         'KL 2000–2089', 'KL 2025–2089', 'Crossings', 'RMSE 2000–2089', 'RMSE 2025–2089']
lines = ['# Supplement Table S4 replacement (newsamples scale fix + nsharp optimizer fix, Oscar run)', '',
         '| ' + ' | '.join(heads) + ' |', '| ' + ' | '.join(['---']*len(heads)) + ' |']
lines += ['| ' + ' | '.join(str(v) for v in r) + ' |' for r in fmt.itertuples(index=False, name=None)]
(out/'table-S4.md').write_text('\n'.join(lines) + '\n')
print('wrote', out/'table-S4.md')

# 4. paste-ready table (open the .xlsx in Excel or Google Sheets, copy the range, paste into Google Docs)
# Order: PIVW then CSCM; method count ascending; ensembles in Figure 24 order; blank row between blocks.
FIG24 = ['fut_ESM1-2-LR_SSP370_constVolc', 'fut_ESM1-2-LR_SSP245_constVolc', 'fut_ESM1-2-LR_SSP126_constVolc',
         'fut_NorESM_RCP45_VolcConst', 'fut_NorESM_RCP45_Volc']
docs = fmt.copy()
docs['_c'] = docs['combination'].map({'PIVW': 0, 'CSCM': 1}); docs['_e'] = docs['experiment'].map(FIG24.index)
docs = docs.sort_values(['_c', 'method_count', '_e']).drop(columns=['_c', '_e', 'mode'])
docs['experiment'] = docs['experiment'].map(LABELS)
blocks = [g for _, g in docs.groupby(['combination', 'method_count'], sort=False)]
blank = pd.DataFrame([[''] * docs.shape[1]], columns=docs.columns)
# Individual reference methods, one block per method, same ensemble order.
ind_rows, ind_note, nvalid = [], None, set()
for prefix, names, folder in [('nmes2run', ['EBMKF_ta4', 'GWI_tot_CGWL'], mode2_src), ('nmesFaIRrun', ['FaIR_comb_unB'], fair_src)]:
    m2 = [[folder/f'{prefix}_exp{s}_r{i}.npz' for i in range(n)] for s, n in enumerate(COUNTS)]
    if not all(p.exists() for files in m2 for p in files): continue
    for method, name in enumerate(names):
        for s in [FIG24.index(e) for e in FIG24]:
            sidx = SCENARIOS.index(FIG24[s])
            arrays = []
            for p in m2[sidx]:
                with np.load(p) as z:
                    arrays.append({k: z[k].copy() for k in z.files})
                    if 'n_valid_kl' in z.files: nvalid.add((sidx < 3, int(z['n_valid_kl'])))
            kde = gaussian_kde(np.array([a['firstcross15_diff'][method] for a in arrays]))
            r = {'combination': name, 'method_count': 1, 'experiment': FIG24[s]}
            for radius in [1, 2, 5]: r[f'firstcross{radius}'] = f"{100*kde.integrate_box_1d(-radius, radius):.1f}%"
            r['kl'] = f"{np.mean([a['kl'][method] for a in arrays])/90:.3f}"
            r['kl75'] = f"{np.mean([a['kl75'][method] for a in arrays])/65:.3f}"
            r['ncrosses'] = f"{np.mean([a['ncrosses'][method] for a in arrays]):.3f}"
            r['rmse'] = f"{np.sqrt(np.mean([a['rmse'][method]**2 for a in arrays])):.4f}"
            r['rmse75'] = f"{np.sqrt(np.mean([a['rmse75'][method]**2 for a in arrays])):.4f}"
            ind_rows.append(r)
if ind_rows:
    if any(mpi and n != 90 for mpi, n in nvalid):
        ind_note = ('Note: the EBMKF_ta4 / GWI_tot_CGWL rows come from a run that scored MPI-ESM1-2-LR over 2000–2090 '
                    '(91 years); they will be regenerated with the 2000–2089 window.')
    ind = pd.DataFrame(ind_rows)
    ind[['combination', 'experiment', 'firstcross1', 'firstcross2', 'firstcross5', 'kl', 'kl75', 'ncrosses', 'rmse', 'rmse75']].to_csv(
        out/'table-S4-individual-references.csv', index=False)
    ind['experiment'] = ind['experiment'].map(LABELS)
    ind = ind[docs.columns]
    blocks += [g for _, g in ind.groupby('combination', sort=False)]
docs = pd.concat([x for g in blocks for x in (g, blank)][:-1], ignore_index=True)
if ind_note:
    docs = pd.concat([docs, blank, pd.DataFrame([[ind_note] + [''] * (docs.shape[1] - 1)], columns=docs.columns)], ignore_index=True)
docs.columns = ['Combination', 'Methods', 'Ensemble', 'Within 1 yr', 'Within 2 yr', 'Within 5 yr', 'KL 2000–2089',
                'KL 2025–2089', 'Crossings', 'RMSE 2000–2089 (°C)', 'RMSE 2025–2089 (°C)']
docs.to_csv(out/'table-S4-for-docs.csv', index=False, encoding='utf-8-sig')
docs.to_excel(out/'table-S4-for-docs.xlsx', index=False)

# 5. comparison with the Dec 2025 outputs
old_dir = HERE/'legacy/dec2025_outputs'; rows = []
for mode in [7, 16, 32]:
    p = old_dir/f'final_combined_{mode}methods_wide_ff.csv'
    if not p.exists(): continue
    j = frame[frame['mode'] == mode].merge(pd.read_csv(p), on=['experiment', 'approach'], suffixes=('', '_old'))
    for _, r in j.iterrows():
        rows.append(dict(mode=mode, methods=r.method_count, ensemble=LABELS[r.experiment],
                         combination='PIVW' if r.approach == 'inverse_variance' else 'CSCM',
                         kl_old=r.kl_old, kl_new=r.kl, kl75_old=r.kl75_old, kl75_new=r.kl75, rmse_old=r.rmse_old, rmse_new=r.rmse,
                         within1_old=r.firstcross1_old, within1_new=r.firstcross1, crossings_old=r.ncrosses_old, crossings_new=r.ncrosses))
if rows:
    pd.DataFrame(rows).to_csv(out/'comparison_old_vs_new.csv', index=False)
print('wrote', out/'table-S4-for-docs.xlsx', 'and', out/'comparison_old_vs_new.csv')
