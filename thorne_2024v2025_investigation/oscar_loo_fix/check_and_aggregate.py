"""Check the Oscar leave-one-out outputs and build Table S4 (Mac side).

usage: python check_and_aggregate.py OSCAR_RESULTS3_FIX_DIR OUTDIR

1. Coverage: every (mode in 7,16,32,8) x 270 member file exists and all fields are finite.
2. Cross-check: cases already computed by the stopped local pipeline run must agree.
3. Aggregation: same formulas as combination_pipeline.aggregate (legacy 90/65-year KL normalization and
   per-valid-year normalization); writes final_combined_{mode}methods_wide_ff.csv, table-S4-numeric.csv,
   table-S4-formatted.csv and table-S4.md into OUTDIR.
"""
import sys, json, glob
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import gaussian_kde

SCENARIOS = ['fut_ESM1-2-LR_SSP126_constVolc', 'fut_ESM1-2-LR_SSP245_constVolc', 'fut_ESM1-2-LR_SSP370_constVolc',
             'fut_NorESM_RCP45_Volc', 'fut_NorESM_RCP45_VolcConst']
COUNTS = [50, 50, 50, 60, 60]; MODES = [7, 16, 32, 8]
FIELDS = ['firstcross15_diff', 'ncrosses', 'rmse', 'kl', 'rmse75', 'kl75']
LOCAL = Path('/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/'
             'future-refresh-newsamples-fix-2026-10-19/combination-run')

src, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
selections = json.loads((LOCAL/'selection.json').read_text())

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
worst = {k: 0.0 for k in FIELDS}; n = 0
for p in sorted(glob.glob(str(LOCAL/'run/Results3/nmes*run_exp*_r*.npz'))):
    q = src/Path(p).name
    if not q.exists() or Path(p).name.startswith('nmes2run'): continue
    with np.load(p) as a, np.load(q) as b:
        for k in FIELDS: worst[k] = max(worst[k], float(np.max(np.abs(a[k] - b[k]))))
    n += 1
print(f'cross-check vs local pipeline: {n} overlapping cases; max abs differences:', {k: f'{v:.2e}' for k, v in worst.items()})
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
         'KL 2000–2090', 'KL 2025–2090', 'Crossings', 'RMSE 2000–2090', 'RMSE 2025–2090']
lines = ['# Supplement Table S4 replacement (newsamples scale fix + nsharp optimizer fix, Oscar run)', '',
         '| ' + ' | '.join(heads) + ' |', '| ' + ' | '.join(['---']*len(heads)) + ' |']
lines += ['| ' + ' | '.join(str(v) for v in r) + ' |' for r in fmt.itertuples(index=False, name=None)]
(out/'table-S4.md').write_text('\n'.join(lines) + '\n')
print('wrote', out/'table-S4.md')
