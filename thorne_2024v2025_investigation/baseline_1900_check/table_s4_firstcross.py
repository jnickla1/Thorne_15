"""Table S4 'first crossing within 1/2/5 yr' under the 1850-1900 baseline (read-only check).

Each case in loo_shift/ was produced by LLOffcompute_shift_worker.py, which runs the unchanged LOO fit and then
evaluates the held-out member's first 1.5 °C crossing twice: at 1.5 (stored as firstcross15_diff, old baseline) and at
1.5 + c (firstcross15_diff_shift, new baseline; c = the member's 1850-1900 minus 1850-1899 offset). Shifting the
threshold by +c in the old frame is equivalent to rebaselining target and combined estimates by -c. PIVW and CSCM
weights, the PIVW predictive SD and nsharp are calibration-window statistics; for shift-equivariant inputs a common
constant shift leaves them unchanged, so the fit is reused.

usage (repo root): python thorne_2024v2025_investigation/baseline_1900_check/table_s4_firstcross.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import gaussian_kde
B = Path(__file__).resolve().parent; R = B.parents[1]
LABELS = {0: 'MPI-ESM1-2-LR SSP1-2.6', 1: 'MPI-ESM1-2-LR SSP2-4.5', 2: 'MPI-ESM1-2-LR SSP3-7.0', 3: 'NorESM RCP4.5 (volcanic)', 4: 'NorESM RCP4.5 (constant volcanic)'}
ORDER = [2, 1, 0, 4, 3]; COUNTS = [50, 50, 50, 60, 60]; COUNT = {7: 7, 16: 18, 32: 33}
rows, xcheck = [], []
for mode in [7, 16, 32]:
    for s in ORDER:
        files = [B/f'loo_shift/nmes{mode}run_exp{s}_r{m}.npz' for m in range(COUNTS[s])]
        if not all(f.exists() for f in files):
            continue
        z = [np.load(f) for f in files]
        old = np.array([x['firstcross15_diff'] for x in z]); new = np.array([x['firstcross15_diff_shift'] for x in z])
        stored = np.array([np.load(R/f'fut_LOO/results/Results3_fix_oscar/{f.name}')['firstcross15_diff'] for f in files])
        for j, comb in enumerate(['PIVW', 'CSCM']):
            r = dict(combination=comb, methods=COUNT[mode], ensemble=LABELS[s])
            for tag, arr in [('old', old[:, j]), ('new', new[:, j])]:
                kde = gaussian_kde(arr)
                for w in [1, 2, 5]:
                    r[f'within{w}_{tag}'] = 100*kde.integrate_box_1d(-w, w)
            r['members_changed'] = int((np.abs(new[:, j]-old[:, j]) > 1e-9).sum())
            r['max_abs_member_change_yr'] = float(np.abs(new[:, j]-old[:, j]).max())
            rows.append(r)
            xcheck.append(dict(combination=comb, methods=COUNT[mode], ensemble=LABELS[s],
                               members_differing_from_table_S4=int((np.abs(old[:, j]-stored[:, j]) > 1e-9).sum()),
                               max_abs_diff_yr=float(np.abs(old[:, j]-stored[:, j]).max())))
t = pd.DataFrame(rows)
for w in [1, 2, 5]: t[f'd_within{w}'] = t[f'within{w}_new'] - t[f'within{w}_old']
t = t.sort_values(['combination', 'methods'], key=lambda c: c.map({'PIVW': 0, 'CSCM': 1}) if c.name == 'combination' else c, kind='stable')
t.round(2).to_csv(B/'tableS4_firstcross_old_vs_new.csv', index=False)
pd.DataFrame(xcheck).to_csv(B/'tableS4_old_vs_stored_crosscheck.csv', index=False)
print(t[['combination', 'methods', 'ensemble', 'within1_old', 'within1_new', 'within2_old', 'within2_new', 'within5_old', 'within5_new', 'members_changed', 'max_abs_member_change_yr']].round(1).to_string(index=False))
print('\nmax |change| in percentage points: within1 %.2f, within2 %.2f, within5 %.2f' % tuple(t[f'd_within{w}'].abs().max() for w in [1, 2, 5]))
x = pd.DataFrame(xcheck); print('cross-check of the unshifted values against the current Table S4 inputs: max |diff| %.3f yr in %d member-rows' % (x.max_abs_diff_yr.max(), x.members_differing_from_table_S4.sum()))
