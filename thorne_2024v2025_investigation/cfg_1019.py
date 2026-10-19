"""Shared paths and weights for the *_1019 scripts (19 Oct 2026 vintage). Inputs frozen/recorded in inputs_1019_manifest.csv."""
import numpy as np, pandas as pd
ROOT = '../'
RES = ROOT + 'Results/'            # heads/tails pickles resultshistens_*{0..99}.pickle and future pickles (19 Oct)
RES2 = ROOT + 'Results2/'          # live: weights (var_scale) and comb scales are always read from the most current files
SNAP = 'snap25_1019'               # copies of Results/headstails_method_fineprobs_histens_* and headstails_statistics_histens_*
HIST = 'hist_results_1019.pickle'  # copy of ../hist_results.pickle (19 Oct 01:32)
HISTSTATS = ROOT + '261018_HadCRUT5.1_var_adjusted/current_methods_statistics_261019True.csv'
REAL = ['lowess1dt36wnc', 'GAM_AR1', 'Kal_flexLin', 'FaIR_comb_unB', 'EBMKF_ta4', 'GWI_tot_CGWL', 'CGWL10y_sfUKCP']
ANTH = ['FaIR_nonat_unB', 'GWI_anthro_CGWL']
VS = pd.read_csv(RES2 + 'historical_names_var_scale.csv').set_index('method_name')
def _w(ms):
    iv = 1 / VS.loc[ms, 'err_var100'].to_numpy(); return 100 * iv / iv.sum()
W = {'real': _w(REAL), 'anthro': _w(ANTH)}       # 1/err_var100, reproduces rates_histens PIVW_comb to 2e-6
W['anth'] = W['anthro']
W_PAPER = {'real': np.array([7.5, 7.4, 7.9, 24.6, 16.4, 25.7, 10.5]), 'anthro': np.array([69., 31.])}
SC = {'real': np.load(RES2 + 'comb_real_scales.npy'), 'anthro': np.load(RES2 + 'comb_anthro_scales.npy')}
# Record of what the 19 Oct tables were computed from: results2_1019/ holds copies of the files as of 19 Oct 18:50.
# If the live files have changed since, say so, so a rerun is not mistaken for the 19 Oct record.
import os, filecmp, sys
_changed = [f for f in ['historical_names_var_scale.csv', 'comb_real_scales.npy', 'comb_anthro_scales.npy']
            if os.path.exists('results2_1019/' + f) and not filecmp.cmp(RES2 + f, 'results2_1019/' + f, shallow=False)]
if _changed: print('cfg_1019: NOTE live Results2 differs from the 19 Oct snapshot in', _changed, file=sys.stderr)
