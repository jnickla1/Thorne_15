"""Shared paths and weights for the *_1018 scripts (18 Oct 2026 vintage). Inputs frozen/recorded in inputs_1018_manifest.csv."""
import numpy as np, pandas as pd
ROOT = '../'
RES = ROOT + 'Results/'            # heads/tails pickles resultshistens_*{0..99}.pickle and future pickles (18 Oct)
RES2 = ROOT + 'Results2/'          # historical_names_var_scale.csv, comb_*_scales.npy (17 Oct 23:58)
SNAP = 'snap25_1018'               # copies of Results/headstails_method_fineprobs_histens_* and headstails_statistics_histens_*
HIST = 'hist_results_1018.pickle'  # copy of ../hist_results.pickle (18 Oct 11:29)
HISTSTATS = ROOT + '261018_HadCRUT5.1_var_adjusted/current_methods_statistics_261018True.csv'
REAL = ['lowess1dt36wnc', 'GAM_AR1', 'Kal_flexLin', 'FaIR_comb_unB', 'EBMKF_ta4', 'GWI_tot_CGWL', 'CGWL10y_sfUKCP']
ANTH = ['FaIR_nonat_unB', 'GWI_anthro_CGWL']
VS = pd.read_csv(RES2 + 'historical_names_var_scale.csv').set_index('method_name')
def _w(ms):
    iv = 1 / VS.loc[ms, 'err_var100'].to_numpy(); return 100 * iv / iv.sum()
W = {'real': _w(REAL), 'anthro': _w(ANTH)}       # 1/err_var100, reproduces rates_histens PIVW_comb to 2e-6
W['anth'] = W['anthro']
W_PAPER = {'real': np.array([7.5, 7.4, 7.9, 24.6, 16.4, 25.7, 10.5]), 'anthro': np.array([69., 31.])}
SC = {'real': np.load(RES2 + 'comb_real_scales.npy'), 'anthro': np.load(RES2 + 'comb_anthro_scales.npy')}
