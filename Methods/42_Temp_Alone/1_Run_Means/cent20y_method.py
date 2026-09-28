import numpy as np
import pandas as pd
import os

def forward_rolling(s, window):
    return s[::-1].rolling(window, min_periods=window).mean()[::-1]

def forward_rolling_std_of_means(df, window, ddof=1):
    means = df[::-1].rolling(window, min_periods=window).mean()[::-1]
    return means.std(axis=1, ddof=ddof)

def forward_rolling_sqrt_sumsq(s, window):
    return np.sqrt((s**2)[::-1].rolling(window, min_periods=window).sum()[::-1]) / window



def load_hadcrut_uncertainty(years, cur_path, window=20):
    version = "5.1.0.0" if years[-1] >= 2025 else "5.0.2.0"
    base = f"{cur_path}/../../../Common_Data/HadCRUT.{version}"

    summary   = pd.read_csv(f"{base}.analysis.summary_series.global.annual.csv")
    ensemble  = pd.read_csv(f"{base}.analysis.ensemble_series.global.annual.csv")

    mean_col     = summary.columns[1]     # Excel col B
    member_cols  = ensemble.columns[3:]   # Excel cols D:GU
    coverage_col = ensemble.columns[2]   # Excel col C

    mean_series = forward_rolling(summary[mean_col], window)
    ensemble_uncert = forward_rolling_std_of_means(ensemble[member_cols], window)
    coverage_uncert = forward_rolling_sqrt_sumsq(ensemble[coverage_col], window)

    total_uncert = np.sqrt(ensemble_uncert**2 + coverage_uncert**2)
    return total_uncert  # row-aligned with summary/ensemble/component



def run_method(years, temperature,uncert, model_run, experiment_type):
    avg_len_l=10
    avg_len_u=10 #actually 10 yrs below, that year, 9 yrs after
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    if experiment_type == 'historical':
        cur_path = os.path.dirname(os.path.realpath(__file__))
        hadcrutUncert = load_hadcrut_uncertainty(years, cur_path)

    for i in range(avg_len_l, len(years) - avg_len_u):
        chunka=temperature[i-avg_len_l:i+avg_len_u]
        chunkb=temperature[i-avg_len_l+1:i+avg_len_u+1]
        chunka_uncert=temps_1std[i-avg_len_l:i+avg_len_u]
        chunkb_uncert=temps_1std[i-avg_len_l+1:i+avg_len_u+1]
        means[i] = np.mean([chunka,chunkb])
        if experiment_type == 'historical':
            tot_uncerta = (hadcrutUncert[i-avg_len_l])**2
            tot_uncertb = (hadcrutUncert[i-avg_len_l+1])**2
        else:
            tot_uncerta= (np.var(chunka) + np.mean(chunka_uncert**2))/len(chunka)
            tot_uncertb =(np.var(chunkb) + np.mean(chunkb_uncert**2))/len(chunkb)
        
        ses[i] = np.sqrt((tot_uncerta+tot_uncertb )/2)
    
    return empser.copy(), empser.copy(), means, ses

