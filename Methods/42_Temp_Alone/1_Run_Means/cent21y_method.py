import numpy as np

def run_method(years, temperature, uncert, model_run, experiment_type):
    avg_len_l=10
    avg_len_u=11 #actually 10 yrs below, that year, 10 yrs after
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    for i in range(avg_len_l, len(years) - avg_len_u+1):
        chunk=temperature[i-avg_len_l:i+avg_len_u]
        chunk_uncert=temps_1std[i-avg_len_l:i+avg_len_u]
        means[i] = np.mean(chunk)
        tot_uncert = np.var(chunk) + np.mean(chunk_uncert**2)
        ses[i] = np.sqrt(tot_uncert) / np.sqrt(len(chunk))
    return means, ses
