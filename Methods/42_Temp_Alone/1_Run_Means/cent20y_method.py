import numpy as np

def run_method(years, temperature,uncert, model_run, experiment_type):
    avg_len_l=10
    avg_len_u=10 #actually 10 yrs below, that year, 9 yrs after
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    for i in range(avg_len_l, len(years) - avg_len_u):
        chunka=temperature[i-avg_len_l:i+avg_len_u]
        chunkb=temperature[i-avg_len_l+1:i+avg_len_u+1]
        chunka_uncert=temps_1std[i-avg_len_l:i+avg_len_u]
        chunkb_uncert=temps_1std[i-avg_len_l+1:i+avg_len_u+1]
        means[i] = np.mean([chunka,chunkb])
        tot_uncerta= np.var(chunka) + np.mean(chunka_uncert**2)
        tot_uncertb = np.var(chunkb) + np.mean(chunkb_uncert**2)
        ses[i] = (np.sqrt(tot_uncerta)+np.sqrt(tot_uncertb ))/ 2/np.sqrt(len(chunka))
    return means, ses

