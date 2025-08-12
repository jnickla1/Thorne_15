import numpy as np
import pandas as pd
import os

def run_method(years, temperature,uncert, model_run, experiment_type):
    avg_len_l=5
    avg_len_u=5 #actually 5 yrs below, that year, 4 yrs after
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    if experiment_type == 'historical':
        cur_path = os.path.dirname(os.path.realpath(__file__))
        hadcrutUncertFile = pd.read_csv(cur_path+"/../../../Common_Data/HadCRUT.5.0.2.0_10_yr_run_mean.csv") 
        hadcrutUncert = hadcrutUncertFile["Total_uncert_(1sigma)"]

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

