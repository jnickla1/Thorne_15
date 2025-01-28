import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import matplotlib.pyplot as plt


def run_method(years, temperature, uncert, model_run, experiment_type):
    avg_len_l=9
    avg_len_u=1 #actually 4 yrs below, that year, 0 yrs after
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    sesl = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    if (experiment_type!="historical"):
        return means, ses, means, sesl #forecasts not valid for future runs, return blanks
    
    cur_path = os.path.dirname(os.path.realpath(__file__))
    WMOoffset = 0.88 # for the WMO data 
    forec = pd.read_csv(cur_path+"/GlobalT_WMOLC-ADCPforecast_1991-2020.csv")
    nsamps = len(forec.columns)-2
    samp_cur = np.full((np.shape(years)[0],nsamps) ,np.nan)
        
    
    for i in range(110, len(years) - avg_len_u):
        chunk=temperature[i-avg_len_l:i+avg_len_u]
        chunk_uncert=temps_1std[i-avg_len_l:i+avg_len_u]
        forec_curyear = forec[forec['start_year']==(i+1850+1)].to_numpy()
        forec_samps = forec_curyear[0][2:] + WMOoffset
        #print(i)
        #print(forec_samps)
        means[i] = np.mean(chunk)/2 + np.nanmean(forec_samps)/2
        
        tot_uncert = np.var(chunk)/2 + np.mean(chunk_uncert**2)/2 + np.nanvar(forec_samps)/2
        ses[i] = np.sqrt(tot_uncert/ len(chunk))
        
        #don't need to increase uncertainty further - taking prior 10 years as having no uncertainty
        #tot_uncert0 =  np.nanvar(forec_samps)/2
        #ses[i] = np.sqrt(tot_uncert0/5)
        

    return means, ses, empser.copy(), empser.copy()

