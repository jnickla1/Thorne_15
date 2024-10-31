import numpy as np
from scipy.stats import linregress
import os
import pandas as pd
#import pdb

def run_method(years, temperature, uncert, model_run, experiment_type):

    means_c = np.full(np.shape(years),np.nan)
    ses_c = np.full(np.shape(years),np.nan)
    means_r = np.full(np.shape(years),np.nan)
    ses_r = np.full(np.shape(years),np.nan)
    
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    dir_path = os.path.dirname(os.path.realpath(__file__))
    lfcs_path = dir_path+"/lfca/ERSST_firstLFCA_s1940.npy"

    data_monthly = pd.read_csv(dir_path+"/../../../Common_Data/HadCRUT.5.0.2.0.analysis.summary_series.global.monthly.csv")
    temps_obs_monthly = data_monthly.loc[:,"Anomaly (deg C)"].to_numpy()
    preind_base = np.mean(temps_obs_monthly[0:50*12])
    
    st_idx =1940 - 1850
    lfcs = np.real(np.load(lfcs_path))

    for endi in range(st_idx, (2017-1850)):
        month_end_index = (endi-50)*12
        table_index = endi - st_idx
        temps_cropped = temps_obs_monthly[(1900-1850)*12 : (endi)*12]
        #breakpoint()
        slope, intercept, r_value, p_value, std_err = linregress(lfcs[table_index, 0:month_end_index], temps_cropped)
        #slope 0.3114061701545104
        predicted_temps = np.real(slope) * lfcs[table_index,0:month_end_index] + intercept
        residuals = temps_cropped - predicted_temps
        residual_std_dev = np.std(residuals)

        means_c[endi] = np.mean(predicted_temps[-13:-1])/2 -preind_base/2+ temperature[endi]/2 
        ses_c[endi] = residual_std_dev/np.sqrt(3)

    means_r[(1900-1850):(2016-1850)] = predicted_temps.reshape(-1, 12).mean(axis=1)/2 -preind_base/2 + temperature[(1900-1850):(2016-1850)]/2 
    ses_r[(1900-1850):(2016-1850)] = residual_std_dev/np.sqrt(3)

    return means_c ,ses_c, means_r, ses_r
