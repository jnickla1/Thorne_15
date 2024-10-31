import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
import os
def run_method(years, temperature, uncert, model_run, experiment_type):

    data_orig = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data_orig.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50])
    cur_path = os.path.dirname(os.path.realpath(__file__))
    comput_temps = pd.read_csv(cur_path+"/GAM_AR_Stephenson/gamAR0_fits_"+experiment_type+".csv", sep=',') #has header
    comput_uncert = pd.read_csv(cur_path+"/GAM_AR_Stephenson/gamAR0_se_fits_"+experiment_type+".csv", sep=',')
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)

    stidx = 29
    for endi in range(stidx, len(years)):

        means[endi] = comput_temps.values[endi-stidx, endi] -preind_base 
        ses[endi] = comput_uncert.values[endi-stidx, endi]
    #print(means[-1]+preind_base)
    return means, ses, comput_temps.values[endi-stidx, :] -preind_base, comput_uncert.values[endi-stidx, :]

