import numpy as np
import pandas as pd
import os
from . import EBMKF_Nicklas as ekf
#EBMKF_Nicklas as ekf

def run_method(years, temperature, uncert, model_run, experiment_type):
    #temperature has 2024, dont have all forcings for this yet - must update ekf.n_iters
    data_orig = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data_orig.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50])
    
    #ekf.temps = temperature[0:ekf.n_iters] +ekf.offset
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4
    
    #ekf.R_tvar=np.square(temps_1std[0:ekf.n_iters])

    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    means2 = empser.copy()
    ses2 = empser.copy()

    (trailfilter_xhat,P2s,S2s,trailfilter_xhatm)=ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)

    means[0:ekf.n_iters], ses[0:ekf.n_iters],means2[0:ekf.n_iters],ses2[0:ekf.n_iters] = ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)
    return means - ekf.offset-preind_base, np.sqrt(np.abs(ses)),means2 -ekf.offset-preind_base,np.sqrt(np.abs(ses2))
