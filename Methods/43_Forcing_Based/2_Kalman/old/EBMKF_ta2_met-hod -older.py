import numpy as np
import pandas as pd
from . import EBMKF_Nicklas2 as ekf
#EBMKF_Nicklas as ekf

def run_method(years, temperature, uncert, model_run, experiment_type):
    #temperature has 2024, dont have all forcings for this yet - must update ekf.n_iters
    data_orig = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data_orig.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50])
    
    #ekf.temps = temperature[0:ekf.n_iters] +ekf.offset
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4
    #ekf.Q=(ekf.covM/20) - does worse!
    #ekf.R_tvar=np.square(temps_1std[0:ekf.n_iters])
    #ekf.fdbkA = 0.4
    #ekf.precompute_coeffs(False)
    
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    means2 = empser.copy()
    ses2 = empser.copy()

    wt_opt_depths = 1/(ekf.opt_depth+9.7279)
    N = 30
    nwt_opt_depths=np.empty(len(ekf.opt_depth)); nwt_opt_depths[:]=ekf.involcavg
    cN=int(np.ceil(N/2))
    fN=int(np.floor(N/2))
    for i in range((fN),(len(nwt_opt_depths)-1)):
        lasta=i+cN;firsta=i-fN
        nwt_opt_depths[i] = (np.sum(wt_opt_depths[(firsta):i+1])+ ekf.involcavg*(cN-1))/N
        #computing half-average - future is assumed to be the average 
    nopt_depths=(1/nwt_opt_depths-9.7279)
    ekf.opt_depth=nopt_depths
    #(trailfilter_xhat,P2s,S2s,trailfilter_xhatm)=ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)

    means[0:ekf.n_iters], ses[0:ekf.n_iters],means2[0:ekf.n_iters],ses2[0:ekf.n_iters] = ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)
    return means - ekf.offset-preind_base, np.sqrt(np.abs(ses)),means2 -ekf.offset-preind_base,np.sqrt(np.abs(ses2))
