import numpy as np
from .lowess1dt10wnc_method import lowess


def run_method(years, temperature, uncert, model_run, experiment_type):
    avg_len_l=3

    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    for i in range(avg_len_l, len(years)):
        chunktmps=temperature[0:i+1]
        chunkyrs=years[0:i+1]
        lfit,lsde = lowess(chunkyrs, chunktmps, xwidth=36, degree=1, kernel='tricube')
        means[i] = lfit[-1]
        ses[i] = lsde[-1]
    return means, ses, lfit, lsde

