import numpy as np
from numpy.polynomial import Polynomial
import os
def run_method(years, temperature, uncert, model_run, experiment_type):


    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs

    # Number of bootstrap samples
    n_bootstrap = 1000

    # Degree of polynomial
    degree = 4
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)



    dir_path = os.path.dirname(os.path.realpath(__file__))
    load_file = dir_path+"/boot_quartic/quartic_se_"+experiment_type+str(model_run)+".npy"
    load_exist = False
    st_idx =100
    st_date = st_idx+1850

    if os.path.exists(load_file):
            # Load the weights from file
        ses_full = np.load(load_file)
        load_exist=True
        for endi in range(st_idx, len(years)):
            regX = np.arange(1850,1850+endi+1)
            n=len(regX)
            regY = temperature[regX-1850]
            # Fit the original quartic polynomial
            p_orig = Polynomial.fit(regX, regY, degree)
            y_pred = p_orig(regX)
            means[endi] = y_pred[-1]
            ses[endi] = ses_full[endi-st_idx, endi]
    else:
        ses_full = np.empty((len(years)-st_idx, len(years)))
        
        for endi in range(st_idx, len(years)):
            regX = np.arange(1850,1850+endi+1)
            n=len(regX)
            regY = temperature[regX-1850]
            # Fit the original quartic polynomial
            p_orig = Polynomial.fit(regX, regY, degree)
            y_pred = p_orig(regX)
            means[endi] = y_pred[-1]

            uregY = temps_1std[regX-1850]
            boot_fits = np.zeros((n_bootstrap, n))  # Store bootstrap coefficients
            print(np.shape(boot_fits))        
            for i in range(n_bootstrap):
                # Resample with replacement
                indices = np.random.choice(n, size=n, replace=True)
                x_resample = regX[indices]
                y_resample = np.random.normal(regY[indices],uregY[indices]) #accounts for HadCRUT5 uncertainty when resampling
                                # Fit a quartic polynomial to the resampled data
                p_boots = Polynomial.fit(x_resample, y_resample, degree)
                boot_fits[i, :] = p_boots(regX)

            ses_full[endi-st_idx,0:(endi+1)]= np.std(boot_fits,axis=0)
            ses[endi] = ses_full[endi-st_idx, endi ]

    if not(load_exist):
        np.save(load_file, ses_full)

    return means, ses, y_pred,ses_full[endi-st_idx,:]
