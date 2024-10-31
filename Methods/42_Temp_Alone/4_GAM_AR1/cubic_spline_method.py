#from Vissier 2018 https://doi.org/10.5194/cp-14-139-2018
#using the 0.28 AR directly 
import numpy as np
import patsy
from sklearn.linear_model import LinearRegression
from scipy import linalg

phi = 0.28  # AR(1) coefficient

# Function to simulate AR(1) process
def simulate_ar1_residuals(n, phi, residuals_std):
    residuals = np.zeros(n)
    residuals[0] = np.random.normal(0, residuals_std)
    for t in range(1, n):
        residuals[t] = phi * residuals[t-1] + np.random.normal(0, residuals_std)
    return residuals


import os
def run_method(years, temperature, uncert, model_run, experiment_type):


    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs

    # Number of MC samples
    n_bootstrap = 1000

    # Degree of polynomial
    degree = 4
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)



    dir_path = os.path.dirname(os.path.realpath(__file__))
    load_file = dir_path+"/AR_cubic_spline/cubic_spline_se_"+experiment_type+".npy"
    load_exist = False
    st_idx =30
    st_date = st_idx+1850

    if os.path.exists(load_file):
            # Load the weights from file
        ses_full = np.load(load_file)
        load_exist=True
        for endi in range(st_idx, len(years)):
            regX = np.arange(1850,1850+endi+1)
            n=len(regX)
            regY = temperature[regX-1850]
            # Fit the cubic spline basis
            spline_basis = patsy.dmatrix("bs(regX, df=7, degree=3, include_intercept=False)", {"regX": regX}, return_type='dataframe')
     
        # Fit a linear model to the spline-transformed features
            model = LinearRegression()
            model.fit(spline_basis, regY)
    
        # Predict the fit values
            y_pred = model.predict(spline_basis)
            means[endi] = y_pred[-1]
            ses[endi] = ses_full[endi-st_idx, endi]
    else:
        ses_full = np.empty((len(years)-st_idx, len(years)))
        
        for endi in range(st_idx, len(years)):
            regX = np.arange(1850,1850+endi+1)
            n=len(regX)
            regY = temperature[regX-1850]
        # Fit the cubic spline basis
            spline_basis = patsy.dmatrix("bs(regX, df=7, degree=3, include_intercept=False)", {"regX": regX}, return_type='dataframe')
     
        # Fit a linear model to the spline-transformed features
            model = LinearRegression()
            model.fit(spline_basis, regY)
    
        # Predict the fit values
            y_pred = model.predict(spline_basis)
            means[endi] = y_pred[-1]
            residuals = regY - y_pred
            residuals_std = np.std(residuals)

            boot_fits = np.zeros((n_bootstrap, n))  # Store AR1 montecarlo runs
            print(np.shape(boot_fits))        
            for sim in range(n_bootstrap):
                 # Simulate AR(1) residuals
                simulated_residuals = simulate_ar1_residuals(n, phi, residuals_std)
        
                # Step 4: Generate new surrogate series by adding simulated residuals to original spline fit
                surrogate_Y = y_pred + simulated_residuals
        
                # Step 5: Fit new spline to the surrogate data
                model.fit(spline_basis, surrogate_Y)
                y_surrogate_pred = model.predict(spline_basis)

                # Store the last predicted value for the surrogate series
                boot_fits[sim, :] = y_surrogate_pred

            ses_full[endi-st_idx,0:(endi+1)]= (np.percentile(boot_fits, 84, axis=0)  - np.percentile(boot_fits, 16, axis=0))/2 #np.std(boot_fits,axis=0)
            ses[endi] = ses_full[endi-st_idx, endi ]

    if not(load_exist):
        np.save(load_file, ses_full)

    return means, ses, y_pred,ses_full[endi-st_idx,:]
