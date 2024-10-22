#from Clarke & Richardson 2021, referencing Hausfather 2017
#https://www.science.org/doi/10.1126/sciadv.1601207
import warnings
import numpy as np
from .lowess1dt10wnc_method import lowess

import statsmodels.api as sm
from statsmodels.tsa.arima.model import ARIMA

# Assume residuals is your array of residuals from your trend fitting procedure

def estimate_arma_parameters(residuals, p=1, q=1):
    """
    Estimate ARMA(p, q) parameters using Maximum Likelihood Estimation (MLE).
    
    Parameters:
    - residuals: residual series (1D array-like)
    - p: AR order
    - q: MA order
    
    Returns:
    - phi: AR coefficients (array-like)
    - theta: MA coefficients (array-like)
    - model: fitted ARMA model object (from statsmodels)
    """
    # Fit ARMA model to the residuals
    with warnings.catch_warnings():
        # Ignore the specific UserWarning from statsmodels
        warnings.filterwarnings("ignore") 
                                #message="Non-stationary starting autoregressive parameters")
        
        # Fit ARMA model to the residuals
        model = ARIMA(residuals, order=(p, 0, q)).fit()

    # Extract AR and MA coefficients
    phi = model.arparams if p > 0 else np.array([])  # AR coefficients (phi)
    theta = model.maparams if q > 0 else np.array([])  # MA coefficients (theta)
    
    return phi, theta, model


def run_method(years, temperature, uncert, model_run, experiment_type):
    avg_len_l=25

    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    for i in range(avg_len_l, len(years)):
        chunktmps=temperature[0:i+1]
        chunkyrs=years[0:i+1]
        lfit,lsde, p_ests = lowess(chunkyrs, chunktmps, xwidth=20, degree=1, kernel='tricube', retP =True)

        residuals = chunktmps - lfit
        residuals_std = np.std(residuals)
        phi, theta, model = estimate_arma_parameters(residuals)
        rho1 = (phi+theta)*(1+phi*theta) / (1+2*phi*theta+theta*theta)
        phi_bc = phi + (1+ 4*(2*phi - rho1))/i
        rho1_bc = rho1 + (1+ 4*(2*phi - rho1))/i
        ne = i / (1+ 2*rho1_bc/(1-phi_bc))
        lsde_mod = lsde * np.sqrt(( ne - p_ests) / (i-p_ests))
        means[i] = lfit[-1]
        ses[i] = lsde_mod[-1]
    return means, ses, lfit, lsde_mod

