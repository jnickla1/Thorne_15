#from eq 5,6,&11 of Clarke & Richardson 2021
#https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020EA001082

import numpy as np
from .lowess1dt10wnc_method import lowess

from sklearn.linear_model import LinearRegression

def estimate_phi(residuals):
    """
    Estimate the AR(1) coefficient phi using OLS regression on the lagged residuals.
    
    Parameters:
    - residuals: numpy array of residuals from the original fit
    
    Returns:
    - phi: the estimated AR(1) coefficient
    """
    # Create lagged residuals
    r_t = residuals[1:]    # r_t (time t)
    r_t_minus_1 = residuals[:-1]  # r_(t-1) (time t-1)

    # Fit linear regression r_t ~ r_(t-1)
    model = LinearRegression()
    model.fit(r_t_minus_1.reshape(-1, 1), r_t)

    # The slope of the regression line is the estimate for phi
    phi = model.coef_[0]
    return phi



def run_method(years, temperature, uncert, model_run, experiment_type):
    avg_len_l=3

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
        phi_estimated = estimate_phi(residuals)
        lsde_mod = lsde * np.sqrt((i / (1+ phi_estimated) * (1-phi_estimated) - p_ests) / (i-p_ests))
        means[i] = lfit[-1]
        ses[i] = lsde_mod[-1]
    return means, ses, lfit, lsde_mod

