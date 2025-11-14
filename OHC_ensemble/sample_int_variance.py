import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from record_ens_stems import *
import os

# Set random seed for reproducibility
np.random.seed(42)

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading ensemble arrays...")
filepath = os.path.dirname(os.path.abspath(__file__))
ohca_change = np.load(filepath+'/ohca_change_ensemble.npy')
ohca_uncertainty = np.load(filepath+'/ohca_uncertainty_ensemble.npy')
coverage_string = np.load(filepath+'/coverage_string_ensemble.npy', allow_pickle=True)


N_YEARS = len(YEARS)

start_idx = np.where(YEARS == START_YEAR)[0][0]

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def calculate_ar1_coefficient(deviances):
    """Calculate AR(1) coefficient from deviances"""
    if len(deviances) < 2 or np.all(np.isnan(deviances)):
        return 0.0
    
    valid = ~np.isnan(deviances)
    if np.sum(valid) < 2:
        return 0.0
    
    y = deviances[valid][1:]
    x = deviances[valid][:-1]
    
    if len(x) == 0 or np.std(x) == 0:
        return 0.0
    
    ar1 = np.corrcoef(x, y)[0, 1]
    return ar1 if not np.isnan(ar1) else 0.0


def skewnorm_params_from_errors(pos_err, neg_err):
    """
    Convert positive and negative standard errors to skew-normal parameters.
    Returns (xi, omega, alpha) for scipy.stats.skewnorm
    """
    # Calculate moments
    mu = np.sqrt(2/np.pi) * (pos_err - neg_err)
    sigma_sq = 0.5 * (pos_err**2 + neg_err**2)
    sigma = np.sqrt(sigma_sq)
    
    if sigma < 1e-10:
        return 0.0, 1.0, 0.0  # Degenerate case
    
    gamma = (2/np.sqrt(2*np.pi)) * (pos_err**3 - neg_err**3) / sigma**3
    
    # Solve for delta from skewness equation
    def skewness_equation(delta):
        if abs(delta) >= 1:
            return 1e10
        term = delta * np.sqrt(2/np.pi) / np.sqrt(1 - 2*delta**2/np.pi)
        return ((4 - np.pi)/2) * term**3 - gamma
    
    # Initial guess
    delta_init = np.sign(gamma) * 0.5
    try:
        delta = fsolve(skewness_equation, delta_init)[0]
        delta = np.clip(delta, -0.99, 0.99)
    except:
        delta = 0.0
    
    # Calculate alpha, omega, xi
    alpha = delta / np.sqrt(1 - delta**2) if abs(delta) < 0.99 else 0.0
    omega = sigma / np.sqrt(1 - 2*delta**2/np.pi) if abs(delta) < 0.99 else sigma
    xi = mu - omega * delta * np.sqrt(2/np.pi)
    
    return xi, omega, alpha


def transform_to_skewnorm(z, pos_err, neg_err):
    """Transform standard normal z-score to skewed normal sample"""
    xi, omega, alpha = skewnorm_params_from_errors(pos_err, neg_err)
    
    # CDF of standard normal at z
    p = stats.norm.cdf(z)
    
    # Inverse CDF of skew normal
    return stats.skewnorm.ppf(p, alpha, loc=xi, scale=omega)


def get_average_ohca():
    ohca_change_mean = ohca_change.mean(axis=0)
    ohca_mean2 = np.nancumsum(ohca_change_mean)
    ohca_uncert_mean = np.abs(ohca_uncertainty).mean(axis=0)
    ohca_uncert2 = np.minimum(ohca_uncert_mean,ohca_uncert_mean[0]) + 5 #clamp to max uncertainty of Zanna
    return ohca_mean2, ohca_uncert2


if __name__ == "__main__":
    # ============================================================================
    # CALCULATE AR1 COEFFICIENTS
    # ============================================================================

    print("Calculating AR1 coefficients...")
    ar1_array = np.zeros((N_ENSEMBLE, N_YEARS))

    WINDOW_SIZE = 20

    for row_idx in range(N_ENSEMBLE):
        for year_idx in range(N_YEARS):
            # Define window: -20 to +20 years
            start = max(0, year_idx - WINDOW_SIZE)
            end = min(N_YEARS, year_idx + WINDOW_SIZE + 1)
            
            window_data = ohca_change[row_idx, start:end]
            
            # Remove mean to get deviances
            valid_mask = ~np.isnan(window_data)
            if np.sum(valid_mask) > 0:
                mean_change = np.nanmean(window_data)
                deviances = window_data - mean_change
                
                # Calculate AR1
                ar1_array[row_idx, year_idx] = calculate_ar1_coefficient(deviances)

    print(f"AR1 range: [{np.nanmin(ar1_array):.3f}, {np.nanmax(ar1_array):.3f}]")
    print("Cropping to only positive values")
    ar1_array[ar1_array<0]=0

    # ============================================================================
    # GENERATE AR1-EVOLVED Z-SCORES
    # ============================================================================

    print("Generating AR1-evolved z-scores...")
    z_scores = np.zeros((N_ENSEMBLE, N_YEARS))

    internal_var_MAX_loc = 40


    for row_idx in range(N_ENSEMBLE):
        # Initialize first year
        z_scores[row_idx, 0] = np.random.normal(0, 1)
        
        # Evolve through time
        for year_idx in range(1, N_YEARS):
            ar1 = ar1_array[row_idx, year_idx]
            ar1 = np.clip(ar1, -0.99, 0.99)  # Ensure stationarity
            
            white_noise = np.random.normal(0, 1)
            
            z_scores[row_idx, year_idx] = (ar1 * z_scores[row_idx, year_idx-1] + 
                                           np.sqrt(1 - ar1**2) * white_noise)
            internal_var_MAX = np.random.normal(loc=internal_var_MAX_loc, scale=internal_var_MAX_loc/4)
            #compute change in total
            ohca_uncer_now = ohca_uncertainty[row_idx, year_idx].real
            ohca_change_now = (z_scores[row_idx, year_idx] - z_scores[row_idx, year_idx-1])* ohca_uncer_now
            if abs(ohca_change_now)>internal_var_MAX:
                #slightly imperfect calculation, force some z-scores to evolve more slowly
                #np.sqrt(1 - ar1**2) * white_noise) * ohca_uncertainty[row_idx, year_idx].real / internal_var_MAX = 1
                ar1bsqrt_fact = internal_var_MAX / ohca_uncer_now
                #make z_score even more in direction it's already going / could do np.sign(white_noise)
                z_scores[row_idx, year_idx] = (ar1 * z_scores[row_idx, year_idx-1] + ar1bsqrt_fact*np.sign(z_scores[row_idx, year_idx-1]))
                
            

    print(f"Z-score variance by year: {np.var(z_scores, axis=0)[:5]} (should be ~1)")

    # ============================================================================
    # SCALE BY UNCERTAINTIES AND HANDLE ASYMMETRIC ERRORS
    # ============================================================================

    # Compute cumulative OHCA from changes (anchored at 0 in START_YEAR)
    ohca_cumulative = np.nancumsum(ohca_change, axis=1)
    ohca_cumulative = ohca_cumulative - ohca_cumulative[:, start_idx:start_idx+1]


    print("Scaling by uncertainties...")
    internal_uncertainty = np.zeros((N_ENSEMBLE, N_YEARS))

    for row_idx in range(N_ENSEMBLE):
        for year_idx in range(N_YEARS):
            unc = ohca_uncertainty[row_idx, year_idx]
            z = z_scores[row_idx, year_idx]
            
            if np.isnan(unc) or (isinstance(unc, complex) and np.isnan(unc.real)):
                internal_uncertainty[row_idx, year_idx] = 0
                continue
            
            # Check if asymmetric (complex with imaginary part)
            if isinstance(unc, complex) and abs(unc.imag) > 1e-10:
                # Decode to pos/neg errors
                pos_err , neg_err = decode_asymmetric_uncertainty(unc)
                # Transform z to skewed normal
                internal_uncertainty[row_idx, year_idx] = transform_to_skewnorm(z, pos_err, neg_err)
            else:
                # Symmetric uncertainty
                std = unc.real if isinstance(unc, complex) else unc
                #existing uncertainty from the existing stem ensemble
                block_idx = row_idx // N_COPIES_PER_RECORD
                existing_std = np.std(ohca_cumulative[(block_idx*N_COPIES_PER_RECORD):((block_idx+1)*N_COPIES_PER_RECORD), year_idx])
                new_std = np.sqrt(max(0,std**2 - existing_std**2)) #keep it above 0
                internal_uncertainty[row_idx, year_idx] = z * new_std
                

    # ============================================================================
    # ADD TO CUMULATIVE OHCA
    # ============================================================================

    print("Computing sampled OHCA timeseries...")


    # Add internal uncertainty
    ohca_sampled = ohca_cumulative + internal_uncertainty
    ohca_sampled = ohca_sampled - ohca_sampled[:, start_idx:start_idx+1]

    # Save output
    np.save('ohca_sampled_ensemble.npy', ohca_sampled)
    print("Saved ohca_sampled_ensemble.npy")

    # ============================================================================
    # PLOT
    # ============================================================================

    print("Plotting sampled ensemble...")


    plot_ensemble_stems(ohca_sampled, coverage_string, YEARS, 
                       'Sampled OHCA Ensemble with Internal Uncertainty',
                       'ensemble_sampled.png',changes=False)

    print("\nDone! Check ensemble_sampled.png")
