import numpy as np
import pandas as pd

from scipy.stats import norm

from EBMKF_Nicklas4 import Teq1850

def save_analytic_percentiles_from_mean_se(output_csv, year_range, data_dict, percentiles=[5, 17, 83, 95, 50]):
    """
    For each variable in data_dict, compute specified percentiles from N(mean, se) using z-scores.
    
    Parameters:
        output_csv (str): Filename to write output CSV.
        year_range (array-like): List of years, e.g., np.arange(1850, 2025).
        data_dict (dict): Keys are variable names, values are tuples (mean_vector, se_vector).
        percentiles (list): List of percentiles to compute (default: [5, 17, 50, 83, 95]).
    """
    df = pd.DataFrame({'year': year_range})
    
    z_scores = {p: norm.ppf(p / 100) for p in percentiles}

    for varname, (mean_vec, se_vec) in data_dict.items():
        for p in percentiles:
            z = z_scores[p]
            df[f'{varname}_p{p}'] = mean_vec + z * se_vec - Teq1850

    df.to_csv(output_csv, index=False)


years = np.arange(1850, 2024)
n = len(years)


files = [
    'Nicklas_GMST_timeseriesBCC-CSM2-MR.npz',
    'Nicklas_GMST_timeseriesCESM2.npz',
    'Nicklas_GMST_timeseriesESM1-2-LR.npz'
]

# Load all files
loaded = [np.load(f) for f in files]

# Define variables
var_names = ['aer', 'ghg', 'ohf', 'anthro', 'nat', 'tot']

data = {}

for var in var_names:
    # Collect _vals and _se from all files
    vals_list = [ld[f'{var}_vals'] for ld in loaded]
    se_list   = [ld[f'{var}_se']   for ld in loaded]

    if var=='nat':
        vals_list = [ld[f'{var}_vals'] for ld in loaded[1:]]
        se_list   = [ld[f'{var}_se']   for ld in loaded[1:]]
        #nat BCC-CSM2-MR looks very strange
    # Convert to arrays for math
    vals_array = np.array(vals_list)  # shape (3, n)
    se_array   = np.array(se_list)    # shape (3, n)
    
    # Average of the vals
    mean_vals = np.mean(vals_array, axis=0)
    
    # Internal variance: mean of squared SEs
    internal_var = np.mean(se_array**2, axis=0)
    
    # External variance: variance across vals
    external_var = np.var(vals_array, axis=0, ddof=0)
    
    # Total standard error
    total_se = np.sqrt(internal_var + external_var)
    
    # Add to dictionary
    key = 'ant' if var == 'anthro' else var  # rename 'anthro' -> 'ant'
    data[key] = (mean_vals, total_se)
    
save_analytic_percentiles_from_mean_se("Nicklas_GMST_timeseriesCombined.csv", years, data)
