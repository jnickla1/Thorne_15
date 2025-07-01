import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
import os
from scipy import stats

def run_method(years, temperature, uncert, model_run, experiment_type):

    empser  = np.full(len(years),np.nan)

    
    #data_orig = pd.read_csv("./Common_Data/HadCRUT5.csv")
    #temps_obs = data_orig.loc[:,"Anomaly"].to_numpy()
    #preind_base = np.mean(temps_obs[0:50])
    
    cur_path = os.path.dirname(os.path.realpath(__file__))

    obt_array0 = np.load(cur_path+"/temperature-attribution/output/temp_all.npy")
    Nres = np.shape(obt_array0)[1]

#starts in 1750 so crop to 100th index
    samp_cur = obt_array0[100:,:] 
    
    def empirical_mean(year_idx,k):
        if (k==0):
            return np.mean(samp_cur[year_idx-1850, :], axis = 1)
        else:
            return empser.copy()

    def empirical_se(year_idx,k):
        if (k==0):
            return np.std(samp_cur[year_idx-1850, :], axis = 1)
        else:
            return empser.copy()

    def empirical_pvalue(year_idx, point,k,two_sided=True,tail_threshold=1/Nres*4, num_closest_samples=int(Nres/7)):
        if(k!=0):
            return np.full(np.shape(year_idx),np.nan)
        year_idx = np.atleast_1d(year_idx)
        point = np.atleast_1d(point)
        # Initialize an array to store empirical p-values
        empirical_p = np.full(np.shape(point), np.nan)
        for i in range(len(year_idx)):
            if ~np.isnan(point[i]):
                if (k==0):
                    dist0 = samp_cur[year_idx[i]-1850, :]  # Distribution for the current year
                else:
                    dist0 = samp_ret[year_idx[i]-1850, :]  # Distribution for the current year
                dist = dist0[~np.isnan(dist0)]
                cdist = np.nanmean(dist)
                if(len(dist)==0):
                    continue
                
                empirical_p_count2 = 2* min( np.sum((dist - point[i]) >= 0),  np.sum((- dist + point[i]) >= 0))/ len(dist) #count number of dist samples farther away than point
                empirical_p_count = np.sum((dist - point[i]) >= 0)/ len(dist) #point[i] is greater than how many within dist?
                    
                if empirical_p_count2 > tail_threshold:
                    if two_sided:
                        empirical_p[i] = empirical_p_count2
                    else:
                        empirical_p[i] = empirical_p_count
                else:
                # Apply exponential tail fitting for very small p-values
                # Estimate parameters for the exponential distribution fitted to the tail
                    if cdist>point[i]:
                        sorted_dist = np.sort((dist - point[i]))  # Smallest to largest, large means closer to the center
                    else:
                        sorted_dist = np.sort((- dist + point[i]))
                    
                    closest_samples = sorted_dist[:num_closest_samples]  # Take the n closest, curve now looks like this  / or "r"
                    # Estimate the shift (smallest absolute difference)
                    shift = closest_samples[-1]  # This is the max absolute distance or the farthest to the right, will be positive (samp originally to the inside of point[i]
                    # Estimate the rate parameter using the remaining closest samples (subtract shift)
                    adjusted_samples = shift - closest_samples  # Subtract shift to get differences, now curve looks like |\
                    lambda_param = 1 / np.mean(adjusted_samples)  # Rate parameter (1 / mean of adjusted differences)
                # Use the exponential distribution to estimate the tail p-value
                    tail_p_value = stats.expon.sf(shift, scale=1/lambda_param)
                
                # Combine empirical and tail-based p-value
                    if two_sided:
                        empirical_p[i] = 2*tail_p_value
                    else:
                        if cdist>point[i]: #point evaluating at is far to the left of the distribution, so whe have passed it
                            empirical_p[i] = 1-tail_p_value
                        else:
                            empirical_p[i] = tail_p_value
        return empirical_p

    def empirical_log_likelihood(year_idx, point,k):
        if(k!=0):
            return np.full(np.shape(year_idx),np.nan)
        year_idx = np.atleast_1d(year_idx)
        point = np.atleast_1d(point)
        # Initialize an array to store empirical p-values
        empirical_ll = np.full(np.shape(point), np.nan)
        for i in range(len(year_idx)):
            if (k==0):
                dist = samp_cur[year_idx[i]-1850, :]  # Distribution for the current year
            else:
                return empser.copy()
            if(sum(np.isnan(dist))==0 and ~np.isnan(point[i])):
                epdf = stats.gaussian_kde(dist)
                empirical_ll[i] = epdf.logpdf(point[i])
        return empirical_ll

    return {
        'mean': empirical_mean,
        'se': empirical_se,
        'pvalue': empirical_pvalue,
        'log_likelihood': empirical_log_likelihood,

    }
