import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
import os
from scipy import stats

def run_method(years, temperature, uncert, model_run, experiment_type):
    Nres = 1000
    empser  = np.full((len(years),Nres),np.nan)
    samp_cur=empser.copy()


    
    cur_path = os.path.dirname(os.path.realpath(__file__))

    def obtain_file(year):
        df = pd.read_csv(cur_path+"/KCC/kcc_notebook/kcc_"+experiment_type+"_all_"+str(year)+".csv", sep=',') #has header
        dfnat = pd.read_csv(cur_path+"/KCC/kcc_notebook/kcc_"+experiment_type+"_nat_"+str(year)+".csv", sep=',') #has header
        human_samp =  df.drop(df.columns[0], axis=1).to_numpy() - dfnat.drop(dfnat.columns[0], axis=1).to_numpy() 
        return human_samp #remove first column


    stidx = 100
    for endi in range(stidx, len(years)):
        obt_array = obtain_file(endi+1850)
        #print(obt_array)
        samp_cur[endi,:] = obt_array[endi,:]
        
    samp_ret=obtain_file(len(years)-1+1850)
#load all samples into this module
    
    def empirical_mean(year_idx,k):
        if (k==0):
            return np.mean(samp_cur[year_idx-1850, :], axis = 1)
        else:
            return np.mean(samp_ret[year_idx-1850, :], axis = 1)

    def empirical_se(year_idx,k):
        if (k==0):
            return np.std(samp_cur[year_idx-1850, :], axis = 1)
        else:
            return np.std(samp_ret[year_idx-1850, :], axis = 1)

    def empirical_pvalue(year_idx, point,k,tail_threshold=5e-3, num_closest_samples=100):
        year_idx = np.atleast_1d(year_idx)
        point = np.atleast_1d(point)
        # Initialize an array to store empirical p-values
        empirical_p = np.full(np.shape(point), np.NaN)
        for i in range(len(year_idx)):
            if ~np.isnan(point[i]):
                if (k==0):
                    dist = samp_cur[year_idx[i]-1850, :]  # Distribution for the current year
                else:
                    dist = samp_ret[year_idx[i]-1850, :]  # Distribution for the current year
                empirical_p_count = 2* min( np.sum((dist - point[i]) >= 0),  np.sum((- dist + point[i]) >= 0))/ len(dist) #count number of dist samples farther away than point
                if empirical_p_count > tail_threshold:
                    empirical_p[i] = empirical_p_count
                else:
                # Apply exponential tail fitting for very small p-values
                # Estimate parameters for the exponential distribution fitted to the tail
                    sorted_dist = np.sort(np.abs(dist - point[i]))  # Sort absolute differences
                    closest_samples = sorted_dist[:num_closest_samples]  # Take the 100 closest, curve now looks like this  / or "r"
                    # Estimate the shift (smallest absolute difference)
                    shift = closest_samples[-1]  # This is the max absolute distance
                    # Estimate the rate parameter using the remaining closest samples (subtract shift)
                    adjusted_samples = shift - closest_samples  # Subtract shift to get differences, now curve looks like |\
                    lambda_param = 1 / np.mean(adjusted_samples)  # Rate parameter (1 / mean of adjusted differences)
                # Use the exponential distribution to estimate the tail p-value
                    tail_p_value = stats.expon.sf(np.abs(point[i] - shift), scale=1/lambda_param)
                
                # Combine empirical and tail-based p-value
                    empirical_p[i] = tail_p_value
        return empirical_p

    def empirical_log_likelihood(year_idx, point,k):
        year_idx = np.atleast_1d(year_idx)
        point = np.atleast_1d(point)
        # Initialize an array to store empirical p-values
        empirical_ll = np.full(np.shape(point), np.NaN)
        for i in range(len(year_idx)):
            if ~np.isnan(point[i]):
                if (k==0):
                    dist = samp_cur[year_idx[i]-1850, :]  # Distribution for the current year
                else:
                    dist = samp_ret[year_idx[i]-1850, :]  # Distribution for the current year                        
                empirical_ll[i] = -np.log(np.std(dist)) - 0.5 * ((point[i] - np.mean(dist)) / np.std(dist)) ** 2
        return empirical_ll

    return {
        'mean': empirical_mean,
        'se': empirical_se,
        'pvalue': empirical_pvalue,
        'log_likelihood': empirical_log_likelihood,

    }
