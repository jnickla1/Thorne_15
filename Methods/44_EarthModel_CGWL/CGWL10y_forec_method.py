import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import matplotlib.pyplot as plt


def run_method(years, temperature, uncert, model_run, experiment_type):
    
    
    avg_len_l=9
    avg_len_u=1 #actually 4 yrs below, that year, 0 yrs after
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    sesl = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    if (experiment_type!="historical"):
        return means, ses, means, sesl #forecasts not valude for future runds, return blanks
    
    cur_path = os.path.dirname(os.path.realpath(__file__))
    WMOoffset = 0.88 # for the WMO data 
    forec = pd.read_csv(cur_path+"/GlobalT_WMOLC-ADCPforecast_1991-2020.csv")
    nsamps = len(forec.columns)-2
    samp_cur = np.full((np.shape(years)[0],nsamps) ,np.nan)
        
    
    for i in range(110, len(years) - avg_len_u):
        chunk=temperature[i-avg_len_l:i+avg_len_u]
        chunk_uncert=temps_1std[i-avg_len_l:i+avg_len_u]
        forec_curyear = forec[forec['start_year']==(i+1850+1)].to_numpy()
        forec_samps = forec_curyear[0][2:] + WMOoffset
        #print(i)
        #print(forec_samps)
        means[i] = np.mean(chunk)/2 + np.nanmean(forec_samps)/2
        
        #tot_uncert = np.var(chunk)/2 + np.mean(chunk_uncert**2)/2 
        #ses[i] = np.sqrt(tot_uncert/ len(chunk) + np.nanvar(forec_samps)/2/sum(~np.isnan(forec_samps)))
        #don't need to increase uncertainty further - taking prior 10 years as having no uncertainty
        tot_uncert0 =  np.nanvar(forec_samps)
        ses[i] = np.sqrt(tot_uncert0)/2
        
        samp_cur[i,:] = ( np.mean(chunk)/2 + forec_samps/2 )#+
            #np.random.normal(loc=0, scale=np.sqrt(np.var(chunk)/2/len(chunk)), size=nsamps)+
            #np.random.normal(loc=0, scale=np.sqrt(np.mean(chunk_uncert**2)/2/len(chunk)), size=nsamps) )


    def empirical_mean(year_idx,k):
        if (k==0):
            return means[year_idx-1850]
        else:
            return empser

    def empirical_se(year_idx,k):
        if (k==0):
            return ses[year_idx-1850]
        else:
            return empser

    #73 samples at start
    def empirical_pvalue(year_idx, point,k,two_sided=True,tail_threshold=1/20, num_closest_samples=10):
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
                dist0 = samp_cur[year_idx[i]-1850, :]  # Distribution for the current year
                mask = np.isnan(dist0)
                dist = dist0[~mask]
            else:
                return empser
            if(sum(~np.isnan(dist))>0 and ~np.isnan(point[i])):
                epdf = stats.gaussian_kde(dist)
                empirical_ll[i] = epdf.logpdf(point[i])
##            if False:
##                plt.figure()
##                xfine = np.linspace(0.5,1.5,100)
##                plt.hist(dist,density=True)
##                epdf = stats.gaussian_kde(dist)
##                plt.plot(xfine,epdf.pdf(xfine))
##                plt.plot(xfine,stats.norm.pdf(xfine, loc=means[i], scale=ses[i]))
##                plt.plot(xfine,stats.norm.pdf(xfine, loc=means[i], scale=sesl[i]))
##                plt.show()
                
        return empirical_ll    

    return {
        'mean': empirical_mean,
        'se': empirical_se,
        'pvalue': empirical_pvalue,
        'log_likelihood': empirical_log_likelihood,

    }

    #return means, ses, empser.copy(), empser.copy()

