import numpy as np
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
import netCDF4


np.random.seed(402)

def run_method(years, temperature, uncert, model_run, experiment_type):
    avg_len_l=9
    avg_len_u=1 #actually 4 yrs below, that year, 0 yrs after
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    sesl = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    cur_path = os.path.dirname(os.path.realpath(__file__))
    WMOoffset = 0.88 # for the WMO data 
    exp_attr = experiment_type.split("_")
    
    if (experiment_type=="historical" or exp_attr[2].lower()=="ssp245" or exp_attr[2].lower()=="rcp45"):
        filein = netCDF4.Dataset(cur_path+"/tasAnom_rcp45_land-prob_global_glb_sample_b8100_1y_ann_18591201-20991130.nc",'r')
        comput_temps = np.array(filein.variables['tasAnom'])
    elif (exp_attr[2].lower()=="ssp126"):
        filein = netCDF4.Dataset(cur_path+"/tasAnom_rcp26_land-prob_global_glb_sample_b8100_1y_ann_18591201-20991130.nc",'r')
        comput_temps = np.array(filein.variables['tasAnom'])
    elif (exp_attr[2].lower()=="ssp370"):
        filein60 = netCDF4.Dataset(cur_path+"/tasAnom_rcp60_land-prob_global_glb_sample_b8100_1y_ann_18591201-20991130.nc",'r')
        comput_temps60 = np.array(filein60.variables['tasAnom'])
        filein85 = netCDF4.Dataset(cur_path+"/tasAnom_rcp85_land-prob_global_glb_sample_b8100_1y_ann_18591201-20991130.nc",'r')
        comput_temps85 = np.array(filein85.variables['tasAnom'])
        comput_temps = np.concatenate((comput_temps60.T,comput_temps85.T)).T
    #(240, 3000), 1860 is first year
    nsamps = np.shape(comput_temps)[1]
    cutoff_n = int(nsamps/30)
    comput_temps_baselines = np.mean(comput_temps[0:40,:],axis=0) #1860-1900 baselines
    comput_temps_align = comput_temps - np.tile(comput_temps_baselines, (np.shape(comput_temps)[0], 1))
    
    samp_cur = np.full((np.shape(years)[0],cutoff_n) ,np.nan)
        
    offset_yrs= 1860-years[0]
    last_i = min(len(years) - avg_len_u, np.shape(comput_temps)[0] - 11 + offset_yrs )
    for i in range(avg_len_l+10, last_i):
        chunk=temperature[i-avg_len_l:i+avg_len_u]
        chunk_uncert=temps_1std[i-avg_len_l:i+avg_len_u]
        chunk_avg = np.mean(chunk)
        hindc_samps = np.mean(comput_temps_align[i-avg_len_l-offset_yrs:i+avg_len_u-offset_yrs,:], axis = 0)
        forec_samps0 = np.mean(comput_temps_align[(i+1-offset_yrs):(i+11-offset_yrs),:], axis = 0)

        sort_indices = np.argsort(np.abs(hindc_samps-chunk_avg))
        forec_samps = forec_samps0[sort_indices[:cutoff_n]] #taking the 1/30th of samples with past 10yr mean closest to observed, so 100
        
        #print(i)
        #print(forec_samps)
        means[i] = np.mean(chunk)/2 + np.nanmean(forec_samps)/2
        
        #tot_uncert = np.var(chunk) + np.mean(chunk_uncert**2) #this is going into a standard error
        #don't need to increase uncertainty further - taking prior 10 years as having no uncertainty
        tot_uncert0 =  np.nanvar(forec_samps)
        ses[i] = np.sqrt(tot_uncert0 )/2 #+tot_uncert/ len(chunk) 
        samp_cur[i,:] = ( np.mean(chunk)/2 + forec_samps/2 ) #integer added to vector
            #np.random.normal(loc=0, scale=np.sqrt(np.var(chunk)/2/len(chunk)), size=cutoff_n)+
            #np.random.normal(loc=0, scale=np.sqrt(np.mean(chunk_uncert**2)/2/len(chunk)), size=cutoff_n) )


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

    def empirical_pvalue(year_idx, point,k,two_sided=True, tail_threshold=1/100*4, num_closest_samples=int(100/7)):
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
            if(sum(~np.isnan(dist))>0 and ~np.isnan(point[i]).all()):
                epdf = stats.gaussian_kde(dist)
                empirical_ll[i] = epdf.logpdf(point[i])
##            if False:
##                plt.figure()
##                xfine = np.linspace(0,1.5,150)
##                plt.hist(dist,density=True)
##                epdf = stats.gaussian_kde(dist)
##                plt.plot(xfine,epdf.pdf(xfine))
##                plt.plot(xfine,stats.norm.pdf(xfine, loc=means[i], scale=ses[i]))
##                plt.plot(xfine,stats.norm.pdf(xfine, loc=means[i], scale=sesl[i]))
##                plt.show()
                
        return empirical_ll    

    def kde_resample(year_idx, nsamps,k):
        if(k!=0 and k!=1):
            return np.full(np.shape(year_idx),np.nan)
        year_idx = np.atleast_1d(year_idx)
        resamps = np.full((len(year_idx),nsamps), np.nan)
        
        for i in range(len(year_idx)):
            if (k==0):
                dist0 = samp_cur[year_idx[i]-1850, :]  # Distribution for the current year
            else:
                return resamps
            dist = dist0[~np.isnan(dist0)]
            if( len(dist)!=0 ):
                epdf = stats.gaussian_kde(dist)
                resamps[i] = epdf.resample(nsamps)

        return resamps
    
    return {
        'mean': empirical_mean,
        'se': empirical_se,
        'pvalue': empirical_pvalue,
        'log_likelihood': empirical_log_likelihood,
        'resample': kde_resample,
    }

    #return means, ses, empser.copy(), empser.copy()

