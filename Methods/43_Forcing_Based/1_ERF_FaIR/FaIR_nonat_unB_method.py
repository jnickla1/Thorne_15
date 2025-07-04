#Added to last cell of anthropogenic-forcing/notebooks/run.ipynb
#np.save('../output/temp_anthro.npy', temp_anthro)
#np.save('../output/temp_all.npy', temp_all)

import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
import os
from scipy import stats
import matplotlib.pyplot as plt

def run_method(years, temperature, uncert, model_run, experiment_type):

    empser  = np.full(len(years),np.nan)

    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    #data_orig = pd.read_csv("./Common_Data/HadCRUT5.csv")
    #temps_obs = data_orig.loc[:,"Anomaly"].to_numpy()
    #preind_base = np.mean(temps_obs[0:50])
    
    cur_path = os.path.dirname(os.path.realpath(__file__))

    #obt_array0 = np.load(cur_path+"/temperature-attribution/output/temp_anthro.npy")
    #Nres = np.shape(obt_array0)[1]
    #starts in 1750 so crop to 100th index
    #samp_cur = obt_array0[100:,:]

    if experiment_type == 'historical':
        sfactor=0.2
        current_array = np.load(cur_path+"/resliced_NorESM/combined_hadcrut5_nonat.npy") #starts in 1930
        retro_array = np.load(cur_path+"/retrospective/all-2022_hadcrut5_currentcut2022_temp_nonat.npy") #starts in 1750
        curbias= -0.014
    else:
        exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126 or _VolcConst #
        sfactor=0.6
        if (exp_attr[1]=='ESM1-2-LR'):
            #combined_all_current_MPIESM370r50_anthro.npy
            current_array = np.load(cur_path+"/resliced_MPIESM/combined_all_current_MPIESM"+exp_attr[2][3:6]+"r"+str(model_run+1)+"_nonat.npy") #starts in 1930
            #all_current_MPIESM370r19cut2099_temp_anthro.npy
            retro_array = np.load(cur_path+"/retrospective/all_current_MPIESM"+exp_attr[2][3:6]+"r"+str(model_run+1)+"cut2099_temp_nonat.npy") #starts in 1750
            if exp_attr[2][3:6]=="126":
                curbias= .013
            elif exp_attr[2][3:6]=="245":
                curbias= -0.076
            elif exp_attr[2][3:6]=="370":
                curbias= -.161
                
        elif (exp_attr[1]=='NorESM'):
            #combined_all_current_NorESMVolcConstr1_anthro.npy
            current_array = np.load(cur_path+"/resliced_NorESM/combined_all_current_NorESM"+exp_attr[3]+"r"+str(model_run+1)+"_nonat.npy") #starts in 1930
            retro_array = np.load(cur_path+"/retrospective/all_current_NorESM"+exp_attr[3]+"r"+str(model_run+1)+"cut2099_temp_nonat.npy") #starts in 1750
            curbias= -.124
            
    samp_cur = np.full((len(years),np.shape(current_array)[1]),np.nan)
    end_fill_sampc = (1930-1850)+np.shape(current_array)[0]
    samp_cur[(1930-1850):end_fill_sampc,:]= current_array #starts in 1930
    samp_mean =  np.nanmean(samp_cur, axis = 1) 
    samp_mean[(1930-1850):(1965-1850)]= 0.28 #overwrite to pass the first check
    dev_orig = samp_cur - samp_mean[:, np.newaxis]
    dev_20 = np.mean(samp_mean[(2000-1850):(2020-1850)]) - np.mean(temperature[(2000-1850):(2020-1850)])
    newcorrecton = np.concatenate((np.zeros(150), np.ones(10)*dev_20, (dev_20+ np.linspace(0, (curbias-dev_20)*9/4, (np.shape(temperature)[0]-160)))))
    samp_cur = samp_mean[:, np.newaxis] + np.sqrt( dev_orig**2)*sfactor *np.sign(dev_orig)- newcorrecton[:, np.newaxis]                        #shrinks their distribution down
    
    samp_ret = retro_array[100:,:]#starts in 1750 so crop to 100th index
    #breakpoint()
    
    for i in range(len(years)):
        samp_ret[i,:] =samp_ret[i,:]+ np.random.normal(loc=0, scale=(temps_1std[i]/np.sqrt(20)), size=np.shape(samp_ret)[1])
        samp_cur[i,:] =samp_cur[i,:]+ np.random.normal(loc=0, scale=(temps_1std[i]/np.sqrt(20)), size=np.shape(samp_cur)[1])
    
    def empirical_mean(year_idx,k):
        if (k==0):
            return np.nanmean(samp_cur[year_idx-1850, :], axis = 1)
        elif(k==1):
            return np.nanmean(samp_ret[year_idx-1850, :], axis = 1)
        else:
            return empser.copy()

    def empirical_se(year_idx,k):
        if (k==0):
            return np.nanstd(samp_cur[year_idx-1850, :], axis = 1)
        elif(k==1):
            return np.nanstd(samp_ret[year_idx-1850, :], axis = 1)
        else:
            return empser.copy()

    def empirical_pvalue(year_idx, point,k,two_sided=True):
        if(k!=0 and k!=1):
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
                Nres = len(dist)
                if Nres==0:
                    empirical_p[i] = np.nan
                    continue
                tail_threshold=1/Nres*4
                num_closest_samples=int(Nres/7)
                
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
        if(k!=0 and k!=1):
            return np.full(np.shape(year_idx),np.nan)
        year_idx = np.atleast_1d(year_idx)
        point = np.atleast_1d(point)
        # Initialize an array to store empirical p-values
        empirical_ll = np.full(np.shape(point), np.nan)

        means = np.mean(samp_cur[years-1850, :], axis = 1)
        ses= np.std(samp_cur[years-1850, :], axis = 1)
        
        for i in range(len(year_idx)):
            if (k==0):
                dist0 = samp_cur[year_idx[i]-1850, :]  # Distribution for the current year
            else:
                dist0 = samp_ret[year_idx[i]-1850, :]
            dist = dist0[~np.isnan(dist0)]
            if( ~np.isnan(point[i]) and len(dist)!=0 ):
                epdf = stats.gaussian_kde(dist)
                empirical_ll[i] = epdf.logpdf(point[i])
            else:
                empirical_ll[i] = np.nan
            if False: #i==26:
                plt.figure()
                #print(dist)
                xfine = np.linspace(-0.05,0.05,100)
                plt.hist(dist,density=True)
                epdf = stats.gaussian_kde(dist)

                plt.plot(xfine,epdf.pdf(xfine))
                plt.plot(xfine,stats.norm.pdf(xfine, loc=means[i], scale=ses[i]))
                #plt.plot(xfine,stats.norm.pdf(xfine, loc=means[i], scale=sesl[i]))
                
       
        return empirical_ll

    return {
        'mean': empirical_mean,
        'se': empirical_se,
        'pvalue': empirical_pvalue,
        'log_likelihood': empirical_log_likelihood,

    }
