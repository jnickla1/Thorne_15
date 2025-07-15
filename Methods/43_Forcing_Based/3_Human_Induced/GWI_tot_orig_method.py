import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
from scipy.integrate import quad
from .GWI_anthro_method import CustomSplineDistribution

import pdb;



percentiles = np.array([5, 17, 50, 83, 95]) / 100  # Convert to 0-1 range

    
cur_path = os.path.dirname(os.path.realpath(__file__))
ord_ind = [0,1,4,2,3]
min_fact = 4

def run_method(years, temperature, uncert, model_run, experiment_type):
    empser  = np.full(np.shape(years),np.nan)
    if experiment_type == 'historical':
        gwi_levels_retro0 = pd.read_csv(cur_path+"/GWI_data/GWI_full_info.csv", header=[0, 1])
        gwi_levels_retro =gwi_levels_retro0.iloc[1:,] #blank row
        gwi_r =gwi_levels_retro['Tot'].to_numpy()
        gwi_levels_curr0 = pd.read_csv(cur_path+"/GWI_data/GWI_hist_only.csv", header=[0, 1])
        gwi_levels_curr =gwi_levels_curr0.iloc[1:,]
        gwi_c =gwi_levels_curr['Tot'].to_numpy()
        lyearr = np.shape(gwi_r)[0]
        lyearc = np.shape(gwi_c)[0]
        lyear = lyearr + 1850
        
    else:
        return empser, empser, empser, empser

    rdists = []
    cdists = []
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4
    
    for i in range(0, lyearr):
        rdists.append(CustomSplineDistribution(ordvalues=gwi_r[i,ord_ind], a=-4, b=7))
    for i in range(0, lyearc):
        cdists.append(CustomSplineDistribution(ordvalues=gwi_c[i,ord_ind], a=-4, b=7))#these start in 1950
        #means_r[i] = float(gwi_r.iloc[i,4]) #50th percentile is 4th column but this is median
        #means_c[i] = float(gwi_c.iloc[1,4])



    def empirical_mean(year_idx,k):


        means = np.full(np.shape(year_idx),np.nan)
        for yr in year_idx:
            if (yr >=lyear or yr <1850):
                means[yr - 1850] =  np.nan
            elif (k==0):
                if (yr >=1950):
                    means[yr - 1850] = cdists[yr-1850-100].mean
                else:
                    means[yr - 1850] = np.nan
            else:
                means[yr - 1850] = rdists[yr-1850].mean


        return means

    def empirical_se(year_idx,k):
        ses = np.full(np.shape(year_idx),np.nan)
        for yr in year_idx:
            if (yr >=lyear or yr <1850):
                ses[yr - 1850] =  np.nan
            elif (k==0):
                if (yr >=1950):
                    ses[yr - 1850] = cdists[yr-1850-100].std
                else:
                    ses[yr - 1850] = np.nan
            else:
                ses[yr - 1850] = rdists[yr-1850].std

            if (ses[yr - 1850] < temps_1std[yr - 1850]/min_fact):
                ses[yr - 1850] = temps_1std[yr - 1850]/min_fact

        return ses

    def empirical_pvalue(year_idx, point,k, two_sided=True):
        year_idx = np.atleast_1d(year_idx)
        point = np.atleast_1d(point) #should be the same dimension
        # Initialize an array to store empirical p-values
        empirical_p = np.full(np.shape(point), np.nan)
        for i, yr in enumerate(year_idx):
            #i = yr-1850
            if (yr >=lyear or yr <1850):
                empirical_p[i] =  np.nan
            elif (k==0):
                if (yr >=1950):
                    
                    cdfpt = cdists[yr-1850-100].cdfn0(point[i])
                    
                    if two_sided:
                        empirical_p[i] = 2*min(1 - cdfpt,cdfpt)
                    else:
                        empirical_p[i] =1-cdfpt
                    #if (cdists[i-100].std < temps_1std[i]/8): #error too small, this does not happen with the current only in the retrospective
                    #    empirical_p[i] = stats.norm.sf(abs((point[i]-cdists[i-100].mean)/ temps_1std[i] * 8))*2

                else:
                    empirical_p[i] =  np.nan

            else:
                cdfpt = rdists[yr-1850].cdfn0(point[i])
                if two_sided:
                    empirical_p[i] = 2*min(1 - cdfpt,cdfpt)
                else:
                    empirical_p[i] = 1- cdfpt
                if (rdists[i].std < temps_1std[i]/min_fact): #error too small, something funny is happening, switch back to norm
                    if two_sided:
                        empirical_p[i] = stats.norm.sf(abs((point[i]-rdists[yr-1850].mean)/ temps_1std[yr-1850] * min_fact))*2
                    else:
                        empirical_p[i] = stats.norm.sf(((point[i]-rdists[yr-1850].mean)/ temps_1std[yr-1850] * min_fact))
            #print(yr)
            #print(empirical_p[i])
        return empirical_p

    def empirical_log_likelihood(year_idx, point,k):
        year_idx = np.atleast_1d(year_idx)
        point = np.atleast_1d(point)
        # Initialize an array to store empirical p-values
        empirical_l = np.full(np.shape(point), np.nan)
        for i, yr in enumerate(year_idx):
            #i = yr-1850
            if (yr >=lyear or yr <1850):
                empirical_l[i] =  np.nan
            elif (k==0):
                if (yr >=1950):
                    empirical_l[i] = (cdists[yr-1850-100].pdfn0(point[i]))
                    #if (cdists[i-100].std < temps_1std[i]/8): #error too small, something funny is happening
                    #    stats.norm.logpdf(point[i],loc=cdists[i-100].mean,scale=temps_1std[i]/8)

                else:
                    empirical_l[i] =  np.nan
            else:
                empirical_l[i] = (rdists[yr-1850].pdfn0(point[i]))
                if (rdists[yr-1850].std < temps_1std[yr-1850]/min_fact): #error too small, something funny is happening
                    empirical_l[i] = stats.norm.pdf(point[i],loc=rdists[yr-1850].mean,scale=temps_1std[yr-1850]/min_fact)
            #print(yr)
            #print(empirical_l[i])     
        return np.log(empirical_l )   

    return {
        'mean': empirical_mean,
        'se': empirical_se,
        'pvalue': empirical_pvalue,
        'log_likelihood': empirical_log_likelihood,

    }

    #return means, ses, empser.copy(), empser.copy()

