import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
from scipy.integrate import quad
from .GWI_anthro_method import CustomSplineDistribution, gen_orig_number

import pdb;



percentiles = np.array([5, 17, 50, 83, 95]) / 100  # Convert to 0-1 range

    
cur_path = os.path.dirname(os.path.realpath(__file__))

ord_ind = [0,1,4,2,3]
min_fact = 4

def run_method(years, temperature, uncert, model_run, experiment_type):
    #ONLY DEFINED NOW FOR CURRENT WARMING LEVEL, so k must be 0 or will return just NANs
    syear = 1950
    if experiment_type == 'historical':
        gwi_levels_curr0 = pd.read_csv(cur_path+"/Thorne2025_GWI_Results/AR6_ESM1-2-LR/"+
                "GWI_results_AR6_HISTORICAL-ONLY_SCENARIO--observed-SSP245_ENSEMBLE-MEMBER-"+
                                       "-all_VARIABLES--GHG-Nat-OHF___REGRESSED-YEARS--1850-1950_to_1850-2024.csv", header=[0, 1])
        
    else:
        #future case, grabbing the ANNUAL resutls
        exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126_constVolc #

        if (exp_attr[1]=='ESM1-2-LR'):
            gwi_levels_curr0 = pd.read_csv(cur_path+"/Thorne2025_GWI_Results/AR6_ESM1-2-LR/"+
                "GWI_results_AR6_HISTORICAL-ONLY_SCENARIO--SMILE_ESM-"+exp_attr[2]+"_ENSEMBLE-MEMBER--"+
                                           str(model_run)+"_VARIABLES--GHG-Nat-OHF___REGRESSED-YEARS--1850-1950_to_1850-2100.csv", header=[0, 1])
            
        elif (exp_attr[1]=='NorESM'):
            model_run_noresm = gen_orig_number(model_run,60)
            gwi_levels_curr0 = pd.read_csv(cur_path+"/Thorne2025_GWI_Results/AR6_NorESM/"+
                "GWI_results_AR6_HISTORICAL-ONLY_SCENARIO--NorESM_rcp45-"+exp_attr[3]+"_ENSEMBLE-MEMBER--"+
                                           str(model_run_noresm)+"_VARIABLES--GHG-Nat-OHF___REGRESSED-YEARS--1850-1950_to_1850-2099.csv", header=[0, 1])

    gwi_levels_curr =gwi_levels_curr0.iloc[1:,]
    gwi_c =gwi_levels_curr['Ant'].to_numpy()
    lyearc = np.shape(gwi_c)[0]
    lyear = lyearc + syear

    rdists = []
    cdists = []
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4
    
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
                if (yr >=syear):
                    means[yr - 1850] = cdists[yr-syear].mean
                else:
                    means[yr - 1850] = np.nan

        return means

    def empirical_se(year_idx,k):
        ses = np.full(np.shape(year_idx),np.nan)
        for yr in year_idx:
            if (yr >=lyear or yr <1850):
                ses[yr - 1850] =  np.nan
            elif (k==0):
                if (yr >=syear):
                    ses[yr - 1850] = cdists[yr-syear].std
                else:
                    ses[yr - 1850] = np.nan

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
                if (yr >=syear):
                    
                    cdfpt = cdists[yr-syear].cdfn0(point[i])
                    
                    if two_sided:
                        empirical_p[i] = 2*min(1 - cdfpt,cdfpt)
                    else:
                        empirical_p[i] =1-cdfpt

                else:
                    empirical_p[i] =  np.nan

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
                if (yr >=syear):
                    empirical_l[i] = (cdists[yr-syear].pdfn0(point[i]))
                else:
                    empirical_l[i] =  np.nan

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

