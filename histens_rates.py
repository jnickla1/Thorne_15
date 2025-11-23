import os
import config
import importlib
import numpy as np
import pandas as pd
from itertools import chain
import pickle
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import pdb;
import time
start = time.process_time()
import collections
import sys
import importlib

from datetime import datetime
current_date = datetime.now()
formatted_date = current_date.strftime("%y%m%d")

historical_regen=False #whether we are regenerating data (some from saved intermediates or just reading from pickle.
#must be False for this file - reading in lots of pickles




def closest(lst, yrs, K):
    return yrs[np.nanargmin(abs(lst-K))]



from hist_heads_tails_evaluation_script import collect_data, relthresh
import scipy.stats as stats
import statsmodels.stats.multitest as multi
import matplotlib.pyplot as plt
import collections
import mplcursors
spaglines = []
retIDlabels = ['c','r']

alt_colors = ['black', 'white']
sel_methods = [ "CGWL10y_for_halfU","FaIR_nonat","EBMKF_ta2"] #"EBMKF_ta4" "min_month_proj" "OLS_refit_CO2forc", "CGWL10y_for_halfU","TheilSen_h7075" ,"FaIR_anthroA",,"EBMKF_ta2"  ] #"EBMKF_ta",
#sel_methods2 = [ "EBMKF_ta4","GAM_AR1",
#                 "lowess1dg36wnc","Kal_flexLin","FaIR_comb_unB","FaIR_nonat_unB","GWI_anthro_CGWL","CGWL10y_sfUKCP"]


            

#"CGWL10y_for_halfU",,"EBMKF_ta2",

#ftl = np.argsort(index_mapping) #from to list - where a certain method should be plotted

def local_linear_slope_and_se(y_segment, t_segment):
    """
    Compute OLS slope and its standard error for y vs t in a local window.
    Returns (slope_per_year, se_slope_per_year).
    """
    y_segment = np.asarray(y_segment)
    t_segment = np.asarray(t_segment)

    mask = ~np.isnan(y_segment)
    y = y_segment[mask]
    t = t_segment[mask]

    n = len(y)
    if n < 3:
        return np.nan, np.nan

    t_mean = np.mean(t)
    y_mean = np.mean(y)

    Sxx = np.sum((t - t_mean)**2)
    if Sxx == 0:
        return np.nan, np.nan

    # OLS slope
    beta1 = np.sum((t - t_mean) * (y - y_mean)) / Sxx

    # Residuals and variance
    residuals = y - (beta1 * t + (y_mean - beta1 * t_mean))
    dof = n - 2
    if dof <= 0:
        return beta1, np.nan

    sigma2 = np.sum(residuals**2) / dof
    se_beta1 = np.sqrt(sigma2 / Sxx)  # standard error of the slope

    return beta1, se_beta1


if __name__ == '__main__':

    experiment_type = sys.argv[1] #'fut_ESM1-2-LR_SSP126_constVolc' #fut_NorESM_RCP45_Volc
    model_run = 0
    exp_attr = experiment_type.split("_")
    results_path = f'Results/results{experiment_type}{model_run}.pickle'
    years = np.arange(1850,2025)
#To read it back:

    if exp_attr[2]=="real":
        from plot_fut_results import sel_methods as sel_methods2
    else:
        from plot_fut_results import sel_methodsA as sel_methods2


    with open (results_path, 'rb') as fp:
        results = pickle.load(fp)

##    pcross_path = 'Results/headstails_method_fineprobs_'+experiment_type+str(model_run)+'.npy'
##    with open (pcross_path, 'rb') as fp:
##        pcrossmatrix = np.load(pcross_path)

    #heights0=np.arange(.2,4.2,0.05) #heights to test at thresholds: 0.05°C increments
    #nthresholds = np.size(heights0)
    inum=12 #internal interpolation within years

    #CALCULATE RATES ... need to push cent10y through again

    (sims_tas, stimes) = collect_data(exp_attr)
    simall = sims_tas[:,model_run]
    years = stimes.astype(int)
    data = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data.loc[:,"Anomaly"].to_numpy()
    temps_CIu =data.loc[:,"Upper"].to_numpy() #Lower confidence limit (2.5%)	Upper confidence limit (97.5%)
    temps_CIl =data.loc[:,"Lower"].to_numpy()
    temps_CIl_hist = simall[0:len(years)]+  temps_CIl - temps_obs
    temps_CIu_hist = simall[0:len(years)]+  temps_CIu - temps_obs
    

    if exp_attr[2]=='real':
        standard = results['cent20y']['LT_trend'][2] #retrospective
    else:
        standardmethod_module = importlib.import_module("Methods.43_Forcing_Based.1_ERF_FaIR.FaIR_anthro_unB_method", package = "Thorne_15_codefigurestats")
        result = standardmethod_module.run_method(years, simall,(temps_CIl_hist,temps_CIu_hist),model_run, experiment_type)
        standard = result['mean'](years,1)#retrospective
        standardunc = result['se'](years,1) #retrospective
        
    
    #breakpoint()
    smooth_std = np.nanmean(np.abs(np.diff(np.diff(standard))))

    sthresh = 0. if exp_attr[1]=="rpi" else (relthresh+0.)
    eval05min=1960
    eval10max=2015
    fineyrs = np.arange(eval05min,eval10max+1/inum,1/inum) 
    std_intp0 = np.interp(fineyrs,years,standard)
    std_intp = std_intp0[~np.isnan(std_intp0)] #all true are at the end so still indexing the same thing
    closest_year05 = fineyrs[np.argmin(np.abs(std_intp - 0.5-sthresh))]
    closest_year10 = fineyrs[np.argmin(np.abs(std_intp - 1.0-sthresh))]
    #still have these for the time of crossing
    closest_years = [fineyrs[np.argmin(np.abs(std_intp - i))] for i in np.arange(0.5,1.05,.1)]



    #ax1.hlines(0.5,xmin=eval05min,xmax = eval05max-15, color='black', linestyle='--', linewidth=2 )

    
    central_yr_estimates =[]


    i=0
    ci=0
    labelcolors=[]
    import neworder
    sorted_results = neworder.sort_results(results)


    ncm = 0 #number of current methods
    for method_name, method_data in sorted_results:
        for k in range(2):
            result = method_data['LT_trend']
            central_est=np.full(len(years),np.nan) #make this blank to start
            labelcurr_or_retro="" 
            if isinstance(result, dict):
                # Empirical method: Call the functions at each time point
                central_est = result['mean'](years,k)
            else:
                central_est = result[k*2]
            if (sum(~np.isnan(central_est))>0 ):    #isinstance(result, tuple):
                # Central estimate and SE, implying gaussian distn
                labelcurr_or_retro = retIDlabels[k]
            if(labelcurr_or_retro=="c"):
                ncm=ncm+1
    print(ncm)

    #if (len(index_mapping) != ncm): #not computing all methods, for debugging only
    #    index_mapping = np.arange(ncm)
    #    ftl = np.argsort(index_mapping)
        
    sel2names = []
    sel2data = []
    for method_name, method_data in sorted_results:
        #print(f"Results from {method_name} (Method Class: {method_data['method_class']}):")
        result = method_data['LT_trend']
        labelcurr_or_retro="" #only set if the current method is not blank
        #print(method_name)
        if method_name in sel_methods2:
            can_data = results[method_name]['LT_trend'][2]
            
            if(sum(~np.isnan(can_data))>0):
                sel2data.append(results[method_name]['LT_trend'][2][-15:])
                sel2names.append(method_name)
            else:
                sel2data.append(results[method_name]['LT_trend'][0][-15:]) #appending purely current method instead
                print(method_name+" current only")
                sel2names.append(method_name)
        for k in range(2):
            central_est=np.full(len(years),np.nan) #make this blank to start
            if isinstance(result, dict):
                # Empirical method: Call the functions at each time point
                central_est = result['mean'](years,k)
                se = result['se'](years,k)
                pvals = result['pvalue'](years, standard,k)
                llikelihood = result['log_likelihood'](years, standard,k)
                llikelihood2 = stats.norm.logpdf(standard,loc=central_est,scale=se)
                
                if(sum(~np.isnan(central_est))>0):
                    print(f"{method_name} sampled, d llike: {np.nanmean(llikelihood)-np.nanmean(llikelihood2)}")
                #breakpoint() 
                    
            else:
                central_est = result[k*2]
                se = result[k*2+1]
                # want to save and store the more accurate calculation if regen: else: load calculations
                pvals = stats.norm.sf(abs((standard-central_est)/ se))*2
                llikelihood = stats.norm.logpdf(standard,loc=central_est,scale=se)

    
        #    percentiles = result
           # print(f"  Percentiles: {percentiles}"
    


    if exp_attr[2]=='real':
        standard10method_module = importlib.import_module("Methods.42_Temp_Alone.1_Run_Means.cent10y_method", package = "Thorne_15_codefigurestats")
        result10 = standard10method_module.run_method(years, simall,(temps_CIl_hist,temps_CIu_hist),model_run, experiment_type)
        standard10 = result10[2] #retrospective
        standard10unc = result10[3] #retrospective
        stdname="cent10y_20y"
    else:
        standard10 = standard
        standard10unc = standardunc
        stdname='FaIR_anthro_unB'

    eval05min=1960
    eval10max=2020
    fineyrs = np.arange(eval05min,eval10max+1/inum,1/inum)
    
    std_intp10n = np.interp(fineyrs,years,standard10)
    std_intp10uncn = np.interp(fineyrs,years,standard10unc)
    std_intp10 = std_intp10n[~np.isnan(std_intp10n)] #all nans are at the end so still indexing the same thing
    std_intp10unc = std_intp10uncn[~np.isnan(std_intp10uncn)]
    closest_ind05 = np.argmin(np.abs(std_intp - 0.5-sthresh))
    closest_ind10 = np.argmin(np.abs(std_intp - 1.0-sthresh))

    #rate05a = std_intp10[closest_ind05 + int(5*inum)] - std_intp10[closest_ind05 - int(5*inum)]
    #rate5a = standard10[int(np.round(closest_year05))-1850+5] - standard10[int(np.round(closest_year05))-1850 - 5]
    #uncert05a = 1.96* np.sqrt(std_intp10unc[closest_ind05 + int(5*inum)]**2 + std_intp10unc[closest_ind05 - int(5*inum)]**2)
    #now using local linear OLS standard error instead
    yrs_05wings=5
    yrs_10wings = 5
    try:
        idx10_offset = int(yrs_10wings * inum)
        rate10a = std_intp10[closest_ind10 + idx10_offset] - std_intp10[closest_ind10 - idx10_offset]
    except:
        print("Adjusting 1.0C rate window to += :")
        # there's not been a full 10 years since this threshold has been exceeded
        yrs_10wings = (len(std_intp10) - closest_ind10) / float(inum)
        idx10_offset = int(yrs_10wings * inum)
        print(yrs_10wings)
        rate10a = std_intp10[closest_ind10 + idx10_offset] - std_intp10[closest_ind10 - idx10_offset]
        
    start05 = max(0, closest_ind05 - int(5*inum))
    end05   = min(len(std_intp10) - 1, closest_ind05 + int(5*inum))
    y05 = std_intp10[start05:end05 + 1:inum]
    t05 = fineyrs[start05:end05 + 1:inum]
    slope05, se_slope05 = local_linear_slope_and_se(y05, t05)  # per year
    # Convert slope SE to per decade and use 1.96 for ~95% CI
    rate05a = slope05 *10
    uncert05a = 1.96 * se_slope05 * 10.0

    start10 = max(0, closest_ind10 - idx10_offset)
    end10   = min(len(std_intp10) - 1, closest_ind10 + idx10_offset)
    y10 = std_intp10[start10:end10 + 1:inum]
    t10 = fineyrs[start10:end10 + 1:inum]

    slope10, se_slope10 = local_linear_slope_and_se(y10, t10)  # per year
    uncert10a = 1.96 * se_slope10 * 10.0
    # Convert slope SE to per decade and use 1.96 for ~95% CI
    rate05a = slope05 *10
    uncert05a = 1.96 * se_slope05 * 10.0
    

    #print(f'Decadal rate in {(closest_year05):.1f}  at 0.5°C: {(rate05a):.4f} ± {(uncert05a):.4f}')
##    yrs_10wings = 5
##    try:
##        rate10a = std_intp10[closest_ind10 + int(yrs_10wings*inum)] - std_intp10[closest_ind10 - int(yrs_10wings*inum)]
##        uncert10a = 1.96* np.sqrt(std_intp10unc[closest_ind10 + int(yrs_10wings*inum)]**2 + std_intp10unc[closest_ind10 - int(yrs_10wings*inum)]**2)
##    except:
##        print("Ajusting 1.0C rate window to += :")
##        #there's not been a full 10 years since this threshold has been exceeded
##        yrs_10wings = (len(std_intp10) - closest_ind10)/float(inum)
##        print(years_10wings)
##        rate10a = std_intp10[closest_ind10 + int(yrs_10wings*inum)] - std_intp10[closest_ind10 - int(yrs_10wings*inum)]
##        uncert10a = 1.96* np.sqrt(std_intp10unc[closest_ind10 + int(yrs_10wings*inum)]**2 + std_intp10unc[closest_ind10 - int(yrs_10wings*inum)]**2)
        

    #print(f'Decadal rate in {(closest_year10):.1f} at 1.0°C: {(rate10a):.4f} ± {(uncert10a):.4f}')


    #plt.figure()
    #plt.plot(standard10)
    #plt.show()


    

    tslope = np.arange(15)

    def fit_quad(y):
        mask = ~np.isnan(y)
        X = np.vstack([tslope[mask]**2, tslope[mask], np.ones(np.sum(mask))]).T
        beta, residuals, rank, s = np.linalg.lstsq(X, y[mask], rcond=None)

        dof = np.sum(mask) - 3
        mse = residuals / dof if dof > 0 else 0
        cov = mse * np.linalg.inv(X.T @ X)
        plt.plot(tslope[-2:],(beta[0]*tslope[-2:]**2 + beta[1]*tslope[-2:] + beta[2]), label=None,lw=4,color='k')
        # Slope at t=14 and its variance
        t_last = 14
        slope = 2 * beta[0] * t_last + beta[1]
        slope_var = (2 * t_last)**2 * cov[0, 0] + cov[1, 1] + 2 * (2 * t_last) * cov[0, 1]
        return slope, slope_var

    # Fit all and collect slopes and variances
    slopes_vars = np.array([fit_quad(y) for y in sel2data])
    #breakpoint()
    slopes = slopes_vars[:, 0]
    vars_ = slopes_vars[:, 1]
    
    slope_combined = np.mean(slopes) *10
    # Uncertainty of the combined slope (standard error)
    slope_combined_std = np.sqrt(np.mean(vars_)+np.var(slopes)) *10 *1.96
    #print(f"Decadal rate in 2024: {slope_combined:.4f} ± {slope_combined_std:.4f}")
    
    #df_results.to_csv('method_statistics_results.csv', index=False)
    #ax1.legend(fontsize=6.5,ncol=4)
    
    # ------------------------------------------------------------------
    # Build dataframe for frontend interface with decadal rates
    # Columns: level, method, year_cross, rate, rate_uncert
    # ------------------------------------------------------------------

    frontend_rows = [
        {
            "level": "0.5°C",
            "method": stdname,
            "year_cross": closest_year05,
            "rate": rate05a,
            "rate_uncert": uncert05a,
            "rate_window": "+/-"+str(yrs_05wings),
        },
        {
            "level": "1.0°C",
            "method": stdname,
            "year_cross": closest_year10,
            "rate": rate10a,
            "rate_uncert": uncert10a,
            "rate_window": "+/-"+str(yrs_10wings),
        },
    ]

    df_frontend = pd.DataFrame(
        frontend_rows,
        columns=["level", "method", "year_cross", "rate", "rate_uncert","rate_window"],
    )

    method_rows = []
    sel3data=np.array(sel2data)
    for mm, name in enumerate(sel2names):
        # Reported warming level in 2024 for this method

        row = {
            "level": "2024",
            "method": name,
            "year_cross": sel3data[mm,-1],
            "rate": float(slopes[mm]*10),
            "rate_uncert": np.sqrt(vars_[mm])*10*1.96,
            "rate_window": -15
        }
        method_rows.append(row)

    # ------------------------------------------------------------------
    # Add equal-weight combined row for 2024
    # ------------------------------------------------------------------
    warm_2024_comb = np.nanmean(sel3data[:,-1])

    comb_row = {
        "level": "2024",
        "method": "equal_weight_comb",
        "year_cross": warm_2024_comb,
        "rate": slope_combined ,
        "rate_uncert": slope_combined_std,
        "rate_window": -15
    }

    method_rows.append(comb_row)


    # ------------------------------------------------------------------
    # Add IV-weight combined row for 2024
    # ------------------------------------------------------------------


    df_results = pd.read_csv(f'Results2/historical_names_var_scale.csv', index_col=0)
    all_methods = df_results['method_name']
    err_vars = df_results['err_var'].to_numpy()
    err_vars100 = df_results['err_var100'].to_numpy()
    #best_scales = df_results['best_alter_scale'].to_numpy()
    #central_est = np.load(f'Results2/{outputfilename}_central_est.npy')
    #node_lls = np.load(f'Results2/{outputfilename}_hermguass_lls.npy')
    #newsamples = np.load(f'Results2/{outputfilename}_newsamples.npy')
    #nyears = central_est.shape[1]
    mask_a = all_methods.isin(sel2names).values
    vars100 = err_vars100[mask_a]
    ivweights = 1/vars100
    ivcenter = np.nansum(sel3data[:,-1] * ivweights[:],axis=0) / np.nansum( ~np.isnan(sel3data[:,-1]) * ivweights[:],axis=0)
    ivslope_combined = np.nansum(slopes * ivweights[:],axis=0) / np.nansum( ~np.isnan(slopes) * ivweights[:],axis=0)*10
    # Uncertainty of the combined slope (standard error)
    ivslope_combined_std = np.sqrt( np.nansum(vars_ * ivweights[:],axis=0) / np.nansum( ~np.isnan(vars_) * ivweights[:],axis=0)
                                  +np.nansum(((slopes-ivslope_combined/10)**2 * ivweights),axis=0) / np.nansum( ~np.isnan(slopes) * ivweights[:],axis=0)) *10 *1.96    

    
         
##    warm_2024_comb = np.nanmean(sel3data[:,-1])
##
##    ivcenters = np.nansum(centralsh * ivweights[:,None],axis=0) / np.nansum( ~np.isnan(centralsh) * ivweights[:,None],axis=0)
##    
##
    ivcomb_row = {
        "level": "2024",
        "method": "PIVW_comb",
        "year_cross": ivcenter,
        "rate": ivslope_combined ,
        "rate_uncert": ivslope_combined_std,
        "rate_window": -15
    }

    method_rows.append(ivcomb_row)


#Implementing CSCM is more complex than what I have time for right now, since we need newsamples for each of the histens members
##newsamples = np.load(f'Results2/{outputfilename}_newsamples.npy')
##
##samplesh = newsamples[mask_a]
##
##        alpha = w / np.sqrt(np.pi)
##        infilled_lls = infill_log_likelihoods(lls[:,10:-10,:],penalty=0)
##        stack_weights = fit_stacking_weights_logp(infilled_lls,alpha)
###blended = predict_non_nan(stack_weights, samplesh)
###sharpened_blended = sharpen_samples(blended, best_nsharp, seed = 2)
###sharp_blend_central = np.mean(sharpened_blended,axis=1)

    df_methods_2024 = pd.DataFrame(
        method_rows,
        columns=["level", "method", "year_cross", "rate", "rate_uncert","rate_window"]
    )

    df_frontend = pd.concat([df_frontend, df_methods_2024], ignore_index=True)


    print(df_frontend)
    print(time.process_time() - start)
    
    
   #sorted_df = df_res_cur2.reset_index(drop=True).sort_values(by=['method_class', 'bias50']).reset_index()
    #sorted_df[['index']].to_csv('to_index_mapping.csv', index=False)

