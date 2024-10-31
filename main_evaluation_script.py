import os
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

regen=True #whether to execute all of the methods directly (True) or read from a pickle file to rapidly redraw the figures (False).

import colorsys

def gen_color(ci, dark=False):
    """
    Generate a distinct color based on an integer index, ci.
    
    Parameters:
    - ci: int, the index for generating different colors.
    - dark: bool, if True, returns a darker version of the color; otherwise, returns a lighter version.
    
    Returns:
    - tuple: RGB color as (r, g, b), values between 0 and 1.
    """
    # Set base hue by cycling through color wheel
    hue = (ci * 137.5) % 360  # Golden angle in degrees for distinct hues
    hue /= 360  # Normalize to [0, 1] range
    
    # Set saturation and value (brightness)
    saturation = 0.9
    value = 0.5 if dark else 0.8  # Adjust value for lighter or darker shades
    
    # Convert HSV to RGB
    r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
    return (r, g, b)



def run_methods(years, avg_temperatures, temp_uncert,model_run, experiment_type, methods_folder=('Methods/42_Temp_Alone/1_Run_Means','Methods/42_Temp_Alone/6_Remove_IV')):
    #'Methods/44_EarthModel_CGWL','Methods/42_Temp_Alone/3_ST_Fits' 'Methods/42_Temp_Alone/6_Remove_IV'
#'Methods/43_Forcing_Based/1_ERF_FaIR','Methods/43_Forcing_Based/3_Human_Induced' 'Methods/42_Temp_Alone/1_Run_Means', 'Methods/43_Forcing_Based/3_Human_Induced' 
#'Methods/42_Temp_Alone/1_Run_Means','Methods/44_EarthModel_CGWL')): 
#'Methods/42_Temp_Alone','Methods/43_Forcing_Based'     'Methods/42_Temp_Alone','Methods/43_Forcing_Based','Methods/44_EarthModel_CGWL'
    # Find all *_method.py files in the folder and subfolders
    method_files = []
    for root, _, files in chain.from_iterable(os.walk(path) for path in methods_folder):
        for file in files:
            if file.endswith('_method.py'):
                method_files.append(os.path.join(root, file))

    # Initialize result storage
    results = {}

    # Send data to each method and retrieve results
    for method_path in method_files:
        # Extract method class (folder structure relative to methods_folder)
        current_file_path = os.getcwd()+"/Methods"
        method_class = os.path.relpath(os.path.dirname(method_path), current_file_path) #methods_folder)

        # Dynamically import the module
        print(method_path)
        module_name = method_path.replace('/', '.').replace('.py', '')
        #print(module_name)
        method_module = importlib.import_module(module_name, package = "Thorne_15_codefigurestats")

        # Call the method function (assumed to be named "run_method" in each *_method.py)
        result = method_module.run_method(years, avg_temperatures, temp_uncert,model_run, experiment_type)

        # Store result along with method class (folder structure)
        method_name = os.path.basename(method_path).replace('_method.py', '')
        results[method_name] = {
            'method_class': method_class,
            'LT_trend': result
        }
    return results

def closest(lst, yrs, K):
    return yrs[np.nanargmin(abs(lst-K))]

retIDlabels = ['c','r']

import scipy.stats as stats
import statsmodels.stats.multitest as multi
import matplotlib.pyplot as plt

alt_colors = ['black', 'dimgray']

if __name__ == '__main__':
    # First evaluation
    data = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50]) #preindustrial baseline 1850-1899
    temps_obs = temps_obs - preind_base #remove this baseline
    temps_CIu =data.loc[:,"Upper"].to_numpy() #Lower confidence limit (2.5%)	Upper confidence limit (97.5%)
    temps_CIl =data.loc[:,"Lower"].to_numpy()
    years=data.loc[:,"Time"].to_numpy()
    nyrs = len(years)
    model_run = 'hadcrut5'
    experiment_type = 'historical'

    # Run all the methods

    if regen:
        results = run_methods(years, temps_obs, (temps_CIl, temps_CIu),model_run, experiment_type)
        with open('hist_results.pickle', 'wb') as fp:
            results2 = {}
            for method_name, method_data in results.items():
                result = method_data['LT_trend']
                method_class = method_data['method_class']
                
                if isinstance(result, dict): #remove the dynamically sampled aspect so the whole thing is an array
                    result2=(result['mean'](years,0),result['se'](years,0),result['mean'](years,1),result['se'](years,1))
                else:
                    result2=result
                
                results2[method_name] = {'method_class': method_class,'LT_trend': result2}
            pickle.dump(results2, fp)
#To read it back:
    else:
        with open ('hist_results.pickle', 'rb') as fp:
            results = pickle.load(fp)
    

    heights0=np.arange(.2,4.2,0.05) #heights to test at thresholds: 0.05°C increments
    nthresholds = np.size(heights0)
    inum=12 #internal interpolation within years
    standard = results['cent20y']['LT_trend'][2] #retrospective
    smooth_std = np.nanmean(np.abs(np.diff(np.diff(standard))))

    fig, (ax1,ax4)= plt.subplots(2, 1, figsize=(10,10), gridspec_kw={ "hspace": 0.3})
    ax1.plot(years, standard, 'ko',zorder=3)

    fig2, (ax05, hist_ax)  = plt.subplots(2, 1, figsize=(7,10),gridspec_kw={'height_ratios': [4, 1] , 'left':0.28 , 'right':0.9 , 'top':0.88, 'bottom':0.12}) #, sharex=True)

    eval05min=1970
    eval05max=2015
    eval05yrs = np.arange(eval05min,eval05max)
    fineyrs = np.arange(eval05min,eval05max+1/inum,1/inum) 
    std_intp0 = np.interp(fineyrs,years,standard)
    std_intp = std_intp0[~np.isnan(std_intp0)]
    closest_idx = np.argmin(np.abs(std_intp - 0.5))
    # Retrieve the corresponding year
    closest_year = fineyrs[closest_idx]

    eval05 = {}
    crossing_yrs =[]
    cmethods = []
    ax05.axvline(closest_year, color='black', linestyle='--', linewidth=2 ) #, label="Gold Standard (1986)")
    hist_ax.axvline(closest_year, color='black', linestyle='--', linewidth=2 ) #, label="Gold Standard (1986)")
    central_yr_estimates =[]
    df_results = pd.DataFrame(columns=['method_name', 'method_class','c/r','smooth_r','avg_unc.(1se)','#q<0.5', '#q<0.1', 'q_min', 'q_small5','log-likeli','RMS','bias','tlog-l','100log-l'])
    i=0
    ci=0
    for method_name, method_data in results.items():
        #print(f"Results from {method_name} (Method Class: {method_data['method_class']}):")
        result = method_data['LT_trend']
        labelcurr_or_retro="" #only set if the current method is not blank
        for k in range(2):
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

                
            if (sum(~np.isnan(central_est))>0 ):    #isinstance(result, tuple):
                # Central estimate and SE, implying gaussian distn
                labelcurr_or_retro = retIDlabels[k]
                # Perform FDR adjustment and other calculations...
                smooth_est = np.nanmean(np.abs(np.diff(np.diff(central_est))))
                nnans = ~(np.isnan(pvals))
                qvals = multi.multipletests(pvals[nnans], alpha=0.5, method='fdr_bh')
                qvals_count_yrs05 = np.sum(qvals[0])
                qvals_count_yrs01 = np.sum(qvals[1]<0.1)
                qvals_smallest = np.min(qvals[1])
                qvals_smallest5 = np.sort(qvals[1])[4]
                avg_uncert = np.nanmean(se)
                #print(qvals_count_yrs ,qvals_smallest, qvals_smallest5 )
                
    #Line PLOT TO SHOW WHAT IT'S DOING

                if(labelcurr_or_retro=="c"):
                    ax1.plot(years, central_est, label=method_name, color = gen_color(ci, dark=False))
                    ax4.plot(years, central_est-standard,alpha=0.1, color = gen_color(ci, dark=False))
                    if(method_name in [ "EBMKF_ta", "CGWL10y_for_halfU","TheilSen_h7075"  ]): #removeMEI_volc_refit"): "CGWL10y_sUKCP", "KCC_human", "CGWL10y_for_hU","30y_etrend3CS","CGWL_10y_IPCC","lfca_hadcrut"
                        #ax1.fill_between(years, central_est-se, central_est+se, alpha=0.6)
                        #ax1.plot(years, central_est, lw=5,)
                        ax4.fill_between(years, central_est-se-standard, central_est+se-standard, alpha=0.3,label=method_name, color = gen_color(ci, dark=False))
                        ax4.plot(years, central_est-standard,label=method_name, color = gen_color(ci, dark=True))
                        #print(method_name,np.nanmean(central_est-standard))
                        
     #Histogram plots at 0.5°C
                if(labelcurr_or_retro=="c"):
                    #first compute the p-vals but one-sided
                    this_method_p_steps = np.full(np.shape(eval05yrs),np.nan)
                    if isinstance(result, dict):
                        this_method_p_steps = result['pvalue'](eval05yrs, np.full(np.shape(eval05yrs),0.5),k, two_sided=False)
                    else:
                        this_method_p_steps = stats.norm.cdf(((central_est[eval05yrs-1850]-0.5)/ se[eval05yrs-1850]))
                   # print(method_name)
                   # print(this_method_p_steps)
                    this_method_crossing_p = np.diff(this_method_p_steps) #probability of crossing in a particular year,
                    halfyrs = eval05yrs[1:]+0.5   #will plot on half-years instead
                    fineeval05yrs = np.arange(eval05yrs[0],eval05yrs[-1]+1/inum,1/inum) 
                    psteps_intp = np.interp(fineeval05yrs,eval05yrs,this_method_p_steps)
                    integrated_diffs = np.cumsum(psteps_intp) - np.flip(np.cumsum(1-np.flip(psteps_intp)))
                    crossing_start = fineeval05yrs[np.argmax(psteps_intp>0.1587)] #interpolated first moment that this_method_p_steps exceeds 0.16
                    #crossing_exp_value = np.sum(this_method_crossing_p * halfyrs) #cant just do this due to many negatives
                    crossing_exp_value = fineeval05yrs[np.argmin(abs(integrated_diffs))]
                    crossing_end =  fineeval05yrs[len(psteps_intp) -1 - np.argmax(psteps_intp[::-1]<(1-0.1587))]#interpolated last moment that this_method_p_steps exceeds 0.84

                    crossing_yrs.append(crossing_exp_value)
                    cmethods.append(method_name)
                    crossing_p_pos = this_method_crossing_p * (this_method_crossing_p>0)
                    crossing_p_pos = crossing_p_pos  / np.sum(crossing_p_pos) #normalize to a total area of 4
                    ax05.fill_between(halfyrs, ci - crossing_p_pos, ci + crossing_p_pos,color=gen_color(ci, dark=False), alpha=0.6, edgecolor='none')
                    xerr_arr = np.array([[max(crossing_exp_value  - crossing_start,0), max(crossing_end -crossing_exp_value,0) ]])
                    ax05.errorbar(crossing_exp_value, ci, xerr=xerr_arr.T, fmt='o', color=alt_colors[ci % 2], capsize=3)

                    eval05[method_name] = [this_method_p_steps, this_method_p_steps,crossing_start,crossing_exp_value , crossing_end] #save for later
                    ci = ci+1 #increment current method counter
                    
                    
              #  elif("lowess" in method_name):
              #      plt.fill_between(years, central_est-se, central_est+se, alpha=0.6)

                short_method_class = method_data['method_class'][0:2] +"/"+method_data['method_class'].split('/')[-1]
                df_results.loc[i]= [ method_name,short_method_class,labelcurr_or_retro,smooth_est/smooth_std,avg_uncert,
                                 qvals_count_yrs05,qvals_count_yrs01,  qvals_smallest,qvals_smallest5, np.nanmean(llikelihood), np.sqrt(np.nanmean((central_est-standard)**2)),
                                   np.nanmean(central_est-standard) , np.nansum(llikelihood),np.nansum(llikelihood[-100:-1]) ] #nansum(llikelihood)
                i=i+1
   
       # else:
            # 100 percentiles
        #    percentiles = result
           # print(f"  Percentiles: {percentiles}"
    ax05.set_yticks(range(len(cmethods)))
    ax05.set_yticklabels(cmethods , fontsize=8)
    for label, color in zip(ax05.get_yticklabels(), alt_colors * (len(cmethods) // 2 + 1)):
        label.set_color(color)

    hist_ax.set_xlim([eval05min, eval05max-10])
    ax05.set_xlim([eval05min, eval05max-10])
    
    ax05.set_xlabel("Year")
    ax05.set_xlim([eval05min, eval05max-10])
    ax05.set_title("Crossing Years for 0.5°C Above Preindustrial by Method")
    #ax05.legend(loc='upper right')
    ax05.set_xticks( hist_ax.get_xticks())
    # Histogram of crossing years at the bottom

    hist_ax.hist(crossing_yrs, color="skyblue", edgecolor='grey')
    hist_ax.set_ylabel("Count")
    hist_ax.set_xlabel("Year")


    df_res_show = df_results.drop(columns=['tlog-l','100log-l'])
    df_res_show['smooth_r'] = df_results['smooth_r'].round(3)
    df_res_cur = df_res_show[df_results['c/r']=='c']
    
    print(df_res_cur.sort_values('log-likeli',ascending=False))

    #df_results.to_csv('method_statistics_results.csv', index=False)
    ax1.legend(fontsize=7,ncol=4)
    ax4.set_ylim(bottom=-0.11,top=0.11)
    xmin, xmax = ax1.get_xlim()
    ax4.set_xlim([xmin, xmax])
    ax4.axhline(y=0, color='k', linestyle='--')
    ax4.legend()
    ax1.set_xlabel("Year")
    ax4.set_xlabel("Year")
    ax1.set_ylabel("Temperature °C Anomaly \n relative to 1850-1900")
    ax4.set_ylabel("Temperature °C Difference \n relative to 20-yr running mean")
    print(time.process_time() - start)
    plt.show()
    df_res_show2 = df_results.copy()
    df_res_show2['smooth_r'] = df_results['smooth_r'].round(3)
    df_res_cur2 = df_res_show2[df_results['c/r']=='c']
    dfres2 = df_res_cur2.sort_values('log-likeli',ascending=False)
    #dfres2.to_csv('current_methods_statistics_241031.csv', index=False)








        ##Process along time marginal for each temp rather than along the temperature marginal for each date
##    #make this monotonically increasing so we are not assessing back and forth thresholds
##    
##    std_increas = standard.copy() 
##    srmax = 0
##    for i in range(len(standard)):
##        if standard[i]>srmax:
##            srmax=standard[i]
##        else:
##            std_increas[i]=np.NaN
##   
##    fineyrs = np.arange(1850,years[-1]+1/inum,1/inum) #included 2024
##    std_intp = np.interp(fineyrs,years,std_increas)
##
##    #find closest fineyears this standard method actually crossed
##
##    crossyears = np.zeros(np.shape(heights0))
##    for i in range(nthresholds):
##        h = heights0[i]
##        tryyr = closest(std_intp, fineyrs, h)
##        if (tryyr>1920 and tryyr<(years[-1]-10)):
##            crossyears[i]=tryyr
##        else:
##            crossyears[i]=np.NaN
##    #plt.plot(fineyrs, std_intp)
##    #plt.show()
##    # Example of handling results
##    print(crossyears)
##    
##        pdftable0=np.zeros((nthresholds,(nyrs-1)*inum+1))
##
##           #compute pdftable, allowing us to take °C slices not just timeslices
##            pdftable0[:,-1]=stats.norm.pdf(heights0,loc=central_est[-1],scale=se[-1])
##            #last one computed first
##            for t in range(nyrs-1):
##                for interp in range(inum):
##                    iloc0=(central_est[t]*(1-interp/inum)+central_est[t]*(interp/inum))
##                    iscale0=(se[t]*(1-interp/inum)+se[t+1]*(interp/inum))
##                    pdftable0[:,t*inum+interp]= stats.norm.pdf(heights0,loc=iloc0,scale=iscale0)
