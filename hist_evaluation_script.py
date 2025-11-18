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

import colorsys

from datetime import datetime
current_date = datetime.now()
formatted_date = current_date.strftime("%y%m%d")

historical_regen=False #different variable name, whether we are regenerating data (some from saved intermediates or just reading from pickle.

# ================================
# Define Method Full Names, Colors, Location on Violin plots
# ================================
lo=3
methods_names = [
    ['42_Temp_Alone/1_Run_Means', "4.3.1:\nRunning Means", 2, 1998,2022,.5], #98
    ['42_Temp_Alone/2_LT_Fits', "4.3.2:\nLong Term Fits", 7,  1998,2022,.4], #98
    ['42_Temp_Alone/3_ST_Fits', "4.3.3:\nShort Term Fits", 16 ,1998,2022,.9], #98
    ['42_Temp_Alone/4_GAM_AR1', "4.3.4:\nSmoothing Splines", 25 + lo,1997,2022,1.2], #97
    ['42_Temp_Alone/5_Kalman', "4.3.5:\nSimple State Space", 28.5 + lo,1997,2022,1], #97
    ['42_Temp_Alone/6_Remove_IV', "4.3.6:\nFilter Internal\nVariability",33.5+ lo,1999,2022,0.4], #99
    ['43_Forcing_Based/0_Linear', "4.4.1:\nLinear CO2",               37+ lo,1999,2022,1.3],
    ['43_Forcing_Based/1_ERF_FaIR',"4.4.2:\nConstrained\nEmulators",38+4.5+ lo,1999,2022,1.7],  #99
    ['43_Forcing_Based/2_Kalman', "4.4.3:\nEnergy Balance\nKalman Filter",42+5.5+ lo,1999,2022,1.6],  #99
    ['43_Forcing_Based/3_Human_Induced', "4.4.4:\nAttributable\nWarming", 46.5+11+ lo,1999,2022,0.3], #99
    ['44_EarthModel_CGWL', "4.5:\nCombining\nESM Projections", 54+11+ lo,1998,2022,1.5]   #98
    ] #Last three coordinates define location on violin plots

def gen_color(ci, dark=False):
    # used to have implemented a version that would allow a dark version of each of these, now too many colors and gets confusing
    colors = {
    '42_Temp_Alone/1_Run_Means': "#CC6677",  # pink
    '42_Temp_Alone/2_LT_Fits': "#332288",  # blueberry
    '42_Temp_Alone/3_ST_Fits': "#DDCC77",  # Green
    '42_Temp_Alone/4_GAM_AR1': "#117733",  # yellow
    '42_Temp_Alone/5_Kalman': "#88CCEE",  # light blue
    '42_Temp_Alone/6_Remove_IV': "#882255",  # maroon
    '43_Forcing_Based/1_ERF_FaIR':"#44AA99",  # teal
    '43_Forcing_Based/2_Kalman': "#999933",  # mustard
    '43_Forcing_Based/3_Human_Induced': "#AA4499",  #pink
    '43_Forcing_Based/0_Linear': "#9494CE",  #greyblue
    '44_EarthModel_CGWL': "#CC3311"   # firetruck
    }
    
    return colors[ci]


running_subset = ('Methods/42_Temp_Alone/1_Run_Means','Methods/43_Forcing_Based/2_Kalman') #Methods/43_Forcing_Based','Methods/44_EarthModel_CGWL' )#still working on debugging lfca
#running_subset = ('Methods/42_Temp_Alone/','Methods/43_Forcing_Based','Methods/44_EarthModel_CGWL' )

                #'Methods/42_Temp_Alone,'Methods/43_Forcing_Based/1_ERF_FaIR','Methods/43_Forcing_Based/3_Human_Induced',
                 # ,'Methods/43_Forcing_Based/0_Linear','Methods/44_EarthModel_CGWL')
#running_subset = ('Methods/42_Temp_Alone/1_Run_Means', 'Methods/42_Temp_Alone','Methods/43_Forcing_Based','Methods/44_EarthModel_CGWL') #methods we now want to run, smaller subset for debugging

def get_brightness(hex_color):
    """Calculates the perceived brightness of a CSS color, a float between 0 and 1.
    """
    # Convert hex color to RGB values
    r = int(hex_color[1:3], 16)
    g = int(hex_color[3:5], 16)
    b = int(hex_color[5:7], 16)

    # Calculate brightness using the relative luminance formula
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


# ================================
# the list passed to the methods_folder is essential: tells the script which methods to run. Needs >1 location to work.
# Find all *_method.py files in the folder and subfolders
# ================================

def run_methods(years, avg_temperatures, temp_uncert,model_run, experiment_type, methods_folder=running_subset,              
                completed_methods = set(),give_methods_path =False):
    if give_methods_path == False:
        method_files = []
        for root, _, files in chain.from_iterable(os.walk(path) for path in methods_folder):
            for file in files:
                if file.endswith('_method.py'):
                    method_files.append(os.path.join(root, file))
        incomplete_method_files = [f for f in method_files if os.path.basename(f).replace('_method.py', '') not in completed_methods]
    else:
        incomplete_method_files = methods_folder
    # Initialize result storage
    results = {}

    # Send data to each method and retrieve results
    for method_path in incomplete_method_files:
    
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
import collections
import mplcursors
spaglines = []

alt_colors = ['black', 'white']
sel_methods = [ "GWI_tot_CGWL","EBMKF_ta4"]  #"CGWL10y_sfUKCP","FaIR_comb_unB"
##try:
##    index_mapping_new = pd.read_csv('all_methods_statistics_251110True.csv)
##    def rank2(method_name_in):
##        try:
##            ret = index_mapping_new[index_mapping_new["method_name"]==method_name_in]["bias50"].values[0] #first is always current
##        except:
##            print("not found METHOD")
##            print(method_name_in)
##            ret = 0
##        return ret
##except:
##    print("not found newest DATAFRAME SAVED")
##    def rank2(method_name_in):
##        return method_name_in
#ftl = np.argsort(index_mapping) #from to list - where a certain method should be plotted

select_hist=True #plot additional smaller histograms within error
if select_hist:
    avail_methods_list_7 = [ "CGWL10y_for_halfU","EBMKF_ta4","GAM_AR1",
             "lowess1dg20wnc","Kal_flexLin","FaIR_comb_unB","GWI_tot_CGWL",] #"GWI_tot_SR15","CGWL10y_sfUKCP"
    rmse_df = pd.read_csv("current_methods_statistics_251110True.csv")
    avail_methods_list_30 = list(rmse_df.loc[rmse_df['RMS'] <= 0.06, 'method_name'])
  
hist_option='first'


if __name__ == '__main__':
    regen=historical_regen #False #whether to execute all of the methods directly (True) or read from a pickle file to rapidly redraw the figures (False).

    # First evaluation
    data = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50]) #preindustrial baseline 1850-1899
    temps_obs = temps_obs - preind_base #remove this baseline
    temps_CIu =data.loc[:,"Upper"].to_numpy() #Lower confidence limit (2.5%)	Upper confidence limit (97.5%)
    temps_CIl =data.loc[:,"Lower"].to_numpy()
    years=data.loc[:,"Time"].to_numpy()
    nyrs = len(years)
    model_run = ''
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
    

    #heights0=np.arange(.2,4.2,0.05) #heights to test at thresholds: 0.05°C increments
    #nthresholds = np.size(heights0)
    inum=12 #internal interpolation within years
    standard = results['cent20y']['LT_trend'][2] #retrospective
    #breakpoint()
    smooth_std = np.nanmean(np.abs(np.diff(np.diff(standard))))

    fig, (ax1,ax4)= plt.subplots(2, 1, figsize=(10,10), gridspec_kw={ "hspace": 0.3})
    ax1.plot(years, standard, 'k-',zorder=3,lw=1.5)

    #fig2, (ax05, hist_ax05)  = plt.subplots(2, 1, figsize=(7,10),gridspec_kw={'height_ratios': [4, 1] , 'left':0.28 , 'right':0.9 , 'top':0.88, 'bottom':0.12}) #, sharex=True)
    fig2, ax05  = plt.subplots(1, 1, figsize=(7,10),gridspec_kw={ 'left':0.28 , 'right':0.85 , 'top':0.92, 'bottom':0.08}) #, sharex=Tr
    fig2b    , (hist_ax05 ,hist_ax10) = plt.subplots(2, 1, figsize=(5,7),gridspec_kw={ 'left':0.1 , 'right':0.9 , 'top':0.88, 'bottom':0.12, 'hspace':.9}) #, sharex=Tr
    fig3, ax10   = plt.subplots(1, 1, figsize=(7,10),gridspec_kw={'left':0.28 , 'right':0.85 , 'top':0.92, 'bottom':0.08}) #, sharex=True)

    eval05min=1970
    eval05max=2015
    eval05yrs = np.arange(eval05min,eval05max)
    fineyrs = np.arange(eval05min,eval05max+1/inum,1/inum) 
    std_intp0 = np.interp(fineyrs,years,standard)
    std_intp = std_intp0[~np.isnan(std_intp0)]
    closest_year05 = fineyrs[np.argmin(np.abs(std_intp - 0.5))]
    closest_year10 = fineyrs[np.argmin(np.abs(std_intp - 1.0))]
    closest_years = [fineyrs[np.argmin(np.abs(std_intp - i))] for i in np.arange(0.5,1.05,.1)]

    eval10min=1995
    eval10max=2024
    eval10yrs = np.arange(eval10min,eval10max)
    
    eval05 = {}
    eval10 = {}
    crossing_yrs05 =[]
    crossing_yrs10 =[]
    crossing_yrs05_7 =[]
    crossing_yrs10_7 =[]
    crossing_yrs05_30 =[]
    crossing_yrs10_30 =[]
    cmethods = []

    #ax1.hlines(0.5,xmin=eval05min,xmax = eval05max-15, color='black', linestyle='--', linewidth=2 )
    #ax1.hlines(1.0,xmin=eval10min,xmax = eval10max, color='black', linestyle='--', linewidth=2 )
    
    ax05.axvline(closest_year05, color='black', linestyle='--', linewidth=2 ) #, label="Gold Standard (1986)")
    hist_ax05.axvline(closest_year05, color='black', linestyle='--', linewidth=2 ) #, label="Gold Standard (1986)")
    ax10.axvline(closest_year10, color='black', linestyle='--', linewidth=2 ) #, label="Gold Standard (1986)")
    hist_ax10.axvline(closest_year10, color='black', linestyle='--', linewidth=2 ) #, label="Gold Standard (1986)")
    
    central_yr_estimates =[]
    ax4_handles=[]
    ax4_labels=[]
    df_results = pd.DataFrame(columns=['method_name', 'method_class','c/r','smooth_r','avg_unc.(1se)','#q<0.5', '#q<0.1', 'q_min',
                                       'q_small5','log-likeli','RMS','bias','tlog-l','100log-l','l05','l10','bias50','Edyrs2','Edyrs6','cross_err0.5','cross_err1.0','#yrs'])
    i=0
    ci=0
    labelcolors=[]

    import neworder
    sorted_results = neworder.sort_results(results) #sorted(results.items(), key=lambda item: (item[1]['method_class'], rank2(item[0])))


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
        
    
    for method_name, method_data in sorted_results:
        #print(f"Results from {method_name} (Method Class: {method_data['method_class']}):")
        result = method_data['LT_trend']
        labelcurr_or_retro="" #only set if the current method is not blank
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
                    if (method_name!="raw1y"):
                        this_zorder= [row[5] for row in methods_names if row[0] == method_data['method_class']]

                        sline, =ax1.plot(years, central_est, label=method_name, color = gen_color(method_data['method_class'], dark=False),lw=0.5,zorder = float(this_zorder[0]))
                        spaglines.append(sline)
                    #ax4.plot(years, central_est-standard,alpha=0.1, color = gen_color(ci, dark=False))
                    if(method_name in sel_methods): #removeMEI_volc_refit"): "CGWL10y_sUKCP", "KCC_human", "CGWL10y_for_hU","30y_etrend3CS","CGWL_10y_IPCC","lfca_hadcrut"
                        #ax1.fill_between(years, central_est-se, central_est+se, alpha=0.6)
                        #ax1.plot(years, central_est, lw=5,)

                        if (method_name=="FaIR_anthroA"):
                            patch = ax4.fill_between(years, central_est-se-standard, central_est+se-standard, alpha=0.05, color = gen_color(method_data['method_class'], dark=False))
                            #line, = ax4.plot(years, central_est-standard, color = gen_color(ci, dark=True))
                        elif (method_name=="FaIR_nonat"):
                            patch = ax4.fill_between(years, central_est-se-standard, central_est+se-standard, alpha=0.5, color = gen_color(method_data['method_class'], dark=False),zorder=4)
                            line, = ax4.plot(years, central_est-standard, color = gen_color(method_data['method_class'], dark=True))
                        else:
                            #assuming hist_fitprob.py -historical has already been run
                            scalvar = pd.read_csv('Results2/historical_names_var_scale.csv')
                            best_alter_scale = scalvar[scalvar["method_name"]==method_name]["best_alter_scale"].to_numpy()[0]
                            patch = ax4.fill_between(years, central_est-se*best_alter_scale-standard, central_est+se*best_alter_scale-standard, alpha=0.3, color = gen_color(method_data['method_class'], dark=False))
                            line, = ax4.plot(years, central_est-standard, color = gen_color(method_data['method_class'], dark=True))
                            
                        if (method_name!="FaIR_anthroA"):
                            ax4_handles.append((patch,line))
                            ax4_labels.append(method_name)
                        #print(method_name,np.nanmean(central_est-standard))
                        
     #Histogram plots at 0.5°C
                edyrs=0
                aedyrs=0
                cross05v=0
                cross10v=0
                #breakpoint()
                non_decreasing = np.sum([x>y for x, y in zip(central_est[125:-1], central_est[126:])])<= 10
                if(labelcurr_or_retro=="c"):
                    #first compute the p-vals but one-sided
                    this_method_p_steps = np.full(np.shape(eval05yrs),np.nan)
                    if isinstance(result, dict):
                        this_method_p_steps = result['pvalue'](eval05yrs, np.full(np.shape(eval05yrs),0.5),k, two_sided=False)
                    else:
                        this_method_p_steps = stats.norm.cdf(((central_est[eval05yrs-1850]-0.5)/ se[eval05yrs-1850]))

                    # Replace NaNs at the start with 0s
                    first_non_nan = np.argmax(~np.isnan(this_method_p_steps))
                    this_method_p_steps[:first_non_nan] = 0

                    # Replace NaNs at the end with 1s if they exist
                    last_non_nan = len(this_method_p_steps) - np.argmax(~np.isnan(this_method_p_steps)[::-1]) - 1
                    this_method_p_steps[last_non_nan + 1:] = 1
                    #print(method_name)
                    #print(this_method_p_steps)
                    #breakpoint()
                    this_method_crossing_p = np.diff(this_method_p_steps) #probability of crossing in a particular year,
                    halfyrs = eval05yrs[:-1]+0.5   #will plot on half-years instead
                    fineeval05yrs = np.arange(eval05yrs[0],eval05yrs[-1]+1/inum,1/inum) 
                    psteps_intp = np.interp(fineeval05yrs,eval05yrs,this_method_p_steps)
                    integrated_diffs = np.cumsum(psteps_intp) - np.flip(np.cumsum(1-np.flip(psteps_intp)))
                    crossing_start = fineeval05yrs[np.argmax(psteps_intp>0.1587)] #interpolated first moment that this_method_p_steps exceeds 0.16
                    #crossing_exp_value = np.sum(this_method_crossing_p * halfyrs) #cant just do this due to many negatives
                    crossing_exp_value = fineeval05yrs[np.argmin(abs(integrated_diffs))]
                    crossing_end =  fineeval05yrs[len(psteps_intp) -1 - np.argmax(psteps_intp[::-1]<(1-0.1587))]#interpolated last moment that this_method_p_steps exceeds 0.84

                    #crossing_yrs05.append(crossing_exp_value)# do we want a histogram of the integrated average times - no first crossing time
        
                    cmethods.append(method_name)
                    crossing_p_pos = this_method_crossing_p * (this_method_crossing_p>0)
                    crossing_p_pos = crossing_p_pos  / np.sum(crossing_p_pos) #normalize to a total area
                    ax05.fill_between(halfyrs, ci - crossing_p_pos, ci + crossing_p_pos,color=gen_color(method_data['method_class'], dark=False), alpha=0.6, edgecolor='none')
                    labelcolors.append(gen_color(method_data['method_class'], dark=False))
                    xerr_arr = np.array([[max(crossing_exp_value  - crossing_start,0), max(crossing_end -crossing_exp_value,0) ]])
                    ax05.errorbar(crossing_exp_value, ci, xerr=xerr_arr.T, fmt='o', color='black', capsize=3) #alt_colors[ci % 2]
                    eval05[method_name] = [this_method_p_steps, this_method_p_steps,crossing_start,crossing_exp_value , crossing_end] #save for later
                    
                    now_cross = np.logical_and(psteps_intp[0:-1]<0.5, psteps_intp[1:]>=0.5)
                    fineeval05yrsh=fineeval05yrs[0:-1]/2 + fineeval05yrs[1:]/2
                    if(non_decreasing): #method_name in sel_methods):
                        ax05.plot(fineeval05yrsh[now_cross],[ci]*sum(now_cross),"d",markersize=5,markerfacecolor="white", markeredgecolor="black", zorder=5)
                    diffcross = fineeval05yrsh[now_cross]
                    if hist_option == 'first':
                        cross05v = diffcross[0]
                    elif hist_option == 'median':
                        cross05v = crossing_exp_value
                    crossing_yrs05.append(cross05v)
                    if select_hist:
                        if method_name in avail_methods_list_7:
                            crossing_yrs05_7.append(cross05v)
                        if method_name in avail_methods_list_30:
                            crossing_yrs05_30.append(cross05v)
                    edyrs = edyrs + .5*np.mean(np.abs(fineeval05yrsh[now_cross] - closest_year05))
                    

                 #Histogram plots at 1.0°C
                if(labelcurr_or_retro=="c"):
                    #first compute the p-vals but one-sided
                    this_method_p_steps = np.full(np.shape(eval10yrs),np.nan)
                    if isinstance(result, dict):
                        this_method_p_steps = result['pvalue'](eval10yrs, np.full(np.shape(eval10yrs),1.0),k, two_sided=False)
                    else:
                        this_method_p_steps = stats.norm.cdf(((central_est[eval10yrs-1850]-1.0)/ se[eval10yrs-1850]))

                    # Replace NaNs at the start with 0s
                    first_non_nan = np.argmax(~np.isnan(this_method_p_steps))
                    this_method_p_steps[:first_non_nan] = 0

                    # Replace NaNs at the end with 1s if they exist
                    last_non_nan = len(this_method_p_steps) - np.argmax(~np.isnan(this_method_p_steps)[::-1]) - 1
                    this_method_p_steps[last_non_nan + 1:] = 1
                    #print(method_name)
                    #print(this_method_p_steps)
                    #breakpoint()
                    this_method_crossing_p = np.diff(this_method_p_steps) #probability of crossing in a particular year,
                    halfyrs = eval10yrs[:-1]+0.5   #will plot on half-years instead
                    fineeval10yrs = np.arange(eval10yrs[0],eval10yrs[-1]+1/inum,1/inum) 
                    psteps_intp = np.interp(fineeval10yrs,eval10yrs,this_method_p_steps)
                    integrated_diffs = np.cumsum(psteps_intp) - np.flip(np.cumsum(1-np.flip(psteps_intp)))
                    crossing_start = fineeval10yrs[np.argmax(psteps_intp>0.1587)] #interpolated first moment that this_method_p_steps exceeds 0.16
                    #crossing_exp_value = np.sum(this_method_crossing_p * halfyrs) #cant just do this due to many negatives
                    crossing_exp_value = fineeval10yrs[np.argmin(abs(integrated_diffs))]
                    crossing_end =  fineeval10yrs[len(psteps_intp) -1 - np.argmax(psteps_intp[::-1]<(1-0.1587))]#interpolated last moment that this_method_p_steps exceeds 0.84

                    #crossing_yrs10.append(crossing_exp_value) #do we want a histogram of the integrated average times

                    
                    #cmethods.append(method_name) already done
                    crossing_p_pos = this_method_crossing_p * (this_method_crossing_p>0)
                    crossing_p_pos = crossing_p_pos  / np.sum(crossing_p_pos) #normalize to a total area
                    ax10.fill_between(halfyrs, ci - crossing_p_pos, ci + crossing_p_pos,color=gen_color(method_data['method_class'], dark=False), alpha=0.6, edgecolor='none')
                    xerr_arr = np.array([[max(crossing_exp_value  - crossing_start,0), max(crossing_end -crossing_exp_value,0) ]])
                    ax10.errorbar(crossing_exp_value, ci, xerr=xerr_arr.T, fmt='o', color='black', capsize=3) #
                    eval10[method_name] = [this_method_p_steps, this_method_p_steps,crossing_start,crossing_exp_value , crossing_end] #save for later

                    now_cross = np.logical_and(psteps_intp[0:-1]<0.5, psteps_intp[1:]>=0.5)
                    fineeval10yrsh=fineeval10yrs[0:-1]/2 + fineeval10yrs[1:]/2
                    if(non_decreasing):
                        ax10.plot(fineeval10yrsh[now_cross],[ci]*sum(now_cross),"d",markersize=5,markerfacecolor="white", markeredgecolor="black", zorder=5)
                        #plt.figure()
                        #plt.plot(eval10yrs,this_method_p_steps)
                        #plt.errorbar(crossing_exp_value, 0.5, xerr=xerr_arr.T, fmt='o', color=alt_colors[ci % 2], capsize=3)
                        #plt.fill_between(halfyrs, .5 - crossing_p_pos, .5 + crossing_p_pos,color=gen_color(ci, dark=False), alpha=0.6, edgecolor='none')
                    edyrs = edyrs + .5*np.mean(np.abs(fineeval10yrsh[now_cross] - closest_year10))
                    diffcross = fineeval10yrsh[now_cross] 
                    if hist_option == 'first':
                        cross10v = diffcross[0]
                    elif hist_option == 'median':
                        cross10v = crossing_exp_value
                    crossing_yrs10.append(cross10v)
                    if select_hist:
                        if method_name in avail_methods_list_7:
                            crossing_yrs10_7.append(cross10v)
                        if method_name in avail_methods_list_30:
                            crossing_yrs10_30.append(cross10v)
                    
                    ci = ci+1 #increment current method counter
                aedyrs=edyrs*2/len(closest_years)
                #calculate at intermediary 0.°C
                if(labelcurr_or_retro=="c"):
                    thrshs= np.arange(0.5,1.05,.1)
                    for j in np.arange(1,len(closest_years)-1):
                        evalmin=int(closest_years[j])-15
                        evalmax=min(int(closest_years[j])+15,2024)
                        evalyrs = np.arange(evalmin,evalmax)
                        fineevalyrs = np.arange(evalyrs[0],evalyrs[-1]+1/inum,1/inum) 
                        this_method_p_steps = np.full(np.shape(evalyrs),np.nan)
                        if isinstance(result, dict):
                            this_method_p_steps = result['pvalue'](evalyrs, np.full(np.shape(evalyrs),thrshs[j]),k, two_sided=False)
                        else:
                            this_method_p_steps = stats.norm.cdf(((central_est[evalyrs-1850]-thrshs[j])/ se[evalyrs-1850]))
                        # Replace NaNs at the start with 0s
                        first_non_nan = np.argmax(~np.isnan(this_method_p_steps))
                        this_method_p_steps[:first_non_nan] = 0
                        # Replace NaNs at the end with 1s if they exist
                        last_non_nan = len(this_method_p_steps) - np.argmax(~np.isnan(this_method_p_steps)[::-1]) - 1
                        this_method_p_steps[last_non_nan + 1:] = 1
                        psteps_intp = np.interp(fineevalyrs,evalyrs,this_method_p_steps)

                        now_cross = np.logical_and(psteps_intp[0:-1]<0.5, psteps_intp[1:]>=0.5)
                        fineevalyrsh=fineevalyrs[0:-1]/2 + fineevalyrs[1:]/2
                        evalmean = np.nanmean(np.abs(fineevalyrsh[now_cross] - closest_years[j]))
                        if np.isnan(evalmean):
                            evalmean=15
                        aedyrs = aedyrs + 1/len(closest_years)*evalmean
                    
              #  elif("lowess" in method_name):
              #      plt.fill_between(years, central_est-se, central_est+se, alpha=0.6)

                short_method_class = method_data['method_class'][0:2] +"/"+method_data['method_class'].split('/')[-1]
                df_results.loc[i]= [ method_name,short_method_class,labelcurr_or_retro,smooth_est/smooth_std,avg_uncert,
                                 qvals_count_yrs05,qvals_count_yrs01,  qvals_smallest,qvals_smallest5, np.nanmean(llikelihood), np.sqrt(np.nanmean((central_est-standard)**2)),
                                   np.nanmean(central_est-standard) , np.nansum(llikelihood),np.nansum(llikelihood[-100:-1]), np.exp(llikelihood[int(closest_year05)-1850]),
                                     np.exp(llikelihood[int(closest_year10)-1850]),np.nanmean(central_est[-50:]-standard[-50:]),edyrs,aedyrs, cross05v, cross10v, np.count_nonzero(~np.isnan(central_est))] 
                i=i+1
   
       # else:
            # 100 percentiles
        #    percentiles = result
           # print(f"  Percentiles: {percentiles}"
    ax05.set_yticks(range(len(cmethods)))
    cmethods2=np.array(cmethods)
    cmethodsrs=cmethods2 #[index_mapping]
    ax05.set_yticklabels(cmethodsrs , fontsize=8)
    for label, color,lc in zip(ax05.get_yticklabels(), alt_colors * (len(cmethods) // 2 + 1),labelcolors):
        #print(get_brightness(lc))
        c = 1*(get_brightness(lc)<120)
        label.set_color(alt_colors[c])
        label.set_backgroundcolor(lc)
        label.get_bbox_patch().set_boxstyle("round,pad=0.06")

        
    ax10.set_yticks(range(len(cmethods)))
    ax10.set_yticklabels(cmethodsrs , fontsize=8)
    for label, color,lc in zip(ax10.get_yticklabels(), alt_colors * (len(cmethods) // 2 + 1),labelcolors):
        c = 1*(get_brightness(lc)<120)
        label.set_color(alt_colors[c])
        label.set_backgroundcolor(lc)
        label.get_bbox_patch().set_boxstyle("round,pad=0.06")


    hist_ax05.set_xlim([eval05min, eval05max-15])
    ax05.set_xlim([eval05min, eval05max-15])

    hist_ax10.set_xlim([eval10min, 2025])
    ax10.set_xlim([eval10min, 2025])
    ax05.set_ylim([-0.5, ci-.5])
    ax10.set_ylim([-0.5, ci-.5])
    ax05.set_xlabel("Year")
    
    ax05.set_title("Crossing Years for 0.5°C Above Preindustrial by Method", pad=10)
    ax05.set_xticks( hist_ax05.get_xticks())
    ax05.invert_yaxis()
    hist_ax05.hist(crossing_yrs05, color="skyblue", edgecolor='grey',bins =np.arange(eval05min, eval05max-15,1) ,label = f"all {len(crossing_yrs05)} current methods")
    hist_ax05.set_ylabel("Count")
    hist_ax05.set_xlabel("Year")
    hist_ax05.set_yticks( np.arange(0,20,2))
    hist_ax10.set_yticks( np.arange(0,20,2))
    hist_ax05.set_ylim( [0,16.5])
    hist_ax10.set_ylim( [0,16.5])

    ax10.set_xlabel("Year")
    ax10.set_xlim([eval10min, eval10max])
    ax10.set_title("Crossing Years for 1.0°C Above Preindustrial by Method", pad=10)
    ax10.set_xticks( hist_ax10.get_xticks())
    ax10.invert_yaxis()
    hist_ax10.hist(crossing_yrs10, color="skyblue", edgecolor='grey',bins =np.arange(eval10min, eval10max,1))
    hist_ax10.set_ylabel("Count")
    hist_ax10.set_xlabel("Year")
    if hist_option == "first":
        hist_ax05.set_title("Distribution of First 0.5°C Crossing Time")
        hist_ax10.set_title("Distribution of First 1.0°C Crossing Time")
    if hist_option == "median":
        hist_ax10.set_title("Distribution of Median 1.0°C Crossing Time")
        hist_ax05.set_title("Distribution of Median 0.5°C Crossing Time")
    if select_hist:
        hist_ax05.hist(crossing_yrs05_30, color="orange", edgecolor='yellow',bins =np.arange(eval05min, eval05max-15,1), label= f"{len(crossing_yrs05_30)} that passed historical RMSE")
        hist_ax05.hist(crossing_yrs05_7, color="purple",edgecolor='pink',bins =np.arange(eval05min, eval05max-15,1), label= f"final {len(crossing_yrs05_7)} selected (realized warming)")
        hist_ax10.hist(crossing_yrs10_30, color="orange",edgecolor='yellow',bins =np.arange(eval10min, eval10max,1) )
        hist_ax10.hist(crossing_yrs10_7, color="purple",edgecolor='pink',bins =np.arange(eval10min, eval10max,1) )
        hist_ax05.legend(bbox_to_anchor=(0.5, -0.28), loc='upper center')

    df_res_show = df_results.drop(columns=['tlog-l','100log-l'])
    df_res_show['smooth_r'] = df_results['smooth_r'].round(3)
    df_res_cur = df_res_show[df_results['c/r']=='c']
    
    print(df_res_cur.sort_values('log-likeli',ascending=False).head(20))

    ax1.grid(color='silver',zorder=-1)
    ax4.grid(color='silver',zorder=-1)

    #df_results.to_csv('method_statistics_results.csv', index=False)
    #ax1.legend(fontsize=6.5,ncol=4)
    handles, labels = ax1.get_legend_handles_labels()
    numcurmethods = len(labels)

    for ci in range(len(methods_names)):
        ax1.plot([1800,1801], [0,0], label= methods_names[ci][1].split(":")[1].replace('\n',' ').strip(), color = gen_color(methods_names[ci][0], dark=False),lw=0.5)
    
    handles, labels = ax1.get_legend_handles_labels()
    new_handles = handles[numcurmethods:]
    new_labels = labels[numcurmethods:]
    ax1.legend(new_handles, new_labels,ncol=2)
    
    ax4.set_ylim(bottom=-0.12,top=0.12)
    xmin, xmax = [1850,2025] #ax1.get_xlim()
    ax4.set_xlim([xmin, xmax])
    ax4.set_xticks(np.arange(1850,2026,25))
    ax4.set_yticks(np.arange(-0.1,0.11,0.025))
    ax1.set_xticks(np.arange(1850,2026,25))
    ax1.set_xlim([xmin, xmax])
    line_20yr = ax4.hlines(y=0,xmin=years[np.argmax(~np.isnan(standard))],
            xmax=years[-np.argmax(~np.isnan(np.flip(standard)))],  color='k', linestyle='-', lw=3)
    patch_20yr=ax4.fill_between(years,-results['cent20y']['LT_trend'][3],results['cent20y']['LT_trend'][3],facecolor='lightgray',edgecolor='black',alpha=0.35)
    ax4_handles.insert(0,(patch_20yr,line_20yr))
    ax4_labels.insert(0,"20-yr centred running \n    mean (HadCRUT5)")
    ax4.legend(ax4_handles,ax4_labels,loc='lower left')
    ax1.set_xlabel("Year")
    ax1.set_title("Evaluated Methods to find Current Long-Term Temperature", pad=10)
    ax4.set_title("`Error` of Top-Performing Methods\n (Method) — 20-yr running mean", pad=10)
    ax4.set_xlabel("Year")
    ax1.set_ylabel("Temperature (°C) Anomaly\n relative to 1850-1900")
    ax4.set_ylabel("Temperature (°C) Difference\n relative to 20-yr running mean")



    hai = [2000,2025]

    for i,ax in enumerate([ax05,ax10]):
        ax.spines['right'].set_visible(False)
        for m in methods_names:
            classcolor = gen_color(m[0])
            ax.text(hai[i],m[2],m[1],fontsize=12,horizontalalignment='center',color=classcolor,fontweight='bold',linespacing = 0.8)



    cai = [1985,2010]
    caiht = np.array([16,16])*2/3 #[12,20]
    
    for i,ax in enumerate([hist_ax05,hist_ax10]):
        ax.arrow(cai[i]-8, caiht[i], -4, 0,head_width=caiht[i]*3/12,head_length=2,facecolor='white',lw=2)
        ax.text(cai[i]-8, caiht[i]*.8,"Too Early", fontsize=12,horizontalalignment='center',color='black')
        ax.arrow(cai[i]+8, caiht[i], 4, 0,head_width=caiht[i]*3/12,head_length=2,facecolor='white',lw=2)
        ax.text(cai[i]+8, caiht[i]*.8,"Too Late", fontsize=12,horizontalalignment='center',color='black')

    curcolor =  gen_color('42_Temp_Alone/2_LT_Fits')
    ax1.annotate('OLS_refit',
            xy=(2002,0.65), xycoords='data',
            xytext=(2014, .78-.35), textcoords='data', color =curcolor,
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
            horizontalalignment='center', verticalalignment='bottom')
    ax1.annotate('quartic',
            xy=(1976,.115), xycoords='data',
            xytext=(1985, -.15), textcoords='data',color=curcolor,
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
            horizontalalignment='center', verticalalignment='bottom')
    
    curcolor =  gen_color('42_Temp_Alone/1_Run_Means')
    ax1.annotate('lag10y',
            xy=(2018.5,1.06), xycoords='data',
            xytext=(2018, .7), textcoords='data',
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),color="white",
            horizontalalignment='center', verticalalignment='bottom')
    ax1.annotate('lag10y',
            xy=(2002,0.72), xycoords='data',
            xytext=(2018, .7), textcoords='data',color=curcolor,
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
            horizontalalignment='center', verticalalignment='bottom')
    curcolor =  gen_color('42_Temp_Alone/3_ST_Fits')
    ax1.annotate('11y_offset',
            xy=(1887,0.17), xycoords='data',
            xytext=(1890, .35), textcoords='data',color=curcolor,
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
            horizontalalignment='center', verticalalignment='bottom')
    ax1.annotate('min_month_proj',
            xy=(1878,0.3), xycoords='data',
            xytext=(1878, .5), textcoords='data',color=curcolor,
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
            horizontalalignment='center', verticalalignment='bottom')
    curcolor =  gen_color('43_Forcing_Based/3_Human_Induced')
    ax1.annotate('KCC_all',
            xy=(1911,.16), xycoords='data',
            xytext=(1912, .35), textcoords='data',color=curcolor,
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
            horizontalalignment='center', verticalalignment='bottom')
    curcolor =  gen_color('43_Forcing_Based/0_Linear')
    ax1.annotate('cons_cArrh_CO2forc',
            xy=(1925,.22), xycoords='data',
            xytext=(1927, .52), textcoords='data',color=curcolor,
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
            horizontalalignment='center', verticalalignment='bottom')
    curcolor =  gen_color('43_Forcing_Based/1_ERF_FaIR')
##    ax1.annotate('FaIR_anthro',
##            xy=(1946,.1), xycoords='data',
##            xytext=(1960, -.175), textcoords='data',color=curcolor,
##            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
##            horizontalalignment='center', verticalalignment='bottom')
##    ax1.annotate('FaIR_all',
##            xy=(1994,.338), xycoords='data',
##            xytext=(2004, .128), textcoords='data',color=curcolor,
##            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
##            horizontalalignment='center', verticalalignment='bottom')
    curcolor =  gen_color('42_Temp_Alone/6_Remove_IV')
    ax1.annotate('remove_MEI\nvolc_refit',
            xy=(2022.4,1.585), xycoords='data',
            xytext=(1990,1.13), textcoords='data',
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),color="white",
            horizontalalignment='center', verticalalignment='bottom')
    ax1.annotate('remove_MEI\nvolc_refit',
            xy=(1990,.855), xycoords='data',
            xytext=(1990,1.13), textcoords='data',color=curcolor,
            arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
            horizontalalignment='center', verticalalignment='bottom')





    
    

    if(True):
        # Adding interactive tooltips
        cursor = mplcursors.cursor(spaglines, hover=True)

        # Customize the displayed tooltip
        @cursor.connect("add")
        def on_add(sel):
            # Show the label of the hovered line
            sel.annotation.set_text(sel.artist.get_label())
    
    print(time.process_time() - start)
    plt.show()
    df_res_show2 = df_results.copy()
    df_res_show2['smooth_r'] = df_results['smooth_r'].round(3)
    df_res_cur2 = df_res_show2[df_results['c/r']=='c']
    dfres2 = df_res_cur2.sort_values('log-likeli',ascending=False)
    dfres2.to_csv('current_methods_statistics_'+formatted_date+str(historical_regen)+'.csv', index=False)
    dfres3 = df_res_cur2.sort_values('RMS',ascending=True)
    dfres3.rename(columns={'cross_err0.5': 'cross_year0.5','cross_err1.0': 'cross_year1.0'}, inplace=True)
    columns_to_keep = ["method_name","method_class","#q<0.5","log-likeli","RMS","bias",'cross_year0.5','cross_year1.0','#yrs']
    def increment_second_digit(s):
        if len(s) >= 2 and s[1].isdigit():
            # Convert the second character to an integer, add 1, and convert back to string
            new_digit = str(int(s[1]) + 1)
            # Reconstruct the string with the updated second digit
            return s[0] + new_digit + s[2:]
        return s # Return original string if not applicable
    dfres3["method_class"] = dfres3["method_class"].apply(increment_second_digit)
    dfres3[columns_to_keep].to_csv('current_methods_statistics_cropped'+formatted_date+str(historical_regen)+'.csv', index=False)
    df_results.to_csv('all_methods_statistics_'+formatted_date+str(historical_regen)+'.csv', index=False)
    #sorted_df = df_res_cur2.reset_index(drop=True).sort_values(by=['method_class', 'bias50']).reset_index()
    #sorted_df[['index']].to_csv('to_index_mapping.csv', index=False)


    fig.savefig('spaggheti_plot_v6.png', dpi=500, bbox_inches='tight')
    fig2.savefig('05_threshold_v6.png', dpi=500, bbox_inches='tight')
    fig2b.savefig(f'crossing_pile_v6_{hist_option}.png', dpi=500, bbox_inches='tight')
    fig3.savefig('10_threshold_v6.png', dpi=500, bbox_inches='tight')
    print(sum((crossing_yrs05_30 -closest_year05)<=0.5  ) / len(crossing_yrs05_30))
    print(sum((crossing_yrs05_30 -closest_year05)<=1  ) / len(crossing_yrs05_30))
    print(sum((crossing_yrs05_30 -closest_year05)<=2  ) / len(crossing_yrs05_30))

    print(sum((crossing_yrs10_30 -closest_year10)<=0.5  ) / len(crossing_yrs10_30))
    print(sum((crossing_yrs10_30 -closest_year10)<=1  ) / len(crossing_yrs10_30))
    print(sum((crossing_yrs10_30 -closest_year10)<=2  ) / len(crossing_yrs10_30))
    


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
