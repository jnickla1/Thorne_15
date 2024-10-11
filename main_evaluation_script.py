import os
import importlib
import numpy as np
import pandas as pd

def run_methods(years, avg_temperatures, temp_uncert,model_run, experiment_type, methods_folder='Methods'):
    # Find all *_method.py files in the folder and subfolders
    method_files = []
    for root, _, files in os.walk(methods_folder):
        for file in files:
            if file.endswith('_method.py'):
                method_files.append(os.path.join(root, file))

    # Initialize result storage
    results = {}

    # Send data to each method and retrieve results
    for method_path in method_files:
        # Extract method class (folder structure relative to methods_folder)
        method_class = os.path.relpath(os.path.dirname(method_path), methods_folder)

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


import scipy.stats as stats
import statsmodels.stats.multitest as multi
import matplotlib.pyplot as plt

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
    experiment_type = 'historical_obs_20yr_mean'

    # Run all the methods
    results = run_methods(years, temps_obs, (temps_CIl, temps_CIu),model_run, experiment_type)

    heights0=np.arange(.2,4.2,0.05) #heights to test at thresholds: 0.05°C increments
    nthresholds = np.size(heights0)
    inum=4 #internal interpolation within years
    standard = results['cent20y']['LT_trend'][0]
    var_std = np.nanmean(np.abs(np.diff(np.diff(standard))))
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
            

    df_results = pd.DataFrame(columns=['method_name', 'method_class','smooth_ratio','qvals_#yrs<0.5', 'qvals_#yrs<0.1', 'qvals_smallest', 'qvals_smallest5','log-likelihood'])
    i=0
    for method_name, method_data in results.items():
        #print(f"Results from {method_name} (Method Class: {method_data['method_class']}):")
        result = method_data['LT_trend']



        
      
        if isinstance(result, tuple):
            # Central estimate and SE, implying gaussian distn
            central_est, se = result
            var_est = np.nanmean(np.abs(np.diff(np.diff(central_est))))
            pvals = stats.norm.sf(abs((standard-central_est)/ se))*2
            nnans = ~(np.isnan(pvals))
            qvals = multi.multipletests(pvals[nnans], alpha=0.5, method='fdr_bh')
            qvals_count_yrs05 = np.sum(qvals[0])
            qvals_count_yrs01 = np.sum(qvals[1]<0.1)
            qvals_smallest = np.min(qvals[1])
            qvals_smallest5 = np.sort(qvals[1])[4]
            llikelihood = stats.norm.logpdf(standard,loc=central_est,scale=se)
            #print(qvals_count_yrs ,qvals_smallest, qvals_smallest5 )
            
#PLOT TO SHOW WHAT IT'S DOING
            if(method_name == "quartic"):
                plt.plot(years, standard, 'ko')
                plt.fill_between(years, central_est-se, central_est+se)
            
           
        else:
            # 100 percentiles
            percentiles = result
           # print(f"  Percentiles: {percentiles}"
           
        df_results.loc[i]= [ method_name, method_data['method_class'],var_est/var_std, qvals_count_yrs05,qvals_count_yrs01,  qvals_smallest,qvals_smallest5, np.nansum(llikelihood)]
        i=i+1
    df_res_show = df_results
    df_res_show['smooth_ratio'] = df_results['smooth_ratio'].round(3)
    print(df_res_show)
    #df_results.to_csv('method_statistics_results.csv', index=False)
    plt.show()
