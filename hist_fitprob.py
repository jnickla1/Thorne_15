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
from scipy.optimize import minimize
import sys
from netCDF4 import Dataset
from fut_evaluation_gen_ensemble import eval_standard
from datetime import datetime
current_date = datetime.now()
formatted_date = current_date.strftime("%y%m%d")

historical_regen=True #MUST BE TRUE IN THIS FILE

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])


running_subset = ('Methods/42_Temp_Alone', 'Methods/43_Forcing_Based','Methods/44_EarthModel_CGWL')
    #'Methods/42_Temp_Alone/1_Run_Means','Methods/43_Forcing_Based/3_Human_Induced')
    #
                #'Methods/42_Temp_Alone/1_Run_Means','Methods/43_Forcing_Based/1_ERF_FaIR','Methods/43_Forcing_Based/3_Human_Induced',
                 # ,'Methods/43_Forcing_Based/0_Linear','Methods/44_EarthModel_CGWL')


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
sel_methods = [ "CGWL10y_for_halfU","FaIR_nonat","EBMKF_ta2"] #"EBMKF_ta4" "min_month_proj" "OLS_refit_CO2forc", "CGWL10y_for_halfU","TheilSen_h7075" ,"FaIR_anthroA",,"EBMKF_ta2"  ] #"EBMKF_ta",

try:
    index_mapping_new = pd.read_csv('all_methods_statistics_250616.csv')
    def rank2(method_name_in):
        try:
            ret = index_mapping_new[index_mapping_new["method_name"]==method_name_in]["bias50"].values[0] #first is always current
        except:
            print("not found METHOD")
            print(method_name_in)
            ret = 0
        return ret
except:
    print("not found newest DATAFRAME SAVED")
    def rank2(method_name_in):
        return method_name_in
#ftl = np.argsort(index_mapping) #from to list - where a certain method should be plotted


from numpy.polynomial.hermite import hermgauss

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
    experiment_type = sys.argv[1] #'fut_ESM1-2-LR_SSP126_constVolc' #fut_NorESM_RCP45_Volc

    if (experiment_type =='historical'):
        outputfilename = experiment_type
        # Run all the methods

        if regen:
            results = run_methods(years, temps_obs, (temps_CIl, temps_CIu),model_run, experiment_type)
            #no need to pickle

    #To read it back:
        else:
            print("The whole purpose of this python script is to regenerate all of them")
            exit()


#####START FUTURE PROCESSING CODE
    else:
        temps_CIu_past =data.loc[:,"Upper"].to_numpy() #Lower confidence limit (2.5%)	Upper confidence limit (97.5%)
        temps_CIl_past =data.loc[:,"Lower"].to_numpy()
        temps_obs_past=temps_obs #already read in
        years_past=data.loc[:,"Time"].to_numpy()
        
        start_run = int(sys.argv[2])
        outputfilename = f'{experiment_type}{model_run}'
            
        exp_attr = experiment_type.split("_")
        if (exp_attr[1]=='ESM1-2-LR'):
            max_runs = 10+start_run #50  #5
            fut_data_loc = config.CLIMATE_DATA_PATH+'/'+exp_attr[1]+'/combined/'+exp_attr[2].lower()+'_aave_tas.nc'
            hist_data_loc = config.CLIMATE_DATA_PATH+'/'+exp_attr[1]+'/combined/historical_aave_tas.nc'
            
        elif (exp_attr[1]=='NorESM'):
            max_runs =  10+start_run #60
            fut_data_loc = config.CLIMATE_DATA_PATH+'/'+exp_attr[1]+'_volc/BethkeEtAl2017/'+exp_attr[2].lower()+exp_attr[3]+'_aave_tas.nc'
            
            if (exp_attr[3]=='NoVolc'):  #options NoVolc VolcConst Volc
                hist_data_loc = config.CLIMATE_DATA_PATH+'/'+exp_attr[1]+'_volc/BethkeEtAl2017/historicalNoVolc_aave_tas.nc'
            else:
                hist_data_loc = config.CLIMATE_DATA_PATH+'/'+exp_attr[1]+'_volc/BethkeEtAl2017/historicalVolc_aave_tas.nc'
            
        else:
            print("Error: unknown model to this eval script")
            sys.exit(1)
                  
        

            
        dataset = Dataset(fut_data_loc, 'r')
        variable = dataset.variables['tas']
        sims_tas = variable[:].__array__()
        stimes = dataset.variables['time']
        stime_mon = stimes[:].__array__()/365+1850

        dataset_hist = Dataset(hist_data_loc, 'r')
        variable_hist= dataset_hist.variables['tas']
        sims_tas_hist = variable_hist[:].__array__()
        stimes_hist = dataset_hist.variables['time']
        stime_mon_hist = stimes_hist[:].__array__()/365+1850

        print("starting computation for "+experiment_type)

        ens_standard = eval_standard(experiment_type)
            
        model_run =start_run #run only one simulation at a time
        print("Model number:")
        print(model_run)
        print("\n\n\n")
        
        this_sim_yr = average_every_n(sims_tas[model_run,:], 12) #converting monthly to yearly
        this_hsim_yr = average_every_n(sims_tas_hist[model_run,:], 12)
        stime_yrs = np.floor(average_every_n(stime_mon, 12)).astype(int)
        shtime_yrs = np.floor(average_every_n(stime_mon_hist, 12)).astype(int)
        
        start_sim = years_past[-1] - stime_yrs[0] +1 #year that we should switch from observations uncertainty to const simulation uncertainty
        #futCIl = np.full((len(stime_yrs) - start_sim),temps_CIl_past[-1])

        #offset_sim = np.mean( this_sim_yr[0:start_sim] - temps_obs_past[-start_sim:])
        #must decide how to offset the simulation - can do it so the first 50 yrs are 0 as we did for obs

        if (exp_attr[1]=='ESM1-2-LR'):
            offset_sim = np.mean(this_hsim_yr[0:50]) #preindustrial baseline 1850-1899
            sim_corrected = this_sim_yr -offset_sim
            simh_corrected = this_hsim_yr -offset_sim
            simall = np.concatenate((simh_corrected,sim_corrected))
            years= np.concatenate((shtime_yrs,stime_yrs))
            
        elif (exp_attr[1]=='NorESM'):
            #replacing temps_obs_past
            long_past_index = ((gen_orig_number(model_run,np.shape(sims_tas)[0])-1) // 20) #either 1, 2, or 3, still in right order
            long_past_data_loc = config.CLIMATE_DATA_PATH+'/NorESM_volc/NorESM1-M-historical/hist_aave_tas.nc'
            variable = Dataset(long_past_data_loc, 'r').variables['tas']
            long_past_tas_array = variable[:].__array__()
            long_past_tas = average_every_n(long_past_tas_array[long_past_index,:],12)
        
            offset_sim = np.mean( this_hsim_yr[0:20] - temps_obs_past[(-1850+1980):(-1850+1980+20)]) # baseline 1980-2000 match
            sim_corrected = this_sim_yr -offset_sim
            simh_corrected = this_hsim_yr -offset_sim
            long_past_tas_corrected = long_past_tas -offset_sim
        
            simall = np.concatenate((long_past_tas_corrected[0:(-1850+1980)],simh_corrected,sim_corrected))
            fixyrs= 1980-3829
            years= np.concatenate((years_past[0:(-1850+1980)],shtime_yrs+fixyrs,stime_yrs+fixyrs))
            #print(years)
            start_sim= start_sim-fixyrs
            
            
            
        
        temps_CIl_hist = simall[0:len(years_past)]+  temps_CIl_past - temps_obs_past
        temps_CIu_hist = simall[0:len(years_past)]+  temps_CIu_past - temps_obs_past
        futCIl = sim_corrected[start_sim:] - np.mean(temps_obs_past[-10:]- temps_CIl_past[-10:]) #constant small uncertainty
        futCIu = sim_corrected[start_sim:] - np.mean(temps_obs_past[-10:]- temps_CIu_past[-10:])
        # Run all the methods
        #print(len(np.concatenate((temps_CIl_hist,futCIl))))
        #breakpoint()
        results_path = f'Results/results{experiment_type}{model_run}.pickle'

        if regen:
            results = run_methods(years, simall,
                                  (np.concatenate((temps_CIl_hist,futCIl)), np.concatenate((temps_CIu_hist,futCIu)) ),
                                  model_run, experiment_type)
            #no need to pickle
        else:
            print("The whole purpose of this python script is to regenerate all of them")
            exit()      

###END FUTURE PROCESSING
#now still have everything saved in results


    #heights0=np.arange(.2,4.2,0.05) #heights to test at thresholds: 0.05°C increments
    #nthresholds = np.size(heights0)
    inum=12 #internal interpolation within years
    standard = results['cent20y']['LT_trend'][2] #retrospective
    standard_se=results['cent20y']['LT_trend'][3]
    np.save(f'Results2/{outputfilename}_standard.npy', standard)
    np.save(f'Results2/{outputfilename}_standard_se.npy', standard_se)
    
    smooth_std = np.nanmean(np.abs(np.diff(np.diff(standard))))



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

    
    central_yr_estimates =[]

    df_results = pd.DataFrame(columns=['method_name', 'short_method_class','err_var', 'err_var100', 'best_alter_scale'])
    i=0
    all_central_est = []
    all_sampled_lls = []
    all_newsamples = []
    ci=0
    labelcolors=[]
    sorted_results = sorted(results.items(), key=lambda item: (item[1]['method_class'], rank2(item[0])))


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
        print(f"Results from {method_name} (Method Class: {method_data['method_class']}):")
        result = method_data['LT_trend']
        labelcurr_or_retro=""
        isgauss = True #keep track of which kind of method it is
        for k in range(1): #only care about current methods here
            central_est=np.full(len(years),np.nan) #make this blank to start
            if isinstance(result, dict):
                # Empirical method: Call the functions at each time point
                central_est = result['mean'](years,k)
                se = result['se'](years,k)
                pvals = result['pvalue'](years, standard,k)
                llikelihood = result['log_likelihood'](years, standard,k)
                llikelihood2 = stats.norm.logpdf(standard,loc=central_est,scale=se)
                isgauss = False
                if(sum(~np.isnan(central_est))>0):
                    print(f"{method_name} sampled, d llike: {np.nanmean(llikelihood)-np.nanmean(llikelihood2)}")

                def rescale_log_likelihood(scale):
                    deviance = standard - central_est
                    resc_standard = central_est + deviance/scale #larger scale makes deviance appear smaller
                    log_lik=result['log_likelihood'](years, resc_standard,k) - np.log(scale) #pdf expands outward with small scale, so must compress vertically
                    return -np.nansum(log_lik) #minimize this quanity
                    
            else:
                central_est = result[k*2]
                se = result[k*2+1]
                se = np.maximum(se, 0.0005)
                # want to save and store the more accurate calculation if regen: else: load calculations
                pvals = stats.norm.sf(abs((standard-central_est)/ se))*2
                llikelihood = stats.norm.logpdf(standard,loc=central_est,scale=se)
                
                def rescale_log_likelihood(scale_alt):
                    log_lik=stats.norm.logpdf(standard,loc=central_est,scale=se*scale_alt)
                    return -np.nansum(log_lik)

                
            if (sum(~np.isnan(central_est))>0 ):    #isinstance(result, tuple):
                # Central estimate and SE, implying gaussian distn
                labelcurr_or_retro = retIDlabels[k]
                # Perform FDR adjustment and other calculations...
                smooth_est = np.nanmean(np.abs(np.diff(np.diff(central_est))))

                short_method_class = method_data['method_class'][0:2] +"/"+method_data['method_class'].split('/')[-1]

                #CALCULATE things necessary for both combination approaches
                err_var = np.nanmean((central_est-standard)**2)
                err_var100 = np.nanmean((central_est[-100:] -standard[-100:])**2)
                best_alter_scale = minimize(rescale_log_likelihood, x0=1, bounds=[(0.01, 100.0)]).x[0]
                nsamples = 1000
                nnodes = 100
                nyrs = len(years)
                #standard_samples= np.random.normal(loc=standard[:, np.newaxis], scale=standard_se[:, np.newaxis], size=(nyrs,nsamples)) #dimensions
                x, w = hermgauss(nnodes)
                #alpha = w / np.sqrt(np.pi)  # normalized weights - will regenerate in next script
                # Precompute y nodes for the target distribution
                standard_samples = standard[:, None] + standard_se[:, None] * np.sqrt(2.0) * x[None, :] #size=(nyrs,nnodes)

                if(isgauss):
                    sampled_lls = stats.norm.logpdf(standard_samples,loc=central_est[:, np.newaxis],scale=se[:, np.newaxis]*best_alter_scale)
                    newsamples = np.random.normal(loc=central_est[:, np.newaxis], scale=se[:, np.newaxis]*best_alter_scale, size=(nyrs,nsamples))
                else:
                    deviances = standard_samples - central_est[:, np.newaxis] #relative to this method's central estimate
                    resc_standard_samples = central_est[:, np.newaxis] + deviances/best_alter_scale
                    sampled_lls = result['log_likelihood'](years, resc_standard_samples ,k) - np.log(best_alter_scale) #need to ensure dimensions work
                    newsamples = (result['resample'](years, nsamples,k) - central_est[:, np.newaxis])/best_alter_scale + central_est[:, np.newaxis]


                #if(method_name=="CGWL10y_for_halfU"):
                #    breakpoint()
                #if(method_name=="EBMKF_ta2"):
                #    breakpoint()
                df_results.loc[i]= [ method_name,short_method_class,err_var, err_var100, best_alter_scale]
                all_central_est.append(central_est)
                all_sampled_lls.append(sampled_lls)      # shape (nyears,n_samples)
                all_newsamples.append(newsamples)  
                i=i+1
   


    df_results.to_csv(f'Results2/{outputfilename}_names_var_scale.csv', index=True)
    np.save(f'Results2/{outputfilename}_central_est.npy', np.stack(all_central_est))
    np.save(f'Results2/{outputfilename}_hermguass_lls.npy', np.stack(all_sampled_lls))
    np.save(f'Results2/{outputfilename}_newsamples.npy', np.stack(all_newsamples))
    

    print(time.process_time() - start)







