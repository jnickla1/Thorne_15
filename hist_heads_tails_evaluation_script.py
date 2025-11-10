from hist_evaluation_script import *
regen=True
annotate_fig=False
crossing_figs=False


from netCDF4 import Dataset
import sys
#making all paths relative to ~/
from os.path import expanduser
cdataprefix = config.CLIMATE_DATA_PATH+'/'

relthresh=-0.7

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])

def gen_orig_number(new_member_number,sz_ens):
    #fix lexicographic reshuffling
    nums = np.arange(1, sz_ens+1)
    reshuffled = sorted([f"{n}|" for n in nums])
    recovered_order = [int(s.rstrip("|")) for s in reshuffled]
    if new_member_number==-1:
        return recovered_order
    else:
        return recovered_order[new_member_number]

regen = 1 #0 no regen #1 regen completely #2 overwrite regen to allow for computed methods to not need to be redone!
evalmins = np.array([1970,1990, 2005])
evalmaxs = np.array([2000,2024,2024])

data = pd.read_csv("./Common_Data/HadCRUT5.csv")
temps_obs = data.loc[:,"Anomaly"].to_numpy()
preind_base = np.mean(temps_obs[0:50]) #preindustrial baseline 1850-1899
temps_obs_past = temps_obs - preind_base #remove this baseline
temps_CIu_past =data.loc[:,"Upper"].to_numpy() #Lower confidence limit (2.5%)	Upper confidence limit (97.5%)
temps_CIl_past =data.loc[:,"Lower"].to_numpy()
years_past=data.loc[:,"Time"].to_numpy()

def add_dash_dot(sspstr):
    return sspstr[0:4] + "-" + sspstr[4] + "." + sspstr[5:]

def collect_data(exp_attr ):
    if (exp_attr[1]=='rpi'):
        st_array = np.loadtxt('Hist_heads_tails/sst_pseudo_rpi.csv', delimiter=',')
        
    elif (exp_attr[1]=='satcal'):
        st_array = np.loadtxt('Hist_heads_tails/sst_pseudo_satcal.csv', delimiter=',')
        
    else:
        print("Error: unknown model to this eval script "+ exp_attr)
        sys.exit(1)

    return (st_array[:,1:], st_array[:,0])

# "EBMKF_ta4","GAM_AR1",
#"lowess1dg20wnc","Kal_flexLin","FaIR_comb_unB","GWI_tot_CGWL","CGWL10y_sfUKCP" ## "GWI_tot_SR15",


                 
#sel_methods = ["CGWL10y_sfUKCP","FaIR_comb_unB","EBMKF_ta4"]  
sel_methods_list_real = ["Methods/42_Temp_Alone/1_Run_Means/cent20y_method.py",
                             "Methods/42_Temp_Alone/3_ST_Fits/lowess1dg20wnc_method.py",
                             "Methods/42_Temp_Alone/5_Kalman/Kal_flexLin_method.py",
                             "Methods/42_Temp_Alone/4_GAM_AR1/GAM_AR1_method.py",
                             "Methods/44_EarthModel_CGWL/CGWL10y_sfUKCP_method.py",
                             "Methods/43_Forcing_Based/3_Human_Induced/GWI_tot_CGWL_method.py",
                             "Methods/43_Forcing_Based/3_Human_Induced/GWI_tot_SR15_method.py",
                             "Methods/43_Forcing_Based/2_Kalman/EBMKF_ta4_method.py",
                             "Methods/43_Forcing_Based/1_ERF_FaIR/FaIR_comb_unB_method.py" ]

sel_methods_list_anthro = ["Methods/42_Temp_Alone/1_Run_Means/cent20y_method.py",
                           "Methods/43_Forcing_Based/3_Human_Induced/GWI_anthro_CGWL_method.py",
                           "Methods/43_Forcing_Based/3_Human_Induced/GWI_anthro_SR15_method.py",
                           "Methods/43_Forcing_Based/3_Human_Induced/GWI_anthro_method.py",
                           "Methods/43_Forcing_Based/1_ERF_FaIR/FaIR_nonat_method.py"
                           ]



    
def run_one_single_ens_member(plotting_figs, experiment_type, start_run, ax1, ax4, colorraw=None):

    # First evaluation
    #histens_rpi/satcal_real/anthro

    exp_attr = experiment_type.split("_")
    
    sthresh = 0.5 if exp_attr[1]=="rpi" else (relthresh+0.5)
    
    if start_run < 0:
        start_run = -start_run
        max_runs = 1+start_run
        plotting_figs= True
    else:
        max_runs = 20+start_run
        plotting_figs= False
    methods_folder= running_subset

    if exp_attr[2]=="real":
        this_sel_methods_list=sel_methods_list_real
    elif exp_attr[2]=="anthro":
        this_sel_methods_list=sel_methods_list_anthro

    print("starting computation for "+experiment_type)
    print("max_runs" + str(max_runs))


    (sims_tas, stimes) = collect_data(exp_attr)
        
    for model_run in range(start_run,max_runs):
        print("Model number:")
        print(model_run)
        print("\n\n\n")
        

        simall = sims_tas[:,model_run]
        years = stimes.astype(int)

        #internally provided uncertainty from HadCRUT - not relevant to most
        data = pd.read_csv("./Common_Data/HadCRUT5.csv")
        temps_obs = data.loc[:,"Anomaly"].to_numpy()
        temps_CIu =data.loc[:,"Upper"].to_numpy() #Lower confidence limit (2.5%)	Upper confidence limit (97.5%)
        temps_CIl =data.loc[:,"Lower"].to_numpy()
        
        temps_CIl_hist = simall[0:len(years)]+  temps_CIl - temps_obs
        temps_CIu_hist = simall[0:len(years)]+  temps_CIu - temps_obs

        results_path = f'Results/results{experiment_type}{model_run}.pickle'
    
        if regen==1:
            results = run_methods(years, simall,
                                  (temps_CIl_hist, temps_CIu_hist ),
                                  model_run, experiment_type,this_sel_methods_list, completed_methods = set(), give_methods_path = True)
            with open(results_path, 'wb') as fp:
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

        elif regen==2:
            # Load existing results if file exists
            if os.path.exists(results_path):
                with open(results_path, 'rb') as fp:
                    existing_results = pickle.load(fp)
                    #existing_results.pop("EBMKF_ta4",None) #redo EMBKF_ta4
                    #existing_results.pop("lowess1dt30wnc",None)

                    

                completed_methods = set(existing_results.keys())
            else:
                existing_results = {}
                completed_methods = set()
    
        # Run all the incompleted methods
            results = run_methods(years, simall,
                                      (temps_CIl_hist, temps_CIu_hist ),
                                      model_run, experiment_type, methods_folder, completed_methods)
            
            # Process results for saving
            results2 = {}
            for method_name, method_data in results.items():
                result = method_data['LT_trend']
                method_class = method_data['method_class']
                
                if isinstance(result, dict):
                    result2 = (result['mean'](years, 0),result['se'](years, 0),result['mean'](years, 1),result['se'](years, 1))
                else:
                    result2 = result
                
                results2[method_name] = {'method_class': method_class, 'LT_trend': result2}
            
            # Combine old and new results
            combined_results = {**existing_results, **results2}
            
            # Save combined results
            with open(results_path, 'wb') as fp:
                pickle.dump(combined_results, fp)

            results = combined_results #what the code below is expecting
                
    #To just read it back, don't compute anything new:
        else:
            with open (results_path, 'rb') as fp:
                results = pickle.load(fp)
        

        inum=12 #internal interpolation within years

  
        standard = results['cent20y']['LT_trend'][2] #retrospective

        #extend standard 5 yrs into the future
        extlength=6
        m, b = np.polyfit(years[(-15-10):(-10)], standard[(-15-10):(-10)], 1)
        new_years = np.arange(years[-10], years[-10] + extlength)
        new_vals  = m * new_years + b
        standard[(-10):(-10+extlength)] =new_vals
        
        smooth_std = np.nanmean(np.abs(np.diff(np.diff(standard))))
        sfineyr = 1972
        while standard[sfineyr-1850]> sthresh-0.03: #failsafe should any members not be larger than 0.5 in 1972
            sfineyr=sfineyr-1
        fineyrs_all = np.arange(sfineyr,years[-1]+1/inum,1/inum)
        std_intp0 = np.interp(fineyrs_all,years,standard)
        std_intp = std_intp0[~np.isnan(std_intp0)]
        fineyrs_c0 = fineyrs_all[~np.isnan(std_intp0)]
        fineyrs_c = fineyrs_c0[1:]
        thrshs= np.arange(sthresh,np.nanmax(standard) ,.1) #thresholds
        print(f"evaluating {len(thrshs)} of 0.1°C thresholds, starting at 0.5°C")
        closest_years = [-1/inum/2+fineyrs_c[np.logical_and(std_intp[0:-1]<i, std_intp[1:]>=i)][0] for i in thrshs]
                   #will have a variable number of steps, at least including 0.5
        closest_yrs_rnd = np.round(closest_years)

        scalvar = pd.read_csv('Results2/historical_names_var_scale.csv')
        
        central_yr_estimates =[]

        df_results = pd.DataFrame(columns=['method_name', 'method_class','c/r','smooth_r','avg_unc.(1se)','#q<0.5', '#q<0.1', 'q_min', 'q_small5',
                                           'log-likeli','RMS','bias','tlog-l',
                                           '100log-l','100RMS','75RMS','100bias',
                                           'bias50', 'xyear0.5' , 'xyear1.0', 'xyear1.5'
                                           ])
        i=0
        lhund=-100
        ncm = 0 #number of current methods
        for method_name, method_data in results.items():
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

        maxllen = np.max((evalmaxs-evalmins)) #*inum)
        pcrossmatrix = np.full((ncm+1,maxllen, 3), np.nan)  #method, fineyear, [0.5, 1.0, or 1.5 threshold]
        
        for method_name, method_data in results.items():
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

                        
                else:
                    central_est = result[k*2]
                    se = result[k*2+1]
                    # want to save and store the more accurate calculation if regen: else: load calculations
                    pvals = stats.norm.sf(abs((standard-central_est)/ se))*2
                    llikelihood = stats.norm.logpdf(standard,loc=central_est,scale=se)


                if ((sum(~np.isnan(central_est))>0 and k==0) or (method_name=='cent20y' and k==1)):    #isinstance(result, tuple):
                    
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


                    

                    if(labelcurr_or_retro=="c"):
                        #calculate at intermediary 0.?°C

                            
                        method_intp0 = np.interp(fineyrs_all,years,central_est)
                        method_intp = method_intp0[~np.isnan(method_intp0)]
                        fineyrs_m0 = fineyrs_all[~np.isnan(method_intp0)]
                        fineyrs_m = fineyrs_m0[1:]
                        thrshs= np.arange(sthresh,np.nanmax(central_est) ,.1) #threshold
                        try:
                            closest_yearsM = [-1/inum/2+fineyrs_m[np.logical_and(method_intp[0:-1]<i, method_intp[1:]>=i)][0] for i in thrshs]
                        except:
                            print(method_name)
                            print(central_est)
                            breakpoint()
                        #this is only picking the first crossing time if there are multiple ones
                        yearfcr05 = closest_yearsM[0]
                        yearfcr10= closest_yearsM[5]
                        if len(closest_yearsM)<=10:
                            yearfcr15= np.nan
                        else:
                            yearfcr15= closest_yearsM[10]
                        maxcut=0

                    else: #cent20yr standard
                        #extend 5 yrs into the future
                        m, b = np.polyfit(years[(-15-10):(-10)], central_est[(-15-10):(-10)], 1)
                        new_years = np.arange(years[-10], years[-10] + extlength)
                        new_vals  = m * new_years + b
                        central_est[(-10):(-10+extlength)] =new_vals
                        yearfcr05 = closest_years[0]
                        if len(closest_years)<=5:
                            import matplotlib.pyplot as plt
                            plt.figure()
                            plt.plot(years,central_est)
                            plt.show()
                            breakpoint()
                            
                        yearfcr10= closest_years[5]
                        if len(closest_years)<=10:
                            yearfcr15= np.nan
                        else:
                            yearfcr15= closest_years[10]
                        maxcut=10

                    thrshs= np.array([0,0.5,1.0])+sthresh

                    for j in range(3):
                        evalmin=evalmins[j]
                        evalmax=evalmaxs[j]-maxcut
                        evalyrs = np.arange(evalmin,evalmax)
                        fineevalyrs = np.arange(evalyrs[0],evalyrs[-1]+1/inum,1/inum) 
                        this_method_p_steps = np.full(np.shape(evalyrs),np.nan)

                        if method_name == "cent20y":
                            best_alter_scale =1
                        else:
                            best_alter_scale = scalvar[scalvar["method_name"]==method_name]["best_alter_scale"].to_numpy()[0]

                        
                        if isinstance(result, dict):
                            #this_method_p_steps = result['pvalue'](evalyrs, np.full(np.shape(evalyrs),thrshs[j]),k, two_sided=False)
                            deviances = -central_est[evalyrs-1850] + thrshs[j, np.newaxis] #relative to this method's central estimate
                            resc_standard_samples = central_est[evalyrs-1850] + deviances/best_alter_scale
                            this_method_p_steps = result['pvalue'](evalyrs, resc_standard_samples,k, two_sided=False)
                        else:
                            this_method_p_steps = stats.norm.cdf(central_est[evalyrs-1850],loc=thrshs[j, np.newaxis],scale=se[evalyrs-1850]*best_alter_scale)
                        # Replace NaNs at the start with 0s

                        first_non_nan = np.argmax(~np.isnan(this_method_p_steps))
                        this_method_p_steps[:first_non_nan] = 0
                        # Replace NaNs at the end with 1s if they exist
                        last_non_nan = len(this_method_p_steps) - np.argmax(~np.isnan(this_method_p_steps)[::-1]) - 1
                        try:
                            #psteps_intp = np.interp(fineevalyrs,evalyrs,this_method_p_steps)
                            pcrossmatrix[i,0:len(this_method_p_steps),j]=this_method_p_steps
                        except:
                            breakpoint()
                        
                    short_method_class = method_data['method_class'][0:2] +"/"+method_data['method_class'].split('/')[-1]
                    #plt.plot(years,central_est,label=method_name)
                    candidate_row= [ method_name,short_method_class,labelcurr_or_retro,smooth_est/smooth_std,avg_uncert,
                                     qvals_count_yrs05,qvals_count_yrs01,  qvals_smallest,qvals_smallest5, np.nanmean(llikelihood), np.sqrt(np.nanmean((central_est-standard)**2)),
                                       np.nanmean(central_est-standard) , np.nansum(llikelihood),
                                         np.nansum(llikelihood[lhund:-1]),np.sqrt(np.nanmean((central_est[lhund:-1]-standard[lhund:-1])**2)),
                                         np.sqrt(np.nanmean((central_est[(lhund+25):-1]-standard[(lhund+25):-1])**2)),
                                         np.nanmean(central_est[lhund:-1]-standard[lhund:-1]),
                                         np.nanmean(central_est[-50:]-standard[-50:]),
                                         yearfcr05, yearfcr10, yearfcr15
                                    ]
                    
                    df_results.loc[i]=candidate_row
                    i=i+1
       
        print(df_results)
        df_results.to_csv('Results/headstails_statistics_'+experiment_type+str(model_run)+'.csv', index=False,mode='w+' )
        np.save('Results/headstails_method_fineprobs_'+experiment_type+str(model_run)+'.npy', pcrossmatrix) 
        
        print(time.process_time() - start)
        

        print("finished this run") 
    return 0
        



if __name__ == '__main__':
    plotting_figs=False
    experiment_type = sys.argv[1] #'fut_ESM1-2-LR_SSP126_constVolc' #fut_NorESM_RCP45_Volc
    start_run = int(sys.argv[2])
    exp_attr = experiment_type.split("_")
    
    if exp_attr[0]=="histens":
        run_one_single_ens_member(plotting_figs, experiment_type, start_run, None, None)
        #plt.legend()
        #plt.show()    
