from hist_evaluation_script import *
regen=True
annotate_fig=False
crossing_figs=False
sel_methods = ["CGWL10y_for_halfU","TheilSen_h7075" ,"FaIR_anthroA","FaIR_anthroA2","EBMKF_ta3"  ]
from netCDF4 import Dataset
import sys

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])


regen = 1 #2 overwrite regen to allow for computed methods to not need to be redone!

if __name__ == '__main__':
    # First evaluation
    data = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50]) #preindustrial baseline 1850-1899
    temps_obs_past = temps_obs - preind_base #remove this baseline
    temps_CIu_past =data.loc[:,"Upper"].to_numpy() #Lower confidence limit (2.5%)	Upper confidence limit (97.5%)
    temps_CIl_past =data.loc[:,"Lower"].to_numpy()
    years_past=data.loc[:,"Time"].to_numpy()
    
    
    experiment_type = sys.argv[1] #'fut_ESM1-2-LR_SSP126_constVolc' #fut_NorESM_RCP45_Volc
    exp_attr = experiment_type.split("_")
    if (exp_attr[1]=='ESM1-2-LR'):
        max_runs = 2#50  #5
        fut_data_loc = '/Users/JohnMatthew/climate_data/'+exp_attr[1]+'/combined/'+exp_attr[2].lower()+'_aave_tas.nc'
        hist_data_loc = '/Users/JohnMatthew/climate_data/'+exp_attr[1]+'/combined/historical_aave_tas.nc'
        
    elif (exp_attr[1]=='NorESM'):
        max_runs = 60
        fut_data_loc = '/Users/JohnMatthew/climate_data/'+exp_attr[1]+'_volc/BethkeEtAl2017/'+exp_attr[2].lower()+exp_attr[3]+'_aave_tas.nc'
        
        if (exp_attr[3]=='NoVolc'):  #options NoVolc VolcConst Volc
            hist_data_loc = '/Users/JohnMatthew/climate_data/'+exp_attr[1]+'_volc/BethkeEtAl2017/historicalNoVolc_aave_tas.nc'
        else:
            hist_data_loc = '/Users/JohnMatthew/climate_data/'+exp_attr[1]+'_volc/BethkeEtAl2017/historicalVolc_aave_tas.nc'
        
    else:
        print("Error: unknown model to this eval script")
        sys.exit(1)
              
    
    
    methods_folder=('Methods/42_Temp_Alone/1_Run_Means', 'Methods/43_Forcing_Based/2_Kalman')

        #'Methods/42_Temp_Alone/1_Run_Means',,'Methods/42_Temp_Alone/3_ST_Fits',  'Methods/42_Temp_Alone/4_GAM_AR1')#

        #'Methods/42_Temp_Alone/1_Run_Means','Methods/42_Temp_Alone/2_LT_Fits','Methods/42_Temp_Alone/3_ST_Fits',
         #           'Methods/42_Temp_Alone/4_GAM_AR1','Methods/43_Forcing_Based/2_Kalman') #,'Methods/43_Forcing_Based','Methods/44_EarthModel_CGWL')
        #big problems with future ENSO index (can't compute MEI) and need to get Bjorn Samset data


        
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
        
    for model_run in range(1,max_runs):
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
            offset_sim = np.mean( this_hsim_yr[0:20] - temps_obs_past[(-1850+1980):(-1850+1980+20)]) # baseline 1980-2000 match
            sim_corrected = this_sim_yr -offset_sim
            simh_corrected = this_hsim_yr -offset_sim
            simall = np.concatenate((temps_obs_past[0:(-1850+1980)],simh_corrected,sim_corrected))
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
    
        if regen==1:
            results = run_methods(years, simall,
                                  (np.concatenate((temps_CIl_hist,futCIl)), np.concatenate((temps_CIu_hist,futCIu)) ),
                                  model_run, experiment_type, methods_folder)
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
                completed_methods = set(existing_results.keys())
            else:
                existing_results = {}
                completed_methods = set()
    
        # Run all the incompleted methods
            results = run_methods(years, simall,
                                      (np.concatenate((temps_CIl_hist,futCIl)), np.concatenate((temps_CIu_hist,futCIu)) ),
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
        

        heights0=np.arange(.2,4.2,0.05) #heights to test at thresholds: 0.05°C increments
        nthresholds = np.size(heights0)
        inum=12 #internal interpolation within years
        standard = results['cent20y']['LT_trend'][2] #retrospective
        smooth_std = np.nanmean(np.abs(np.diff(np.diff(standard))))

        fig, (ax1,ax4)= plt.subplots(2, 1, figsize=(10,10), gridspec_kw={ "hspace": 0.3})
        ax1.plot(years, standard, 'k-',zorder=3,lw=1.5)


        if(crossing_figs):
            #fig2, (ax05, hist_ax05)  = plt.subplots(2, 1, figsize=(7,10),gridspec_kw={'height_ratios': [4, 1] , 'left':0.28 , 'right':0.9 , 'top':0.88, 'bottom':0.12}) #, sharex=True)
            fig2, ax05  = plt.subplots(1, 1, figsize=(7,10),gridspec_kw={ 'left':0.28 , 'right':0.85 , 'top':0.92, 'bottom':0.08}) #, sharex=Tr
            fig2b    , (hist_ax05 ,hist_ax10) = plt.subplots(2, 1, figsize=(5,5),gridspec_kw={ 'left':0.1 , 'right':0.9 , 'top':0.88, 'bottom':0.12, 'hspace':0.5}) #, sharex=Tr
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
        df_results = pd.DataFrame(columns=['method_name', 'method_class','c/r','smooth_r','avg_unc.(1se)','#q<0.5', '#q<0.1', 'q_min', 'q_small5','log-likeli','RMS','bias','tlog-l','100log-l','l05','l10','bias50','Edyrs2','Edyrs6'])
        i=0
        ci=0
        labelcolors=[]
        sorted_results = sorted(results.items(), key=lambda item: item[1]['method_class'])


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
        if (len(index_mapping) != ncm): #not computing all methods, for debugging only
            index_mapping = np.arange(ncm)
            ftl = np.argsort(index_mapping)
            
        
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
                        if (method_name!="raw01y"):
                            
                            sline, =ax1.plot(years, central_est, label=method_name, color = gen_color(method_data['method_class'], dark=False),lw=0.5)
                            spaglines.append(sline)
                        #ax4.plot(years, central_est-standard,alpha=0.1, color = gen_color(ci, dark=False))
                        if(method_name in sel_methods): #removeMEI_volc_refit"): "CGWL10y_sUKCP", "KCC_human", "CGWL10y_for_hU","30y_etrend3CS","CGWL_10y_IPCC","lfca_hadcrut"
                            #ax1.fill_between(years, central_est-se, central_est+se, alpha=0.6)
                            #ax1.plot(years, central_est, lw=5,)

                            if (method_name=="FaIR_anthroA"):
                                patch = ax4.fill_between(years, central_est-se-standard, central_est+se-standard, alpha=0.05, color = gen_color(method_data['method_class'], dark=False))
                                #line, = ax4.plot(years, central_est-standard, color = gen_color(ci, dark=True))
                            elif (method_name=="FaIR_anthroA2"):
                                patch = ax4.fill_between(years, central_est-se-standard, central_est+se-standard, alpha=0.3, color = gen_color(method_data['method_class'], dark=False),zorder=4)
                                line, = ax4.plot(years, central_est-standard, color = gen_color(method_data['method_class'], dark=True))
                            else:
                                patch = ax4.fill_between(years, central_est-se-standard, central_est+se-standard, alpha=0.3, color = gen_color(method_data['method_class'], dark=False))
                                line, = ax4.plot(years, central_est-standard, color = gen_color(method_data['method_class'], dark=True))
                                
                            if (method_name!="FaIR_anthroA"):
                                ax4_handles.append((patch,line))
                                ax4_labels.append(method_name)
                            #print(method_name,np.nanmean(central_est-standard))
                    if(crossing_figs):       
             #Histogram plots at 0.5°C
                        edyrs=0
                        aedyrs=0
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

                            crossing_yrs05.append(crossing_exp_value)
                            cmethods.append(method_name)
                            crossing_p_pos = this_method_crossing_p * (this_method_crossing_p>0)
                            crossing_p_pos = crossing_p_pos  / np.sum(crossing_p_pos) #normalize to a total area
                            ax05.fill_between(halfyrs, ftl[ci] - crossing_p_pos, ftl[ci] + crossing_p_pos,color=gen_color(method_data['method_class'], dark=False), alpha=0.6, edgecolor='none')
                            labelcolors.append(gen_color(method_data['method_class'], dark=False))
                            xerr_arr = np.array([[max(crossing_exp_value  - crossing_start,0), max(crossing_end -crossing_exp_value,0) ]])
                            ax05.errorbar(crossing_exp_value, ftl[ci], xerr=xerr_arr.T, fmt='o', color='black', capsize=3) #alt_colors[ci % 2]
                            eval05[method_name] = [this_method_p_steps, this_method_p_steps,crossing_start,crossing_exp_value , crossing_end] #save for later
                            
                            now_cross = np.logical_and(psteps_intp[0:-1]<0.5, psteps_intp[1:]>=0.5)
                            fineeval05yrsh=fineeval05yrs[0:-1]/2 + fineeval05yrs[1:]/2
                            if(method_name in sel_methods):
                                ax05.plot(fineeval05yrsh[now_cross],ftl[ci],"d",markersize=5,markerfacecolor="white", markeredgecolor="black", zorder=5)
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

                            crossing_yrs10.append(crossing_exp_value)
                            #cmethods.append(method_name) already done
                            crossing_p_pos = this_method_crossing_p * (this_method_crossing_p>0)
                            crossing_p_pos = crossing_p_pos  / np.sum(crossing_p_pos) #normalize to a total area
                            ax10.fill_between(halfyrs, ftl[ci] - crossing_p_pos, ftl[ci] + crossing_p_pos,color=gen_color(method_data['method_class'], dark=False), alpha=0.6, edgecolor='none')
                            xerr_arr = np.array([[max(crossing_exp_value  - crossing_start,0), max(crossing_end -crossing_exp_value,0) ]])
                            ax10.errorbar(crossing_exp_value, ftl[ci], xerr=xerr_arr.T, fmt='o', color='black', capsize=3) #
                            eval10[method_name] = [this_method_p_steps, this_method_p_steps,crossing_start,crossing_exp_value , crossing_end] #save for later

                            now_cross = np.logical_and(psteps_intp[0:-1]<0.5, psteps_intp[1:]>=0.5)
                            fineeval10yrsh=fineeval10yrs[0:-1]/2 + fineeval10yrs[1:]/2
                            if(method_name in sel_methods):
                                ax10.plot(fineeval10yrsh[now_cross],ftl[ci],"d",markersize=5,markerfacecolor="white", markeredgecolor="black", zorder=5)
                                #plt.figure()
                                #plt.plot(eval10yrs,this_method_p_steps)
                                #plt.errorbar(crossing_exp_value, 0.5, xerr=xerr_arr.T, fmt='o', color=alt_colors[ci % 2], capsize=3)
                                #plt.fill_between(halfyrs, .5 - crossing_p_pos, .5 + crossing_p_pos,color=gen_color(ci, dark=False), alpha=0.6, edgecolor='none')
                            edyrs = edyrs + .5*np.mean(np.abs(fineeval10yrsh[now_cross] - closest_year10))
                            
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
                    else: #blank values
                        closest_year05=1985
                        closest_year10=2010
                        edyrs=0
                        aedyrs=0
                    short_method_class = method_data['method_class'][0:2] +"/"+method_data['method_class'].split('/')[-1]
                    df_results.loc[i]= [ method_name,short_method_class,labelcurr_or_retro,smooth_est/smooth_std,avg_uncert,
                                     qvals_count_yrs05,qvals_count_yrs01,  qvals_smallest,qvals_smallest5, np.nanmean(llikelihood), np.sqrt(np.nanmean((central_est-standard)**2)),
                                       np.nanmean(central_est-standard) , np.nansum(llikelihood),np.nansum(llikelihood[-100:-1]), np.exp(llikelihood[int(closest_year05)-1850]),
                                         np.exp(llikelihood[int(closest_year10)-1850]),np.nanmean(central_est[-50:]-standard[-50:]),edyrs,aedyrs] 
                    i=i+1
       
        if(crossing_figs): 
            ax05.set_yticks(range(len(cmethods)))
            cmethods2=np.array(cmethods)
            cmethodsrs=cmethods2[index_mapping]
            ax05.set_yticklabels(cmethodsrs , fontsize=8)
            for label, color,lc in zip(ax05.get_yticklabels(), alt_colors * (len(cmethods) // 2 + 1),labelcolors):
                #print(get_brightness(lc))
                c = 1*(get_brightness(lc)<120)
                label.set_color(alt_colors[c])
                label.set_backgroundcolor(lc)
                label.get_bbox_patch().set_boxstyle("round,pad=0.22")

                
            ax10.set_yticks(range(len(cmethods)))
            ax10.set_yticklabels(cmethodsrs , fontsize=8)
            for label, color,lc in zip(ax10.get_yticklabels(), alt_colors * (len(cmethods) // 2 + 1),labelcolors):
                c = 1*(get_brightness(lc)<120)
                label.set_color(alt_colors[c])
                label.set_backgroundcolor(lc)
                label.get_bbox_patch().set_boxstyle("round,pad=0.22")


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
            hist_ax05.hist(crossing_yrs05, color="skyblue", edgecolor='grey',bins =np.arange(eval05min, eval05max-15,1) )
            hist_ax05.set_ylabel("Count")
            hist_ax05.set_xlabel("Year")
            hist_ax05.set_yticks( np.arange(0,12,2))
            hist_ax10.set_yticks( np.arange(0,20,4))

            ax10.set_xlabel("Year")
            ax10.set_xlim([eval10min, eval10max])
            ax10.set_title("Crossing Years for 1.0°C Above Preindustrial by Method", pad=10)
            ax10.set_xticks( hist_ax10.get_xticks())
            ax10.invert_yaxis()
            hist_ax10.hist(crossing_yrs10, color="skyblue", edgecolor='grey',bins =np.arange(eval10min, eval10max,1))
            hist_ax10.set_ylabel("Count")
            hist_ax10.set_xlabel("Year")
            hist_ax05.set_title("Distribution of Median 0.5°C Crossing Time for All Methods")
            hist_ax10.set_title("Distribution of Median 1.0°C Crossing Time for All Methods")

           

            hai = [2000,2025] 

            for i,ax in enumerate([ax05,ax10]):
                ax.spines['right'].set_visible(False)
                for m in methods_names:
                    classcolor = gen_color(m[0])
                    ax.text(hai[i],m[2],m[1],fontsize=12,horizontalalignment='center',color=classcolor,fontweight='bold')



            cai = [1985,2010]
            caiht = np.array([12,20])*2/3
            
            for i,ax in enumerate([hist_ax05,hist_ax10]):
                ax.arrow(cai[i]-8, caiht[i], -4, 0,head_width=caiht[i]*3/12,head_length=2,facecolor='white',lw=2)
                ax.text(cai[i]-8, caiht[i]*.8,"Too Early", fontsize=12,horizontalalignment='center',color='black')
                ax.arrow(cai[i]+8, caiht[i], 4, 0,head_width=caiht[i]*3/12,head_length=2,facecolor='white',lw=2)
                ax.text(cai[i]+8, caiht[i]*.8,"Too Late", fontsize=12,horizontalalignment='center',color='black')
                
        df_res_show = df_results.drop(columns=['tlog-l','100log-l'])
        df_res_show['smooth_r'] = df_results['smooth_r'].round(3)
        df_res_cur = df_res_show[df_results['c/r']=='c']
        
        print(df_res_cur.sort_values('log-likeli',ascending=False))

        ax1.grid(color='silver',zorder=-1)
        ax4.grid(color='silver',zorder=-1)

        #df_results.to_csv('method_statistics_results.csv', index=False)
        ax1.legend(fontsize=6.5,ncol=4)
        ax4.set_ylim(bottom=-0.11,top=0.11)
        xmin, xmax = [1850,2100] #ax1.get_xlim() #2000
        ax4.set_xlim([xmin, xmax])
        ax4.set_xticks(np.arange(2000,2105,25))
        ax1.set_xticks(np.arange(2000,2105,25))
        ax1.set_xlim([xmin, xmax])
        ax4.hlines(y=0,xmin=years[np.argmax(~np.isnan(standard))],xmax=years[-np.argmax(~np.isnan(np.flip(standard)))],  color='k', linestyle='-', lw=3)
        ax4.legend(ax4_handles,ax4_labels)
        ax1.set_xlabel("Year")
        ax1.set_title("Evaluated Methods to find Current "+exp_attr[1]+" "+exp_attr[2], pad=10)
        ax4.set_title("`Error` of Top-Performing Methods\n versus 20-yr running mean", pad=10)
        ax4.set_xlabel("Year")
        ax1.set_ylabel("Temperature (°C) Anomaly\n relative to 1850-1900")
        ax4.set_ylabel("Temperature (°C) Difference\n relative to 20-yr running mean")


        if(annotate_fig): 
            curcolor =  gen_color('42_Temp_Alone/2_LT_Fits')
            ax1.annotate('OLS_refit',
                    xy=(2002,0.65), xycoords='data',
                    xytext=(2018, .78-.35), textcoords='data', color =curcolor,
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
                    xy=(2005,0.76), xycoords='data',
                    xytext=(2018, .7), textcoords='data',color=curcolor,
                    arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
                    horizontalalignment='center', verticalalignment='bottom')
            curcolor =  gen_color('42_Temp_Alone/3_ST_Fits')
            ax1.annotate('11y_offset',
                    xy=(1885,0.17), xycoords='data',
                    xytext=(1887, .35), textcoords='data',color=curcolor,
                    arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
                    horizontalalignment='center', verticalalignment='bottom')
            curcolor =  gen_color('43_Forcing_Based/3_Human_Induced')
            ax1.annotate('KCC_all',
                    xy=(1920,.22), xycoords='data',
                    xytext=(1922, .4), textcoords='data',color=curcolor,
                    arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
                    horizontalalignment='center', verticalalignment='bottom')
            curcolor =  gen_color('43_Forcing_Based/1_ERF_FaIR')
            ax1.annotate('FaIR_anthro',
                    xy=(1946,.1), xycoords='data',
                    xytext=(1960, -.175), textcoords='data',color=curcolor,
                    arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
                    horizontalalignment='center', verticalalignment='bottom')
            ax1.annotate('FaIR_all',
                    xy=(1994,.338), xycoords='data',
                    xytext=(2004, .128), textcoords='data',color=curcolor,
                    arrowprops=dict(facecolor=curcolor,shrink=0,width=1.5,headwidth=5),
                    horizontalalignment='center', verticalalignment='bottom')
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
        dfres2.to_csv('Results/current_fut_statistics_'+experiment_type+str(model_run)+'.csv', index=False)
        df_results.to_csv('Results/all_fut_statistics_'+experiment_type+str(model_run)+'.csv', index=False)
        
        if(False):
            #don't change the order on this future set of experiments yet
            sorted_df = df_res_cur2.reset_index(drop=True).sort_values(by=['method_class', 'bias50']).reset_index()
            sorted_df[['index']].to_csv('to_index_mapping.csv', index=False)




