from hist_evaluation_script import *
regen=True
annotate_fig=False
crossing_figs=False
sel_methods = ["CGWL10y_sfUKCP","FaIR_comb_unB","EBMKF_ta4"]
from netCDF4 import Dataset
import sys
from fut_evaluation_gen_ensemble import eval_standard
#making all paths relative to ~/
from os.path import expanduser
cdataprefix = config.CLIMATE_DATA_PATH+'/'

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
    if (exp_attr[1]=='ESM1-2-LR'):
       # max_runs = 10+start_run #50  #5
        fut_data_loc = cdataprefix +exp_attr[1]+'/combined/'+exp_attr[2].lower()+'_aave_tas.nc'
        hist_data_loc =cdataprefix +exp_attr[1]+'/combined/historical_aave_tas.nc'
        
    elif (exp_attr[1]=='NorESM'):
       # max_runs =  10+start_run #60
        fut_data_loc = cdataprefix+exp_attr[1]+'_volc/BethkeEtAl2017/'+exp_attr[2].lower()+exp_attr[3]+'_aave_tas.nc'
        
        if (exp_attr[3]=='NoVolc'):  #options NoVolc VolcConst Volc
            hist_data_loc = cdataprefix+exp_attr[1]+'_volc/BethkeEtAl2017/historicalNoVolc_aave_tas.nc'
        else:
            hist_data_loc = cdataprefix+exp_attr[1]+'_volc/BethkeEtAl2017/historicalVolc_aave_tas.nc'
        
    else:
        print("Error: unknown model to this eval script "+ exp_attr)
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

    return (sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist)

def select_data(styr, model_run,sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist):
    this_sim_yr = average_every_n(sims_tas[model_run,:], 12) #converting monthly to yearly
    this_hsim_yr = average_every_n(sims_tas_hist[model_run,:], 12)
    stime_yrs = np.floor(average_every_n(stime_mon, 12)).astype(int)
    shtime_yrs = np.floor(average_every_n(stime_mon_hist, 12)).astype(int)
    
    start_sim = styr - stime_yrs[0] +1 #year that we should switch from observations uncertainty to const simulation uncertainty
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
        long_past_index = ((gen_orig_number(model_run,np.shape(sims_tas)[0])-1 )// 20) #either 1, 2, or 3, still in right order
        long_past_data_loc = cdataprefix+'NorESM_volc/NorESM1-M-historical/hist_aave_tas.nc'
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

    return (simall,start_sim,years, sim_corrected, simh_corrected)
    




    
def run_one_single_ens_member(plotting_figs, experiment_type, start_run, ax1, ax4, colorraw=None):
    from hist_evaluation_script import rank2
    # First evaluation

    exp_attr = experiment_type.split("_")
    
              
    
    if start_run < 0:
        start_run = -start_run
        max_runs = 1+start_run
        plotting_figs= True
    else:
        max_runs = 2+start_run
        plotting_figs= False
    methods_folder= running_subset



    print("starting computation for "+experiment_type)
    print("max_runs" + str(max_runs))
    ens_standard = eval_standard(experiment_type)

    (sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist) = collect_data(exp_attr)
        
    for model_run in range(start_run,max_runs):
        print("Model number:")
        print(model_run)
        print("\n\n\n")
        

        (simall,start_sim,years, sim_corrected, simh_corrected)= select_data(years_past[-1],model_run,sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist)
            
        
        temps_CIl_hist = simall[0:len(years_past)]+  temps_CIl_past - temps_obs_past
        temps_CIu_hist = simall[0:len(years_past)]+  temps_CIu_past - temps_obs_past
        futCIl = sim_corrected[start_sim:] - np.mean(temps_obs_past[-10:]- temps_CIl_past[-10:]) #constant small uncertainty
        futCIu = sim_corrected[start_sim:] - np.mean(temps_obs_past[-10:]- temps_CIu_past[-10:])

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
                   # existing_results.pop("EBMKF_ta4",None) #redo EMBKF_ta3
                   # existing_results.pop("EBMKF_ta",None)
                   # existing_results.pop("FaIR_all_unB",None) #redo EMBKF_ta3
                   # existing_results.pop("FaIR_anthro_unB",None)
                   # existing_results.pop("FaIR_nonat_unB",None)
                   # existing_results.pop("FaIR_comb_unB",None)
                    existing_results.pop("GWI_anthro",None) 
                    existing_results.pop("GWI_anthro_AR6",None)
                    existing_results.pop("GWI_anthro_CGWL",None)
                    existing_results.pop("GWI_anthro_SR15",None)
                    existing_results.pop("GWI_tot",None)     
                    existing_results.pop("GWI_tot_AR6",None)
                    existing_results.pop("GWI_tot_CGWL",None)
                    existing_results.pop("GWI_tot_SR15",None)
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
        

        inum=12 #internal interpolation within years
        standard = results['cent20y']['LT_trend'][2] #retrospective
        smooth_std = np.nanmean(np.abs(np.diff(np.diff(standard))))

        if plotting_figs:
            if ax1==None:
                fig, (ax1,ax4)= plt.subplots(2, 1, figsize=(10,10), gridspec_kw={ "hspace": 0.3})
            ax1.plot(years, standard, 'k-',zorder=3,lw=1.5)
            if colorraw!=None:
                ax1.plot(years,simall,'-', color = colorraw, linewidth=2.5,zorder=3)
                ax1.plot(years,simall,'-', color = 'k', linewidth=.5,zorder=3)

        fineyrs_all = np.arange(years[0],years[-1]+1/inum,1/inum)
        std_intp0 = np.interp(fineyrs_all,years,standard)
        std_intp = std_intp0[~np.isnan(std_intp0)]
        fineyrs_c0 = fineyrs_all[~np.isnan(std_intp0)]
        fineyrs_c = fineyrs_c0[1:]
        thrshs= np.arange(1.1,np.nanmax(standard) ,.1) #thresholds
        print(f"evaluating {len(thrshs)} of 0.1°C thresholds, starting at 1.1°C")
        closest_years = [-1/inum/2+fineyrs_c[np.logical_and(std_intp[0:-1]<i, std_intp[1:]>=i)][0] for i in thrshs]
                   #will have a variable number of steps, at least including 0.5
        closest_yrs_rnd = np.round(closest_years)


        #repeat with ensemble standards
        std_intp0e = np.interp(fineyrs_all,years,ens_standard)
        std_intpE = std_intp0e[~np.isnan(std_intp0e)]
        fineyrs_ce0 = fineyrs_all[~np.isnan(std_intp0e)]
        fineyrs_ce = fineyrs_ce0[1:]
        thrshsE= np.arange(1.1,np.nanmax(ens_standard) ,.1) #thresholds
        closest_yearsE = [-1/inum/2+fineyrs_ce[np.logical_and(std_intpE[0:-1]<i, std_intpE[1:]>=i)][0] for i in thrshsE]
                    #will have a variable number of steps, at least including 0.5
        closest_yrs_rndE = np.round(closest_yearsE)


        
        #breakpoint()
        
        
        central_yr_estimates =[]
        ax4_handles=[]
        ax4_labels=[]
        df_results = pd.DataFrame(columns=['method_name', 'method_class','c/r','smooth_r','avg_unc.(1se)','#q<0.5', '#q<0.1', 'q_min', 'q_small5',
                                           'log-likeli','RMS','bias','tlog-l',
                                           '100log-l','100RMS','75RMS','100bias',
                                           'l15','l20',
                                           'bias50','Edyrs15','Edyrs20','Mdyrs15','Mdyrs20','Fdyrs15','Fdyrs20',
                                           'EdyrsA','RMSyrsA','ncEdyrs',
                                            'e100log-l','e100RMS','e75RMS','e100bias',
                                           'el15','el20',
                                           'ebias50','eEdyrs15','eEdyrs20', 'eMdyrs15','eMdyrs20','eFdyrs15','eFdyrs20',
                                           'eEdyrsA','eRMSyrsA','nceEdyrs'

                                           ])
        i=0
        ci=0
        labelcolors=[]
        sorted_results = sorted(results.items(), key=lambda item: (item[1]['method_class'], rank2(item[0])))

        lhund=-100
        if (exp_attr[1]=='NorESM'):
            lhund=-99
            
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

                    ellikelihood = result['log_likelihood'](years, ens_standard,k)
                    
                    if(sum(~np.isnan(central_est))>0):
                        print(f"{method_name} sampled, d llike: {np.nanmean(llikelihood)-np.nanmean(llikelihood2)}")

                        
                else:
                    central_est = result[k*2]
                    se = result[k*2+1]
                    # want to save and store the more accurate calculation if regen: else: load calculations
                    pvals = stats.norm.sf(abs((standard-central_est)/ se))*2
                    llikelihood = stats.norm.logpdf(standard,loc=central_est,scale=se)
                    ellikelihood = stats.norm.logpdf(ens_standard,loc=central_est,scale=se)

                    
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
                    aedyrs=np.zeros(len(thrshs))
                    maedyrs=np.zeros(len(thrshs))
                    faedyrs=np.zeros(len(thrshs))
                    ncross = 0
                    
                    aedyrsE=np.zeros(len(thrshsE))
                    ncrossE = 0
                    maedyrsE=np.zeros(len(thrshsE))
                    faedyrsE=np.zeros(len(thrshsE))

                    if(labelcurr_or_retro=="c"):
                        if (method_name!="raw01y" and plotting_figs):
                            sline, =ax1.plot(years, central_est, label=method_name, color = gen_color(method_data['method_class'], dark=False),lw=0.5)
                            spaglines.append(sline)

                        if((method_name in sel_methods) and plotting_figs):
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


                   
                        #calculate at intermediary 0.?°C
                        win_sz=25
                        for j in range(len(closest_years)):
                            evalmin=int(closest_years[j])-win_sz
                            evalmax=min(int(closest_years[j])+win_sz,2000-lhund)
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
                            try:
                                psteps_intp = np.interp(fineevalyrs,evalyrs,this_method_p_steps)
                            except:
                                breakpoint()
                            now_cross = np.logical_and(psteps_intp[0:-1]<0.5, psteps_intp[1:]>=0.5)
                            ncross= ncross + np.sum(now_cross)
                            if j==4 and np.sum(now_cross)>1 :
                                ncross = ncross + 0.1*np.sum(now_cross)
                            fineevalyrsh=fineevalyrs[0:-1]/2 + fineevalyrs[1:]/2
                            diffcross = (fineevalyrsh[now_cross] - closest_years[j])
                            evalmean = np.nanmean(diffcross) #if crossing multiple times take the mean
                            if np.isnan(evalmean):
                                evalmean=win_sz #just put 15 yrs difference
                            aedyrs[j] = evalmean
                            if len(diffcross)>0:
                                mevalmean = diffcross[np.nanargmax(abs(diffcross))] #if crossing multiple times take the worst one
                                feval = diffcross[0]
                            else:
                                mevalmean=15 #just put 15 yrs difference
                                feval = 15
                            maedyrs[j] = mevalmean
                            faedyrs[j]=feval
                            

                        for j in range(len(closest_yearsE)):
                            evalmin=int(closest_yearsE[j])-win_sz
                            evalmax=min(int(closest_yearsE[j])+win_sz,2000-lhund)
                            evalyrs = np.arange(evalmin,evalmax)
                            fineevalyrs = np.arange(evalyrs[0],evalyrs[-1]+1/inum,1/inum) 
                            this_method_p_steps = np.full(np.shape(evalyrs),np.nan)
                            if isinstance(result, dict):
                                this_method_p_steps = result['pvalue'](evalyrs, np.full(np.shape(evalyrs),thrshsE[j]),k, two_sided=False)
                            else:
                                this_method_p_steps = stats.norm.cdf(((central_est[evalyrs-1850]-thrshsE[j])/ se[evalyrs-1850]))
                            # Replace NaNs at the start with 0s
                            first_non_nan = np.argmax(~np.isnan(this_method_p_steps))
                            this_method_p_steps[:first_non_nan] = 0
                            # Replace NaNs at the end with 1s if they exist
                            last_non_nan = len(this_method_p_steps) - np.argmax(~np.isnan(this_method_p_steps)[::-1]) - 1
                            this_method_p_steps[last_non_nan + 1:] = 1
                            psteps_intp = np.interp(fineevalyrs,evalyrs,this_method_p_steps)
                            now_cross = np.logical_and(psteps_intp[0:-1]<0.5, psteps_intp[1:]>=0.5)
                            ncrossE= ncrossE + np.sum(now_cross)
                            if j==4 and np.sum(now_cross)>1 :
                                ncrossE = ncrossE + 0.1*np.sum(now_cross)
                            fineevalyrsh=fineevalyrs[0:-1]/2 + fineevalyrs[1:]/2
                            diffcrossE = (fineevalyrsh[now_cross] - closest_yearsE[j])
                            evalmean = np.nanmean(diffcrossE) #if crossing multiple times take the mean
                            if np.isnan(evalmean):
                                evalmean=win_sz #just put 15 yrs difference
                            aedyrsE[j] = evalmean
                            if len(diffcrossE)>0:
                                mevalmean = diffcrossE[np.nanargmax(abs(diffcrossE))] #if crossing multiple times take the worst one
                                feval = diffcrossE[0]
                            else:
                                mevalmean=15 #just put 15 yrs difference
                                feval = 15
                            maedyrsE[j] = mevalmean
                            faedyrsE[j]=feval
                            
                      #  elif("lowess" in method_name):
                      #      plt.fill_between(years, central_est-se, central_est+se, alpha=0.6)
                    else: #blank values
                        aedyrs=np.zeros(len(thrshs))
                    short_method_class = method_data['method_class'][0:2] +"/"+method_data['method_class'].split('/')[-1]

                    
                    detailLL = np.zeros(5)
                    detaileLL = np.zeros(5)
                    for dll in range(4):
                        idll = (dll+1)*5-1 #particular thresholds
                        if idll<len(closest_yrs_rnd):
                            detailLL[dll] =  np.exp(llikelihood[int(closest_yrs_rnd[idll])-1850])
                        if idll<len(closest_yrs_rndE):
                            detaileLL[dll] =  np.exp(ellikelihood[int(closest_yrs_rndE[idll])-1850])
                    
                    candidate_row= [ method_name,short_method_class,labelcurr_or_retro,smooth_est/smooth_std,avg_uncert,
                                     qvals_count_yrs05,qvals_count_yrs01,  qvals_smallest,qvals_smallest5, np.nanmean(llikelihood), np.sqrt(np.nanmean((central_est-standard)**2)),
                                       np.nanmean(central_est-standard) , np.nansum(llikelihood),
                                         np.nansum(llikelihood[lhund:-1]),np.sqrt(np.nanmean((central_est[lhund:-1]-standard[lhund:-1])**2)),
                                         np.sqrt(np.nanmean((central_est[(lhund+25):-1]-standard[(lhund+25):-1])**2)),
                                         np.nanmean(central_est[lhund:-1]-standard[lhund:-1]),
                                         detailLL[0], detailLL[1],
                                         np.nanmean(central_est[-50:]-standard[-50:]),aedyrs[4], (aedyrs[9] if (len(aedyrs)>9) else -1),
                                         maedyrs[4], (maedyrs[9] if (len(maedyrs)>9) else -1),faedyrs[4], (faedyrs[9] if (len(faedyrs)>9) else -1),
                                         np.mean(aedyrs),np.sqrt(np.mean(aedyrs**2)),ncross,
                                         np.nansum(ellikelihood[lhund:-1]),np.sqrt(np.nanmean((central_est[lhund:-1]-ens_standard[lhund:-1])**2)),
                                         np.sqrt(np.nanmean((central_est[(lhund+25):-1]-ens_standard[(lhund+25):-1])**2)),
                                         np.nanmean(central_est[lhund:-1]-ens_standard[lhund:-1]),
                                         detaileLL[0], detaileLL[1],
                                         np.nanmean(central_est[-50:]-ens_standard[-50:]),aedyrsE[4], (aedyrsE[9] if (len(aedyrsE)>9) else -1),
                                         maedyrsE[4], (maedyrsE[9] if (len(maedyrsE)>9) else -1), faedyrsE[4], (faedyrsE[9] if (len(faedyrsE)>9) else -1),
                                         np.mean(aedyrsE),np.sqrt(np.mean(aedyrsE)**2),ncrossE]
                    #breakpoint()
                    df_results.loc[i]=candidate_row
                    i=i+1
       
        if(crossing_figs):
            #extra stuff to nmke the figures look  nice
            ax05.set_yticks(range(len(cmethods)))
            cmethods2=np.array(cmethods)
            cmethodsrs=cmethods2 #[index_mapping]
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
                
        df_res_show = df_results.drop(columns=['q_min','q_small5'])
        df_res_show['smooth_r'] = df_results['smooth_r'].round(3)
        df_res_cur = df_res_show[df_results['c/r']=='c']
        
        print(df_res_cur.sort_values('100log-l',ascending=False))

        if(plotting_figs):
            ax1.grid(color='silver',zorder=-1)
            ax4.grid(color='silver',zorder=-1)
            #df_results.to_csv('method_statistics_results.csv', index=False)
            ax1.legend(fontsize=6.5,ncol=4)
            ax4.set_ylim(bottom=-0.11,top=0.11)
            xmin, xmax = [2000,2100] #ax1.get_xlim() #2000
            ax4.set_xlim([xmin, xmax])
            ax4.set_xticks(np.arange(2000,2105,25))
            ax1.set_xticks(np.arange(2000,2105,25))
            ax1.set_xlim([xmin, xmax])
            ax4.hlines(y=0,xmin=years[np.argmax(~np.isnan(standard))],xmax=years[-np.argmax(~np.isnan(np.flip(standard)))],  color='k', linestyle='-', lw=3)
            ax4.legend(ax4_handles,ax4_labels)
            ax1.set_xlabel("Year")
            ax1.set_title("Evaluated Methods to find Current "+exp_attr[1]+" "+add_dash_dot(exp_attr[2]), pad=10)
            ax4.set_title("`Error` of Top-Performing Methods\n versus 20-yr running mean", pad=10)
            ax4.set_xlabel("Year")
            ax1.set_ylabel("Temperature (°C) Anomaly\n relative to 1850-1900")
            ax4.set_ylabel("Temperature (°C) Difference\n relative to 20-yr running mean")

            # Adding interactive tooltips
            cursor = mplcursors.cursor(spaglines, hover=True)

            # Customize the displayed tooltip
            @cursor.connect("add")
            def on_add(sel):
                # Show the label of the hovered line
                sel.annotation.set_text(sel.artist.get_label())
        

                
        


        print(time.process_time() - start)
        
        df_res_show2 = df_results.copy()
        df_res_show2['smooth_r'] = df_results['smooth_r'].round(3)
        df_res_cur2 = df_res_show2[df_results['c/r']=='c']
        dfres2 = df_res_cur2.sort_values('100log-l',ascending=False) #do not need to sort saved output results - we will have to read and reference anyway
        dfres2.to_csv('Results/current_fut_statistics_'+experiment_type+str(model_run)+'.csv', index=False,mode='w+' )
        df_results.to_csv('Results/all_fut_statistics_'+experiment_type+str(model_run)+'.csv', index=False,mode='w+' )
        print("finished this run") 
    return closest_years, closest_yearsE #pass stuff back so we can draw more, thresholds at 1.1 1.2 ....
        
def colorgen(c,step,const):
    y=(((c+1)%3-1)*step[0], ((c//3+1)%3-1)*step[1], ((c//9+1)%3-1)*step[2])
    color=tuple( max(min(sum(t), 255),0) / 255 for t in zip(y,const))
    return color #

def colorcomb(col1,const):
    color=tuple( max(min(sum(t)/2*256, 255),0) / 255 for t in zip(col1,const))
    return color #


if __name__ == '__main__':
    plotting_figs=False
    experiment_type = sys.argv[1] #'fut_ESM1-2-LR_SSP126_constVolc' #fut_NorESM_RCP45_Volc
    start_run = int(sys.argv[2])
    exp_attr = experiment_type.split("_")
    
    if exp_attr[0]=="fut":
        run_one_single_ens_member(plotting_figs, experiment_type, start_run, None, None)
        plt.show()
        
    elif exp_attr[0]=="futplotcomb" and exp_attr[1]=="ESM1-2-LR":
        fig1, (axens, ax1a,ax1b)= plt.subplots(3, 1, figsize=(10,10), gridspec_kw={ "height_ratios" :[.6,1,1],"hspace": 0.35})
        fig4, (ax4a,ax4b)= plt.subplots(2, 1, figsize=(10,10), gridspec_kw={  "hspace": 0.3})
        exp1 = "fut_ESM1-2-LR_SSP126_constVolc"
        exp2 = "fut_ESM1-2-LR_SSP370_constVolc"
        col370 = (194./255, 52./255, 109./255)
        col126 = (52./255,120./255, 130./255)
        (sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist) = collect_data(exp2.split("_"))
        for model_run in range(0,50):
            (simall,start_sim,years, sim_corrected, simh_corrected)= select_data(2024,model_run,sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist)
            if model_run == -(start_run):
                axens.plot(years,simall,'-', color = col370, linewidth=2.5,zorder=3)
                axens.plot(years,simall,'-', color = 'k', linewidth=.5,zorder=3)
            else:
                axens.plot(years,simall,'-', color=colorgen(model_run,(30,70,30),(194, 52, 109)), linewidth=0.5,zorder=1) #(52./255, 235./255, 235./255)
            
        (sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist) = collect_data(exp1.split("_"))
        for model_run in range(0,50):
            (simall,start_sim,years, sim_corrected, simh_corrected)= select_data(2024,model_run,sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist)
            if model_run == -(start_run):
                axens.plot(years,simall,'-', color = col126, linewidth=2.5,zorder=3)
                axens.plot(years,simall,'-', color = 'k', linewidth=.5,zorder=3)
            else:
                axens.plot(years,simall,'-', color=colorgen(model_run,(70,30,30),(52,120,130)), linewidth=0.5,zorder=1) #(52./255, 235./255, 235./255)
            
        axens.set_yticks(np.arange(0, 4, step=.5), minor=True)
        axens.set_xticks(np.arange(2000,2105,25))
        axens.grid(color='silver',zorder=-1,which='both')
       
        nud = 0.5
        closest_years126, closest_yearsE126 = run_one_single_ens_member(plotting_figs, exp1, start_run, ax1a, ax4a, col126)
        closest_years370, closest_yearsE370 =run_one_single_ens_member(plotting_figs, exp2, start_run, ax1b, ax4b, col370)
        axens.plot([closest_yearsE126[4],closest_yearsE126[4]],[0.9,1.5],color = col126)
        axens.plot([closest_yearsE126[4]+nud,closest_yearsE126[4]+nud],[0.9,1.5],color = col126, lw=0.5)
        axens.plot([closest_yearsE126[4]-nud,closest_yearsE126[4]-nud],[0.9,1.5],color = col126, lw=0.5)
        axens.text(closest_yearsE126[4], 0.6, str(round(closest_yearsE126[4], ndigits=1)), color = col126,horizontalalignment='center')
        axens.plot([closest_yearsE370[4],closest_yearsE370[4]],[0.9,1.5],color = col370)
        axens.plot([closest_yearsE370[4]+nud,closest_yearsE370[4]+nud],[0.9,1.5],color = col370, lw=0.5)
        axens.plot([closest_yearsE370[4]-nud,closest_yearsE370[4]-nud],[0.9,1.5],color = col370, lw=0.5)
        axens.text(closest_yearsE370[4], 0.6, str(round(closest_yearsE370[4], ndigits=1)), color = col370,horizontalalignment='center')

        axens.plot([closest_years126[4],closest_years126[4]],[1.5,3],color = col126,lw=2.5)
        axens.plot([closest_years126[4],closest_years126[4]],[2.5,3],color = "black", lw=0.5)
        axens.text(closest_years126[4], 3.1, str(round(closest_years126[4], ndigits=1)), color = col126,horizontalalignment='center')
        axens.plot([closest_years370[4],closest_years370[4]],[1.5,3],color = col370,lw=2.5)
        axens.plot([closest_years370[4],closest_years370[4]],[2.5,3],color = "black", lw=0.5)
        axens.text(closest_years370[4], 3.1, str(round(closest_years370[4], ndigits=1)), color = col370,horizontalalignment='center')
        
        axens.set_ylabel(ax1a.get_ylabel()) #matching ylabels
        ax1a.set_ylim([0.5,3.5])
        ax1a.set_xlabel("")
        ax1b.set_ylim([0.5,3.5])
        axens.set_ylim([0.5,3.5])
        axens.set_xlim(ax1a.get_xlim())
        ax1b.get_legend().remove()
        plt.show()

    elif exp_attr[0]=="futplotcomb" and exp_attr[1]=="NorESM":
        fig1, (axens, ax1a,ax1b)= plt.subplots(3, 1, figsize=(10,10), gridspec_kw={ "height_ratios" :[.6,1,1],"hspace": 0.35})
        fig4, (ax4a,ax4b)= plt.subplots(2, 1, figsize=(10,10), gridspec_kw={  "hspace": 0.3})
        exp1 = "fut_NorESM_RCP45_VolcConst"
        exp2 = "fut_NorESM_RCP45_Volc"
        
        file_path="Results/ensemble_mean_"+exp2+".csv"
                #if os.path.exists(file_path): #will raise error if this file isn't here
        df = pd.read_csv(file_path, header=None)
        standards_np = df.to_numpy()
        mask = standards_np > 1.5
        result = np.argmax(mask, axis=1)
        min_run = np.argmin(result)#11
        max_run = np.argmax(result)
        
        colnovolc = (245./255, 212./255, 66./255)
        colvolc = (16./255, 133./255, 28./255)
        (sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist) = collect_data(exp2.split("_"))
        for model_run in range(0,50):
            (simall,start_sim,years, sim_corrected, simh_corrected)= select_data(2024,model_run,sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist)
            if model_run == min_run:
                min_line,=axens.plot(years,simall,'-', color = colorcomb(colvolc,colnovolc),zorder=3, linewidth=2.5) #,label = "Volc. Sample Min")
                black_line, = axens.plot(years,simall,'-', color = 'k', linewidth=.5,zorder=3)
            elif model_run == max_run :# -(start_run):
                max_line,=axens.plot(years,simall,'-', color = colvolc,zorder=3,linewidth=2.5) #, label = "Volc. Sample Max")
                axens.plot(years,simall,'-', color = 'k', linewidth=.5,zorder=3)
            else:
                axens.plot(years,simall,'-', color=colorgen(model_run,(30,70,30),(16,133,28)), linewidth=0.5,zorder=1, alpha=0.75) #(52./255, 235./255, 235./255)



        (sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist) = collect_data(exp1.split("_"))
        for model_run in range(0,50):
            (simall,start_sim,years, sim_corrected, simh_corrected)= select_data(2024,model_run,sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist)
            if model_run == 0:
                const_novolc_line,= axens.plot(years,simall,'-', color = colnovolc, linewidth=0.5,zorder=3,alpha=0.5) #,label = "Volc. Const.")
                #black_line = axens.plot(years,simall,'-', color = 'k', linewidth=.5,zorder=3)
            else:
                axens.plot(years,simall,'-', color=colorgen(model_run,(70,30,30),(245, 212, 66)), linewidth=0.5,zorder=1,alpha=0.5) #(52./255, 235./255, 235./255)
            
        axens.set_yticks(np.arange(0, 4, step=.5), minor=True)
        axens.set_xticks(np.arange(2000,2105,25))
        axens.set_xlim([2000,2100])
        axens.grid(color='silver',zorder=-1,which='both')
       
        nud = 0.5
        closest_years126, closest_yearsE126 = run_one_single_ens_member(plotting_figs, exp2, -min_run, ax1a, ax4a, colorcomb(colvolc,colnovolc))
        closest_years370, closest_yearsE370 = run_one_single_ens_member(plotting_figs, exp2, -max_run, ax1b, ax4b, colvolc)
##        axens.plot([closest_yearsE126[4],closest_yearsE126[4]],[0.9,1.5],color = colorcomb(colvolc,colnovolc))
##        axens.plot([closest_yearsE126[4]+nud,closest_yearsE126[4]+nud],[0.9,1.5],color = colorcomb(colvolc,colnovolc), lw=0.5)
##        axens.plot([closest_yearsE126[4]-nud,closest_yearsE126[4]-nud],[0.9,1.5],color = colorcomb(colvolc,colnovolc), lw=0.5)
##        axens.text(closest_yearsE126[4], 0.6, str(round(closest_yearsE126[4], ndigits=1)), color = colorcomb(colvolc,colnovolc),horizontalalignment='center')
        axens.plot([closest_yearsE370[4],closest_yearsE370[4]],[0.3,1.5],color =colvolc)
        axens.plot([closest_yearsE370[4]+nud,closest_yearsE370[4]+nud],[0.3,1.5],color = colvolc, lw=0.5)
        axens.plot([closest_yearsE370[4]-nud,closest_yearsE370[4]-nud],[0.3,1.5],color = colvolc, lw=0.5)
        axens.text(closest_yearsE370[4], 0.05, str(round(closest_yearsE370[4], ndigits=1)), color = colvolc,horizontalalignment='center')

        axens.plot([closest_years126[4],closest_years126[4]],[1.5,2.4],color = colorcomb(colvolc,colnovolc),lw=2.5)
        axens.plot([closest_years126[4],closest_years126[4]],[2.,2.4],color = "black", lw=0.5)
        axens.text(closest_years126[4], 2.45, str(round(closest_years126[4], ndigits=1)), color = colorcomb(colvolc,colnovolc),horizontalalignment='center')
        axens.plot([closest_years370[4],closest_years370[4]],[1.5,2.4],color = colvolc,lw=2.5)
        axens.plot([closest_years370[4],closest_years370[4]],[2.,2.4],color = "black", lw=0.5)
        axens.text(closest_years370[4], 2.45, str(round(closest_years370[4], ndigits=1)), color = colvolc,horizontalalignment='center')
        
        axens.set_ylabel(ax1a.get_ylabel()) #matching ylabels
        ax1a.set_ylim([0,2.75])
        ax1a.set_title("Evaluated Methods to find Current NorESM RCP 4.5 - Min. Volcanic Activity")
        ax1b.set_title("Evaluated Methods to find Current NorESM RCP 4.5 - Max. Volcanic Activity")
        axens.set_title("Data from Bethke et. al. 2017: \nPotential volcanic impacts on future climate variability")
        ax1a.set_xlabel("")
        ax1b.set_ylim([0,2.75])
        axens.set_ylim([0,2.75])
        axens.legend([(min_line,black_line),(max_line,black_line), (const_novolc_line,)], ["Volc. Sample Min","Volc. Sample Max", "Volc. Const."])
        
        ax1b.get_legend().remove()
        plt.show()    

