from hist_evaluation_script import *
regen=True
annotate_fig=False
crossing_figs=False
sel_methods = ["GWI_tot_CGWL","EBMKF_ta4"]
from netCDF4 import Dataset
import sys
from scipy.optimize import minimize
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

def select_data(styr, model_run,sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist,exp_attr):
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
    


from numpy.polynomial.hermite import hermgauss


    
def run_one_single_ens_member(plotting_figs, experiment_type, start_run, ax1, ax4, colorraw=None,exp_index=0):
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
        

        (simall,start_sim,years, sim_corrected, simh_corrected)= select_data(years_past[-1],model_run,sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist,exp_attr)
            
        
        temps_CIl_hist = simall[0:len(years_past)]+  temps_CIl_past - temps_obs_past
        temps_CIu_hist = simall[0:len(years_past)]+  temps_CIu_past - temps_obs_past
        futCIl = sim_corrected[start_sim:] - np.mean(temps_obs_past[-10:]- temps_CIl_past[-10:]) #constant small uncertainty
        futCIu = sim_corrected[start_sim:] - np.mean(temps_obs_past[-10:]- temps_CIu_past[-10:])

        results_path = f'Results/results{experiment_type}{model_run}.pickle'
    
        if regen==1:
            results = run_methods(years, simall,
                                  (np.concatenate((temps_CIl_hist,futCIl)), np.concatenate((temps_CIu_hist,futCIu)) ),
                                  model_run, experiment_type, methods_folder  =
                                  ['Methods/43_Forcing_Based/2_Kalman/EBMKF_ta4_method.py','Methods/43_Forcing_Based/3_Human_Induced/GWI_tot_CGWL_method.py',
                                   'Methods/42_Temp_Alone/1_Run_Means/cent20y_method.py'],give_methods_path=True)


        inum=12 #internal interpolation within years
        standard = results['cent20y']['LT_trend'][2] #retrospective
        standard_se = results['cent20y']['LT_trend'][3] #retrospective
        smooth_std = np.nanmean(np.abs(np.diff(np.diff(standard))))

        import neworder
        sorted_results = neworder.sort_results(results) #sorted(results.items(), key=lambda item: (item[1]['method_class'], rank2(item[0])))

        fineyrs_all = np.arange(years[0],years[-1]+1/inum,1/inum)
        std_intp0 = np.interp(fineyrs_all,years,standard)
        std_intp = std_intp0[~np.isnan(std_intp0)]
        fineyrs_c0 = fineyrs_all[~np.isnan(std_intp0)]
        fineyrs_c = fineyrs_c0[1:]
        crit_j =0
        thrshs= np.arange(1.5,np.nanmax(standard) ,.1) #thresholds
        print(f"evaluating {len(thrshs)} of 0.1°C thresholds, starting at {thrshs[0]}°C")
        closest_years = [-1/inum/2+fineyrs_c[np.logical_and(std_intp[0:-1]<i, std_intp[1:]>=i)][0] for i in thrshs]
                   #will have a variable number of steps, at least including 0.5
        closest_yrs_rnd = np.round(closest_years)

        
        
        central_yr_estimates =[]
        ax4_handles=[]
        ax4_labels=[]

        i=0
        ci=0
        labelcolors=[]
        #print(len(years))
        lhund=-101
        l75=-76
        lten=-10
        if (exp_attr[1]=='NorESM'):
            lhund=-100
            l75=-75
            lten = -10
            #so there are fewer points to average over in the NorESM case
            

        aedyrs=np.zeros((2,len(thrshs)))
        maedyrs=np.zeros((2,len(thrshs)))
        faedyrs=np.zeros((2,len(thrshs)))
        ncrosses = np.zeros(2)
        kls = np.zeros(2)
        kls75 = np.zeros(2)
        rmses = np.zeros(2)
        rmses75 = np.zeros(2)
                    
        for method_name, method_data in sorted_results:
            if method_name == 'cent20y':
                continue
            print(f"Results from {method_name} (Method Class: {method_data['method_class']}):")
            result = method_data['LT_trend']
            labelcurr_or_retro="" #only set if the current method is not blank
            isgauss = True 
            for k in [0]:
                central_est=np.full(len(years),np.nan) #make this blank to start
                if isinstance(result, dict):
                    # Empirical method: Call the functions at each time point
                    central_est = result['mean'](years,k)
                    se = result['se'](years,k)
                    pvals = result['pvalue'](years, standard,k)
                    llikelihood = result['log_likelihood'](years, standard,k)
                    llikelihood2 = stats.norm.logpdf(standard,loc=central_est,scale=se)

                    ellikelihood = result['log_likelihood'](years, ens_standard,k)
                    isgauss = False
                    if(sum(~np.isnan(central_est))>0):
                        print(f"{method_name} sampled, d llike: {np.nanmean(llikelihood)-np.nanmean(llikelihood2)}")

                    def rescale_log_likelihood(scale):
                        deviance = standard[lhund:lten] - central_est[lhund:lten]
                        resc_standard = central_est[lhund:lten] + deviance/scale #larger scale makes deviance appear smaller
                        log_lik=result['log_likelihood'](years[lhund:lten], resc_standard,k) - np.log(scale) #pdf expands outward with small scale, so must compress vertically
                        return -np.nansum(log_lik) #minimize this quanity


                    x, w = hermgauss(100)
                    z = np.sqrt(2.0) * x
                    alpha = w / np.sqrt(np.pi)
                    def rescale_KL(scale, stand0, stand_se0, lasts, laste,eps=1e-300,eeoffset=0):
                        """
                        Sum_y KL( P_y || Q_y ),  P_y = N(standard[y], standard_se[y]^2),
                        Q_y = empirical via scipy.stats.gaussian_kde on sharp[y, :], using Gauss–Hermite quadrature.
                        """
                        stand = stand0
                        stand_se= stand_se0
                        
                        total = 0.0
                        for y in range(len(stand)+lasts,len(stand)+laste-eeoffset):
                            if np.isnan(stand[y]):
                                continue
                            mu, sd = float(stand[y]), float(max(stand_se[y], 1e-12))
                            y_nodes = mu + sd * z                         # transform nodes to P_y support
                            log_p = stats.norm.logpdf(y_nodes, loc=mu, scale=sd)
                            deviances = y_nodes - central_est[y]
                            ynodes2 = central_est[y] + deviances/scale
                            log_q = result['log_likelihood']([years[y]]*100, ynodes2,0) - np.log(scale)
                            total += float(np.sum(alpha * (log_p - log_q)))
                        return total


                        
                else:
                    central_est = result[k*2]
                    se = result[k*2+1]
                    # want to save and store the more accurate calculation if regen: else: load calculations
                    pvals = stats.norm.sf(abs((standard-central_est)/ se))*2
                    llikelihood = stats.norm.logpdf(standard,loc=central_est,scale=se)
                    ellikelihood = stats.norm.logpdf(ens_standard,loc=central_est,scale=se)
                    
                    def rescale_log_likelihood(scale_alt):
                        log_lik=stats.norm.logpdf(standard[lhund:lten],loc=central_est[lhund:lten],scale=se[lhund:lten]*scale_alt)
                        return -np.nansum(log_lik)
                    
                    def rescale_KL(scale,standardi0,standard_sei0,lasts, laste, eeoffset=0):
                        ivcenters=central_est[lasts:laste]
                        ivse =se[lasts:laste]
                        standardi=standardi0[lasts:laste]
                        standard_sei=standard_sei0[lasts:laste]
                        yearKL = np.log(scale*ivse/standard_sei) + (standard_sei**2 + (standardi-ivcenters)**2 )/2/(scale*ivse)**2 - 0.5
                        return np.nansum(yearKL,axis=None)
                    
                    
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

                    best_alter_scale = minimize(rescale_log_likelihood, x0=1, bounds=[(0.01, 100.0)]).x[0]
                    kls[i] = rescale_KL(best_alter_scale,standard,standard_se,lhund,lten )
                    kls75[i] = rescale_KL(best_alter_scale,standard,standard_se,l75,lten )
                    m=i
                    if(labelcurr_or_retro=="c"):
  
                        #calculate at intermediary 0.?°C
                        win_sz=25
                        ncross=0
                        for j in range(len(closest_years)):
                            evalmin=int(closest_years[j])-win_sz
                            evalmax=min(int(closest_years[j])+win_sz,1850+len(standard))
                            evalyrs = np.arange(evalmin,evalmax)
                            fineevalyrs = np.arange(evalyrs[0],evalyrs[-1]+1/inum,1/inum) 
                            this_method_p_steps = np.full(np.shape(evalyrs),np.nan)
                            if isinstance(result, dict):
                                evalpts = np.full(np.shape(evalyrs),thrshs[j])
                                deviances = evalpts - central_est[evalyrs-1850] #relative to this method's central estimate
                                resc_standard_samples = central_est[evalyrs-1850] + deviances/best_alter_scale
                                this_method_p_steps = result['pvalue'](evalyrs, resc_standard_samples,k, two_sided=False)
                            else:
                                this_method_p_steps = stats.norm.cdf((central_est[evalyrs-1850]-thrshs[j])/ (se[evalyrs-1850] *best_alter_scale))
                            # Replace NaNs at the start with 0s

                            first_non_nan = np.argmax(~np.isnan(this_method_p_steps))
                            this_method_p_steps[:first_non_nan] = 0
                            # Replace NaNs at the end with 1s if they exist
                            last_non_nan = len(this_method_p_steps) - np.argmax(~np.isnan(this_method_p_steps)[::-1]) - 1
                            this_method_p_steps[last_non_nan + 1:] = 1
                            psteps_intp = np.interp(fineevalyrs,evalyrs,this_method_p_steps)


                            now_cross = np.logical_and(psteps_intp[0:-1]<0.5, psteps_intp[1:]>=0.5)
                            #ncross= ncross + np.sum(now_cross)
                            if j==crit_j and np.sum(now_cross)>0 : #special tag at 0.5 or 1.5
                                ncross = ncross + np.sum(now_cross)
                            ncrosses[m]=ncross
                            fineevalyrsh=fineevalyrs[0:-1]/2 + fineevalyrs[1:]/2
                            diffcross = (fineevalyrsh[now_cross] - closest_years[j])
                            evalmean = np.nanmean(diffcross) #if crossing multiple times take the mean
                            if np.isnan(evalmean):
                                evalmean=win_sz #just put 15 yrs difference
                            aedyrs[m,j] = evalmean
                            if len(diffcross)>0:
                                mevalmean = diffcross[np.nanargmax(abs(diffcross))] #if crossing multiple times take the worst one
                                feval = diffcross[0]
                            else:
                                mevalmean=15 #just put 15 yrs difference
                                feval = 15
                            maedyrs[m,j] = mevalmean
                            faedyrs[m,j]=feval
                            


                    else: #blank values
                        aedyrs=np.zeros(len(thrshs))
                        
                    short_method_class = method_data['method_class'][0:2] +"/"+method_data['method_class'].split('/')[-1]

                    
                    
##                    candidate_row= [ method_name,short_method_class,labelcurr_or_retro,smooth_est/smooth_std,avg_uncert,
##                                     qvals_count_yrs05,qvals_count_yrs01,  qvals_smallest,qvals_smallest5, np.nanmean(llikelihood), np.sqrt(np.nanmean((central_est-standard)**2)),
##                                       np.nanmean(central_est-standard) , np.nansum(llikelihood),
##                                         np.nansum(llikelihood[lhund:-1]),np.sqrt(np.nanmean((central_est[lhund:-1]-standard[lhund:-1])**2)),
##                                         np.sqrt(np.nanmean((central_est[(lhund+25):-1]-standard[(lhund+25):-1])**2)),
##                                         np.nanmean(central_est[lhund:-1]-standard[lhund:-1]),
##                                         detailLL[0], detailLL[1],
##                                         np.nanmean(central_est[-50:]-standard[-50:]),aedyrs[4], (aedyrs[9] if (len(aedyrs)>9) else -1),
##                                         maedyrs[4], (maedyrs[9] if (len(maedyrs)>9) else -1),faedyrs[4], (faedyrs[9] if (len(faedyrs)>9) else -1),
##                                         np.mean(aedyrs),np.sqrt(np.mean(aedyrs**2)),ncross,
##                                         np.nansum(ellikelihood[lhund:-1]),np.sqrt(np.nanmean((central_est[lhund:-1]-ens_standard[lhund:-1])**2)),
##                                         np.sqrt(np.nanmean((central_est[(lhund+25):-1]-ens_standard[(lhund+25):-1])**2)),
##                                         np.nanmean(central_est[lhund:-1]-ens_standard[lhund:-1]),
##                                         detaileLL[0], detaileLL[1],
##                                         np.nanmean(central_est[-50:]-ens_standard[-50:]),aedyrsE[4], (aedyrsE[9] if (len(aedyrsE)>9) else -1),
##                                         maedyrsE[4], (maedyrsE[9] if (len(maedyrsE)>9) else -1), faedyrsE[4], (faedyrsE[9] if (len(faedyrsE)>9) else -1),
##                                         np.mean(aedyrsE),np.sqrt(np.mean(aedyrsE)**2),ncrossE]
                    #breakpoint()

                    #rescale variance

                    rmses[m] = np.sqrt(np.nanmean((central_est[lhund:lten]-standard[lhund:lten])**2))
                    rmses75[m] = np.sqrt(np.nanmean((central_est[l75:lten]-standard[l75:lten])**2))
                    

                    #output parameters
                    #np.set_printoptions(precision=9)
                    #print(central_est[lhund])
                    #print(standard[lhund])
                    #print(central_est[lten-1])
                    #print(standard[lten-1:])
                    #print(len(standard))
                    #print(method_name)
                    #print(rmses75[m])
                    #breakpoint()
                    #print(i)
                    i=i+1


        np.savez_compressed(f"Results3/nmes2run_exp{exp_index}_r{model_run}.npz",
                        firstcross15_diff=faedyrs[:,crit_j],
                        ncrosses=ncrosses,
                        rmse=rmses,
                        rmse75 = rmses75,
                        kl=kls,
                        kl75 = kls75)

        #np.set_printoptions(precision=9)
        #print(rmses75) 
        #print(rmses) 


        print(time.process_time() - start)
        

        print("finished this run") 

        


if __name__ == '__main__':
    plotting_figs=False
    exp_index = int(sys.argv[1]) #'fut_ESM1-2-LR_SSP126_constVolc' #fut_NorESM_RCP45_Volc
    start_run = int(sys.argv[2])
    scenarios=["fut_ESM1-2-LR_SSP126_constVolc","fut_ESM1-2-LR_SSP245_constVolc","fut_ESM1-2-LR_SSP370_constVolc",
                   "fut_NorESM_RCP45_Volc","fut_NorESM_RCP45_VolcConst"]
    experiment_type=scenarios[exp_index]
    
    run_one_single_ens_member(plotting_figs, experiment_type, start_run, None, None,exp_index=exp_index)
        #plt.show()
        
    

