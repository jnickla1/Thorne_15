from hist_evaluation_script import *
regen=True
annotate_fig=False
crossing_figs=False
from netCDF4 import Dataset
import sys

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])


regen = 1 #0 no regen #1 regen completely #2 overwrite regen to allow for computed methods to not need to be redone!

def eval_standard(experiment_type):
    file_path="Results/ensemble_mean_"+experiment_type+".csv"
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, header=None)
        standards_np = df.to_numpy()
        return np.nanmean(standards_np,axis=0)
        
    
    # First evaluation
    data = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50]) #preindustrial baseline 1850-1899
    temps_obs_past = temps_obs - preind_base #remove this baseline
    temps_CIu_past =data.loc[:,"Upper"].to_numpy() #Lower confidence limit (2.5%)	Upper confidence limit (97.5%)
    temps_CIl_past =data.loc[:,"Lower"].to_numpy()
    years_past=data.loc[:,"Time"].to_numpy()
    
    
    #experiment_type = sys.argv[1] #'fut_ESM1-2-LR_SSP126_constVolc' #fut_NorESM_RCP45_Volc
    exp_attr = experiment_type.split("_")
    if (exp_attr[1]=='ESM1-2-LR'):
        max_runs = 50
        fut_data_loc = config.CLIMATE_DATA_PATH+'/'+exp_attr[1]+'/combined/'+exp_attr[2].lower()+'_aave_tas.nc'
        hist_data_loc = config.CLIMATE_DATA_PATH+'/'+exp_attr[1]+'/combined/historical_aave_tas.nc'
        
    elif (exp_attr[1]=='NorESM'):
        max_runs = 60
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

    standards=[]
    print("starting computation for "+experiment_type)
        
    for model_run in range(0,max_runs):
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
        #results_path = f'Results/results{experiment_type}{model_run}.pickle'
    
        results = run_methods(years, simall,
                                  (np.concatenate((temps_CIl_hist,futCIl)), np.concatenate((temps_CIu_hist,futCIu)) ),
                                  model_run, experiment_type, ["Methods/42_Temp_Alone/1_Run_Means/cent20y_method.py"], completed_methods = set(), give_methods_path = True)

        

        standards.append(results['cent20y']['LT_trend'][2]) #retrospective

    
    standards_np = np.array(standards)
    mean_cross = np.nanmean(standards_np,axis=0)
    df = pd.DataFrame(standards_np)
    df.to_csv(file_path, header=False, index=False)
    
    return mean_cross
        





