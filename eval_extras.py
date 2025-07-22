### Combining real observations with future simulations


# First evaluation
    data = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50]) #preindustrial baseline 1850-1899
    temps_obs_past = temps_obs - preind_base #remove this baseline
    temps_CIu_past =data.loc[:,"Upper"].to_numpy() #Lower confidence limit (2.5%)	Upper confidence limit (97.5%)
    temps_CIl_past =data.loc[:,"Lower"].to_numpy()
    years_past=data.loc[:,"Time"].to_numpy()
    
    model_run = 0
    experiment_type = 'fut-ESM1-2-LR-SSP126-constVolc'
    methods_folder=('Methods/42_Temp_Alone/1_Run_Means','Methods/42_Temp_Alone/2_LT_Fits','Methods/42_Temp_Alone/3_ST_Fits') #,'Methods/43_Forcing_Based','Methods/44_EarthModel_CGWL')
    #big problems with future ENSO index (can't compute MEI) and need to get Bjorn Samset data 
    dataset = Dataset(config.CLIMATE_DATA_PATH+"/ESM1-2-LR/combined/ssp126_aave_tas.nc', 'r')
    variable = dataset.variables['tas']
    sims_tas = variable[:].__array__()
    stimes = dataset.variables['time']
    stime_mon = stimes[:].__array__()/365+1850

    this_sim_yr = average_every_n(sims_tas[0,:], 12)
    stime_yrs = np.floor(average_every_n(stime_mon, 12)).astype(int)
    start_sim = years_past[-1] - stime_yrs[0] +1
    #futCIl = np.full((len(stime_yrs) - start_sim),temps_CIl_past[-1])
    offset_sim = np.mean( this_sim_yr[0:start_sim] - temps_obs_past[-start_sim:])
    sim_corrected = this_sim_yr[start_sim:]-offset_sim
    futCIl = sim_corrected - np.mean(temps_obs_past[-10:]- temps_CIl_past[-10:])
    futCIu = sim_corrected - np.mean(temps_obs_past[-10:]- temps_CIu_past[-10:])
    # Run all the methods

    if regen:
        results = run_methods(np.concatenate((years_past,stime_yrs[start_sim:])), np.concatenate((temps_obs_past,sim_corrected)),
                              (np.concatenate((temps_CIl_past,futCIl)), np.concatenate((temps_CIu_past,futCIu)) ),
                              model_run, experiment_type, methods_folder)
        
