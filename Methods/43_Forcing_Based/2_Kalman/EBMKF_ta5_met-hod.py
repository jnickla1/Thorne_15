import numpy as np
import pandas as pd
import os
#EBMKF_Nicklas as ekf
from netCDF4 import Dataset
import pdb
import xarray as xr

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


def ta_smooth(orig_opt_depths,fill_value,optical=True):
    if(optical):
        wt_opt_depths = 1/(orig_opt_depths+9.7279)
    else:
        wt_opt_depths =orig_opt_depths
    N = 30
    nwt_opt_depths=np.full(len(orig_opt_depths),fill_value)
    cN=int(np.ceil(N/2))
    fN=int(np.floor(N/2))
    for i in range((fN),(len(nwt_opt_depths)-1)):
        lasta=i+cN;firsta=i-fN
        nwt_opt_depths[i] = (np.sum(wt_opt_depths[(firsta):i+1])+ fill_value*(cN-1))/N
    #computing half-average - future is assumed to be the average
    if(optical):
        return(1/nwt_opt_depths-9.7279)
    else:
        return(nwt_opt_depths)

def run_method(years, temperature, uncert, model_run, experiment_type):
    #temperature has 2024, dont have all forcings for this yet - must update ekf.n_iters
    data_orig = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data_orig.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50])
    
    #ekf.temps = temperature[0:ekf.n_iters] +ekf.offset
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4
    
    #ekf.R_tvar=np.square(temps_1std[0:ekf.n_iters])
    #ekf.fdbkA = 0.35 #original paper value for this estimate
    #ekf.precompute_coeffs(False)
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    means2 = empser.copy()
    ses2 = empser.copy()
    from . import EBMKF_Nicklas5 as ekf #only change cloud feedback, but dynamically
    
    if experiment_type == "historical":
        
        ekf.n_iters = 174 
        ekf.opt_depth=ta_smooth(ekf.data[:,3]*0.001,ekf.involcavg)
        means[0:ekf.n_iters], ses[0:ekf.n_iters],means2[0:ekf.n_iters],ses2[0:ekf.n_iters] = ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)
        return means - ekf.offset-preind_base, np.sqrt(np.abs(ses)),means2 -ekf.offset-preind_base,np.sqrt(np.abs(ses2))

    else:
        
        exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126_constVolc #
        new_iter=len(years)
        given_preind_base = np.mean(temperature[0:50])
        

        if (exp_attr[1]=='ESM1-2-LR'):
            ekf.n_iters = 174 
            unf_new_opt_depth= np.full(new_iter, (1/ekf.involcavg-9.7279))
            unf_new_opt_depth[0:ekf.n_iters]=ekf.data[:,3]*0.001 #infill old optical depths based on existing historical AOD record, later ta-smoothed
            
            ohca_later =Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/ESM1-2-LR/opottempmint/"+exp_attr[2].lower()+"_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_l = ohca_later[:].__array__()
            ohca_ly = average_every_n(ohca_l[model_run,:], 12)
            ohca_earlier = Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/ESM1-2-LR/opottempmint/historical_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_e = ohca_earlier[:].__array__()
            ohca_ey = average_every_n(ohca_e[model_run,:], 12)
            ohca_meas = np.concatenate((ohca_ey,ohca_ly))
            ohca_meas = ohca_meas - ohca_meas[0] #start at 0 #supplies full record

            rt_later =Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/ESM1-2-LR/rt/global_rt_"+exp_attr[2].lower()+".nc", 'r').variables['rt']
            rt_l = rt_later[:].__array__()
            rt_ly = average_every_n(rt_l[model_run,:], 12)
            rt_earlier = Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/ESM1-2-LR/rt/global_rt_historical.nc", 'r').variables['rt']
            rt_e = rt_earlier[:].__array__()
            rt_ey = average_every_n(rt_e[model_run,:], 12)
            new_rtTOA = np.concatenate((rt_ey,rt_ly))

        
        elif (exp_attr[1]=='NorESM'):
            # need to change thermal conductivity ??
            ekf.n_iters = 174 
            unf_new_opt_depth= np.full(new_iter, (1/ekf.involcavg-9.7279))
            #unf_new_opt_depth[0:ekf.n_iters]=ekf.data[:,3]*0.001 #infill old optical depths based on existing historical AOD record, later ta-smoothed DONE - replace with NorESM simulation
            aod_spinup = Dataset(os.path.expanduser('~/')+"/data/jnickla1/climate_data/NorESM_volc/NorESM1-M-historical/hist_aod.nc", 'r').variables['__xarray_dataarray_variable__']
            aod_s = aod_spinup[:].__array__()
            long_past_index = (gen_orig_number(model_run,60) // 20)
            unf_new_opt_depth[(1850-1850):(2006-1850)]= average_every_n(aod_s[long_past_index,:], 12)/1000
            
            ohca_meas=np.zeros(new_iter)
            #ohca_meas[0:ekf.n_iters] = ekf.ocean_heat_measured/ekf.zJ_from_W   #DONE - replace with NorESM simulation
            ohca_spinup = Dataset(os.path.expanduser('~/')+"/data/jnickla1/climate_data/NorESM_volc/NorESM1-M-historical/hist_ohca_mon.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_s = ohca_spinup[:].__array__()
            ohca_s_yr = average_every_n(ohca_s[long_past_index,:], 12)
            ohca_meas[(1850-1850):(2006-1850)]= ohca_s_yr - ohca_s_yr[0]
            
            #ohca_earlier = Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/NorESM_volc/OHCA/historicalVolc_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_earlier = Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/NorESM_volc/OHCA/historicalVolc_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_e = ohca_earlier[:].__array__()
            #ohca_meas[(1980-1850):(2006-1850)] = average_every_n(ohca_e[model_run,:], 12) + ohca_meas[(1979-1850)]
            ohca_meas[(1980-1850):(2006-1850)] = ohca_e[(model_run if model_run< 37 else model_run-1) ,:] + 2*ohca_meas[(1979-1850)] - ohca_meas[(1978-1850)] #missing one model_run ohca, so shift them

            aod_earlier = Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/NorESM_volc/BethkeEtAl2017/historicalVolc_aod.nc", 'r').variables['__xarray_dataarray_variable__']
            aod_e = aod_earlier[:].__array__()
            unf_new_opt_depth[(1980-1850):(2006-1850)]= average_every_n(aod_e[model_run,:], 12)/1000
            
            new_rtTOA = np.zeros(new_iter)
            rt_spinup = Dataset(os.path.expanduser('~/')+"/data/jnickla1/climate_data/NorESM_volc/NorESM1-M-historical/rtmt/rtmt_NorESM1-M_historical.nc", 'r').variables['rtmt']
            rt_s = rt_spinup[:].__array__()
            new_rtTOA[(1850-1850):(2006-1850)]=average_every_n(rt_s[long_past_index,:], 12)
            rt_earlier = Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/NorESM_volc/BethkeEtAl2017/historicalVolc_rtmt/rtmt_NorESM1-M_historicalVolc.nc", 'r').variables['rtmt']
            rt_e = rt_earlier[:].__array__()
            new_rtTOA[(1980-1850):(2006-1850)]=average_every_n(rt_e[model_run,:], 12)
            if exp_attr[3]=='Volc':
                #opt_depth        
                aod_later = Dataset(os.path.expanduser('~/')+"/data/jnickla1/climate_data/NorESM_volc/BethkeEtAl2017/rcp45Volc_aod.nc", 'r').variables['__xarray_dataarray_variable__']
                aod_l = aod_later[:].__array__()
                unf_new_opt_depth[(2006-1850):(2100-1850)]= average_every_n(aod_l[model_run,:], 12)/1000
                #ohca
                ohca_later = Dataset(os.path.expanduser('~/')+"/data/jnickla1/climate_data/NorESM_volc/OHCA/rcp45Volc_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
                ohca_l = ohca_later[:].__array__()
                ohca_meas[(2006-1850):(2100-1850)]= ohca_l[model_run,:] + 2*ohca_meas[(2005-1850)]-ohca_meas[(2004-1850)]
                rt_later =Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/NorESM_volc/BethkeEtAl2017/rcp45"+exp_attr[3]+"_rtmt/rtmt_NorESM1-M_rcp45"+exp_attr[3]+".nc", 'r').variables['rtmt']
                rt_l = rt_later[:].__array__()
                new_rtTOA[(2006-1850):(2100-1850)]= average_every_n(rt_l[model_run,:], 12)
            

            elif exp_attr[3]=='VolcConst':
                aod_later = Dataset(os.path.expanduser('~/')+"/data/jnickla1/climate_data/NorESM_volc/BethkeEtAl2017/rcp45VolcConst_partial20_aod.nc", 'r').variables['__xarray_dataarray_variable__']
                aod_l = aod_later[:].__array__()
                unf_new_opt_depth[(2006-1850):(2100-1850)]= average_every_n(aod_l[model_run%20,:], 12)/1000
                ohca_later = Dataset(os.path.expanduser('~/')+"/data/jnickla1/climate_data/NorESM_volc/OHCA/rcp45VolcConst_partial20_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
                ohca_l = ohca_later[:].__array__()
                ohca_meas[(2006-1850):(2100-1850)]=ohca_l[model_run if (model_run>2 and model_run<14) else model_run%14,:]  + 2*ohca_meas[(2005-1850)]-ohca_meas[(2004-1850)]
                rt_later =Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/NorESM_volc/BethkeEtAl2017/rcp45"+exp_attr[3]+"_rtmt/rtmt_NorESM1-M_rcp45"+exp_attr[3]+"_20.nc", 'r').variables['rtmt']
                rt_l = rt_later[:].__array__()
                new_rtTOA[(2006-1850):(2100-1850)]= average_every_n(rt_l[model_run%20,:], 12)
                
            erf_data = pd.read_csv(os.path.expanduser(f"~/data/jnickla1/climate_data/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv"))
            contrails = erf_data['contrails'][(1980-1750):(2100-1750)].values
            unf_new_opt_depth[(1980-1850):] = unf_new_opt_depth[(1980-1850):] + (contrails-0.015)/.18*.04  # get rid of gradual decline in the baseline over 21st century


                
        new_opt_depth = ta_smooth(unf_new_opt_depth,ekf.involcavg)
        temps = temperature-given_preind_base+preind_base + ekf.offset #ensure this matches expected starting temperature in 1850

        frac_blocked = 1-(5.5518/(unf_new_opt_depth+9.9735)) #(9.068/(comb_recAOD_cleaned+9.7279))
        erf_data_solar = pd.read_csv(os.path.expanduser(f"~/data/jnickla1/climate_data/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv"))['solar'].to_numpy()[(1850-1750):]
        solar_full = erf_data_solar[:len(temps)]+ 340.4099428 - 0.108214
        tot_volc_erf = - solar_full*frac_blocked + 151.22
        erf_trapez=(1*tot_volc_erf[0:-2]+2*tot_volc_erf[1:-1]+1*tot_volc_erf[2:] )/4
        temps[2:] = temps[2:] - erf_trapez/(ekf.heatCp - ekf.Cs)
        ohca_meas[2:] = ohca_meas[2:] - (erf_trapez/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*200 #not perfect but good for now, also unclear where the 200 comes from
            
        if ekf.n_iters != new_iter:
            new_tsi = np.full(new_iter, ekf.sw_in)
            new_tsi = pd.read_csv(os.path.expanduser('~/')+"data/jnickla1/climate_data/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp119_1750-2500.csv")['solar'][100:(100+new_iter+1)].values
            new_tsi = new_tsi + ekf.sw_in-np.mean(new_tsi)

            new_R_tvar =np.full(new_iter, np.mean(ekf.R_tvar[150:174]))
            new_R_tvar[0:ekf.n_iters]=ekf.R_tvar[0:ekf.n_iters]
            

            new_Roc_tvar =np.full(new_iter, np.mean(ekf.Roc_tvar[50:174]))
            new_Roc_tvar[0:ekf.n_iters]= ekf.Roc_tvar[0:ekf.n_iters]
            
            data3 =  np.genfromtxt(open(os.path.expanduser('~/')+"data/jnickla1/climate_data/SSP_inputdata/KF6projectionSSP.csv", "rb"),dtype=float, delimiter=',')
            SSPnames=[126,434,245,370,585]
            if exp_attr[2]=='RCP45':
                find_case = 245
            else:
                find_case = int(exp_attr[2][3:])
            rcp = SSPnames.index(find_case)
            handoffyr = 1850+ekf.n_iters
            if (exp_attr[1]=='ESM1-2-LR'):
                new_Co2_df = pd.read_csv(open(os.path.expanduser('~/')+"data/jnickla1/climate_data/SSP_inputdata/eCO2_"+exp_attr[1]+"_"+exp_attr[2].lower()+".csv"),dtype=float, delimiter=',')
                new_lCo2 = np.log10(new_Co2_df['eCO2'].values)
                erf_data = pd.read_csv(os.path.expanduser(f"~/data/jnickla1/climate_data/SSP_inputdata/ERF_"+exp_attr[1]+"_"+exp_attr[2].lower()+".csv"))
                erf_natvolc = erf_data['ERF_natural']
            elif (exp_attr[1]=='NorESM'):
                from . import gen_eCO2
                if exp_attr[3]=='Volc':
                    erf_data = pd.read_csv(os.path.expanduser(f"~/data/jnickla1/Thorne_15_codefigurestats/Methods/43_Forcing_Based/3_Human_Induced/GWI_data/ERFanthro_NorESM_rcp45Volc_full.csv"))
                    #erf_data = pd.read_csv(os.path.expanduser(f"~/data/jnickla1/climate_data/SSP_inputdata/ERF_NorESM_rcp45Volc.csv"))
                    erf_natvolc = pd.read_csv(os.path.expanduser(f"~/data/jnickla1/Thorne_15_codefigurestats/Methods/43_Forcing_Based/3_Human_Induced/GWI_data/ERFnatural_NorESM_rcp45Volc_full.csv")).to_numpy()[:,model_run+1]
                elif exp_attr[3]=='VolcConst':
                    erf_data = pd.read_csv(os.path.expanduser(f"~/data/jnickla1/Thorne_15_codefigurestats/Methods/43_Forcing_Based/3_Human_Induced/GWI_data/ERF_NorESM_rcp45VolcConst_full.csv"))
                    erf_natvolc = erf_data['ERF_natural']
                model_outputlCo2 = gen_eCO2.calculate_equivalent_co2(erf_data['ERF_anthro'].values)
                new_lCo2 = np.log10(model_outputlCo2)
                
                
                #np.concatenate((np.log10(ekf.data[:(1980-1850),2]), np.log10(ekf.data[(1981-1850),2] - model_outputlCo2[0]  + model_outputlCo2)))            #compute_update(xhatf[k-1],0,volcs[k-1],case[startyr-2015+k-1],caseA[startyr-2015+k-1]+1)
            #np.log10(data3[:,1+rcp]),data3[:,6+rcp]
            new_anthro_clouds = np.concatenate((ekf.anthro_clouds,data3[handoffyr-2015: 1850+new_iter-2015,6+rcp]+1))
            #breakpoint()
            
            ekf.opt_depth=new_opt_depth
            ekf.tsi = new_tsi
            ekf.R_tvar = new_R_tvar
            ekf.Roc_tvar=new_Roc_tvar
            ekf.lCo2=new_lCo2
            ekf.anthro_clouds = new_anthro_clouds
            #TOA_means_artif = ta_smooth(new_rtTOA,0.00001,optical=False) + (new_rtTOA - erf_natvolc)/2 #   + np.mean(erf_volc) #half is the prior 30 yrs, half is what we expect the TOA balance to be removing this year's volcanism
            ekf.T02= np.mean(temps[(2000-1850):(2005-1850)]) #ensure critical T02 temperature is correct
            ekf.precompute_coeffs(False)
            TOA_meas_artif0 =  new_rtTOA - (erf_natvolc - erf_data_solar[:len(new_rtTOA )]) * .75  #+ np.mean(erf_volc)  to be removing this year's volcanism
            TOA_meas_artif = TOA_meas_artif0 - np.mean(TOA_meas_artif0[0:70]) #this should be balanced at the start
            
            #ohca_meas = ohca_meas - np.mean(new_rtTOA[0:70])*ekf.zJ_from_W*np.arange(new_iter) #correct for nonequilibrium drift in OHCA - doesn't work
            #Cd= 136.5*deepdepth/1000
            #np.sum(TOA_meas_artif)
            print("CONSISTENCY CHECK")
            
            surf_heat = (temps[-1]-np.mean(temps[0:70]))/(ekf.heatCp-ekf.Cs)
            print(f"OHCA + surf = {ohca_meas[-1]/ekf.zJ_from_W + surf_heat}")
            print(f"sum TOA = {np.sum(TOA_meas_artif)}")
            print("Error:")
            print((np.sum(TOA_meas_artif)-ohca_meas[-1]/ekf.zJ_from_W - surf_heat)/np.sum(TOA_meas_artif))
            #51.969 vs 52.819 - pretty close for historical record
            ohca_scale = (np.sum(TOA_meas_artif)- surf_heat)/ohca_meas[-1]*ekf.zJ_from_W
            print("OHCA_scale:")
            print(ohca_scale)

            #now fix the deep ocean depth / gad which means fitting values for blind model
            #preferably run just up to 2020s
            
            #ekf.precompute_coeffs(False)
            #ekf.n_iters=new_iter

        
        #breakpoint()
        #new_observ = np.transpose(np.array([temperature-given_preind_base+preind_base + ekf.offset,ohca_meas/ekf.zJ_from_W ]))
        ekf.TOA_meas_artif_var = np.full(new_iter,ekf.TOA_crop_var/4) #just put constant variance on this
        ndim =5
        new_observ = np.array([[temps,ohca_meas/ekf.zJ_from_W*ohca_scale ,TOA_meas_artif]]).T
        sz = (new_iter,ndim,1) # size of array
        ekf.sz=sz
        sz2d=(new_iter,ndim,ndim)
        ekf.sz2d=sz2d
        ekf.xhat=np.zeros(sz)      # a posteri estimate of x
        ekf.P=np.zeros(sz2d)         # a posteri error estimate
        ekf.F=np.zeros(sz2d)         # state transitions
        ekf.xhatminus=np.zeros(sz) # a priori estimate of x
        ekf.Pminus=np.zeros(sz2d)    # a priori error estimate
        ekf.K=np.zeros((new_iter,ndim,3))         # gain or blending factor
        ekf.xhathat=np.zeros(sz)   # smoothed a priori estimate of x
        ekf.Phat=np.zeros(sz2d)      # smoothed posteri error estimate 
        ekf.Khat=np.zeros((new_iter,ndim,ndim))      # smoothed gain or blending factor
        ekf.Shat=np.zeros((new_iter,3,3))
        ekf.S=np.zeros((new_iter,3,3))
        ekf.xblind=np.zeros((new_iter,2,1))
        ekf.lml=np.zeros(sz)
        ekf.lsml=0
        ekf.y=np.zeros((new_iter,3,1))
        ekf.qqy=np.zeros((new_iter,3))
        
        #breakpoint()

        means , ses ,means2 ,ses2 = ekf.ekf_run(new_observ,new_iter,retPs=3)
        return means - ekf.offset-preind_base + given_preind_base , np.sqrt(np.abs(ses)),means2 -ekf.offset-preind_base + given_preind_base,np.sqrt(np.abs(ses2))


   
    #(trailfilter_xhat,P2s,S2s,trailfilter_xhatm)=ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)
