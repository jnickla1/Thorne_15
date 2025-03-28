import numpy as np
import pandas as pd
import os
#EBMKF_Nicklas as ekf
from netCDF4 import Dataset
import pdb

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])

def ta_smooth(orig_opt_depths,fill_value):
    wt_opt_depths = 1/(orig_opt_depths+9.7279)
    N = 30
    nwt_opt_depths=np.full(len(orig_opt_depths),fill_value)
    cN=int(np.ceil(N/2))
    fN=int(np.floor(N/2))
    for i in range((fN),(len(nwt_opt_depths)-1)):
        lasta=i+cN;firsta=i-fN
        nwt_opt_depths[i] = (np.sum(wt_opt_depths[(firsta):i+1])+ fill_value*(cN-1))/N
    #computing half-average - future is assumed to be the average 
    return(1/nwt_opt_depths-9.7279)

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

    
    if experiment_type == "historical":
        from . import EBMKF_Nicklas2a as ekf #only changed cloud feedback
        ekf.n_iters = 174 
        ekf.opt_depth=ta_smooth(ekf.data[:,3]*0.001,ekf.involcavg)
        means[0:ekf.n_iters], ses[0:ekf.n_iters],means2[0:ekf.n_iters],ses2[0:ekf.n_iters] = ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)
        return means - ekf.offset-preind_base, np.sqrt(np.abs(ses)),means2 -ekf.offset-preind_base,np.sqrt(np.abs(ses2))

    else:
        
        exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126_constVolc #
        new_iter=len(years)
        given_preind_base = np.mean(temperature[0:50])
        

        if (exp_attr[1]=='ESM1-2-LR'):
            from . import EBMKF_Nicklas2b as ekf #changed cloud feedback and thermal conductivity
            ekf.n_iters = 174 
            unf_new_opt_depth= np.full(new_iter, (1/ekf.involcavg-9.7279))
            unf_new_opt_depth[0:ekf.n_iters]=ekf.data[:,3]*0.001 #infill old optical depths based on existing historical AOD record, later ta-smoothed

           # new_opt_depth[ekf.n_iters:]=0.03
           #dont' change it after what's already been written for this case
            
            import xarray as xr
            ohca_later =Dataset(os.path.expanduser('~/')+"climate_data/ESM1-2-LR/opottempmint/"+exp_attr[2].lower()+"_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_l = ohca_later[:].__array__()
            ohca_ly = average_every_n(ohca_l[model_run,:], 12)
            ohca_earlier = Dataset(os.path.expanduser('~/')+"climate_data/ESM1-2-LR/opottempmint/historical_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_e = ohca_earlier[:].__array__()
            ohca_ey = average_every_n(ohca_e[model_run,:], 12)
            ohca_meas = np.concatenate((ohca_ey,ohca_ly))
            ohca_meas = ohca_meas - ohca_meas[0] #start at 0 #supplies full record

        
        elif (exp_attr[1]=='NorESM'):
            from . import EBMKF_Nicklas2c as ekf #changed cloud feedback and thermal conductivity
            ekf.n_iters = 174 
            unf_new_opt_depth= np.full(new_iter, (1/ekf.involcavg-9.7279))
            unf_new_opt_depth[0:ekf.n_iters]=ekf.data[:,3]*0.001 #infill old optical depths based on existing historical AOD record, later ta-smoothed
            ekf.n_iters = 174 
            ohca_meas=np.zeros(new_iter)
            ohca_meas[0:ekf.n_iters] = ekf.ocean_heat_measured/ekf.zJ_from_W
            
            import xarray as xr
            #ohca_earlier = Dataset(os.path.expanduser('~/')+"climate_data/NorESM_volc/OHCA/historicalVolc_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_earlier = Dataset(os.path.expanduser('~/')+"climate_data/NorESM_volc/OHCA/historicalVolc_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_e = ohca_earlier[:].__array__()
            #ohca_meas[(1980-1850):(2006-1850)] = average_every_n(ohca_e[model_run,:], 12) + ohca_meas[(1979-1850)]
            ohca_meas[(1980-1850):(2006-1850)] = ohca_e[(model_run if model_run< 37 else model_run-1) ,:] + 2*ohca_meas[(1979-1850)] - ohca_meas[(1978-1850)]

            aod_earlier = Dataset(os.path.expanduser('~/')+"climate_data/NorESM_volc/BethkeEtAl2017/historicalVolc_aod.nc", 'r').variables['__xarray_dataarray_variable__']
            aod_e = aod_earlier[:].__array__()
            unf_new_opt_depth[(1980-1850):(2006-1850)]= average_every_n(aod_e[model_run,:], 12)/1000

            if exp_attr[3]=='Volc':
                #opt_depth        
                aod_later = Dataset(os.path.expanduser('~/')+"/climate_data/NorESM_volc/BethkeEtAl2017/rcp45Volc_aod.nc", 'r').variables['__xarray_dataarray_variable__']
                aod_l = aod_later[:].__array__()
                unf_new_opt_depth[(2006-1850):(2100-1850)]= average_every_n(aod_l[model_run,:], 12)/1000
                #ohca
                ohca_later = Dataset(os.path.expanduser('~/')+"/climate_data/NorESM_volc/OHCA/rcp45Volc_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
                ohca_l = ohca_later[:].__array__()
                ohca_meas[(2006-1850):(2100-1850)]= ohca_l[model_run,:] + 2*ohca_meas[(2005-1850)]-ohca_meas[(2004-1850)]
            

            elif exp_attr[3]=='VolcConst':
                aod_later = Dataset(os.path.expanduser('~/')+"/climate_data/NorESM_volc/BethkeEtAl2017/rcp45VolcConst_partial20_aod.nc", 'r').variables['__xarray_dataarray_variable__']
                aod_l = aod_later[:].__array__()
                unf_new_opt_depth[(2006-1850):(2100-1850)]= average_every_n(aod_l[model_run%20,:], 12)/1000
                ohca_later = Dataset(os.path.expanduser('~/')+"/climate_data/NorESM_volc/OHCA/rcp45VolcConst_partial20_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
                ohca_l = ohca_later[:].__array__()
                ohca_meas[(2006-1850):(2100-1850)]=ohca_l[model_run if (model_run>2 and model_run<14) else model_run%14,:]  + 2*ohca_meas[(2005-1850)]-ohca_meas[(2004-1850)]      
            erf_data = pd.read_csv(os.path.expanduser(f"~/climate_data/SSP_inputdata/ERFs-Smith-ar6/ERF_SSP245_1750-2500.csv"))
            contrails = erf_data['contrails'][(1980-1750):(2100-1750)].values
            unf_new_opt_depth[(1980-1850):] = unf_new_opt_depth[(1980-1850):] + (contrails-0.015)/.18*.04  # get rid of gradual decline in the baseline over 21st century

                
        new_opt_depth = ta_smooth(unf_new_opt_depth,ekf.involcavg)
            
        if ekf.n_iters != new_iter:
            new_tsi = np.full(new_iter, ekf.sw_in)
            new_tsi = pd.read_csv(os.path.expanduser('~/')+"climate_data/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp119_1750-2500.csv")['solar'][100:(100+new_iter+1)].values
            new_tsi = new_tsi + ekf.sw_in-np.mean(new_tsi)

            new_R_tvar =np.full(new_iter, np.mean(ekf.R_tvar[150:174]))
            new_R_tvar[0:ekf.n_iters]=ekf.R_tvar[0:ekf.n_iters]
            

            new_Roc_tvar =np.full(new_iter, np.mean(ekf.Roc_tvar[50:174]))
            new_Roc_tvar[0:ekf.n_iters]= ekf.Roc_tvar[0:ekf.n_iters]
            
            data3 =  np.genfromtxt(open(os.path.expanduser('~/')+"climate_data/SSP_inputdata/KF6projectionSSP.csv", "rb"),dtype=float, delimiter=',')
            SSPnames=[126,434,245,370,585]
            if exp_attr[2]=='RCP45':
                find_case = 245
            else:
                find_case = int(exp_attr[2][3:])
            rcp = SSPnames.index(find_case)
            handoffyr = 1850+ekf.n_iters
            if (exp_attr[1]=='ESM1-2-LR'):
                new_Co2_df = pd.read_csv(open(os.path.expanduser('~/')+"climate_data/SSP_inputdata/eCO2_"+exp_attr[1]+"_"+exp_attr[2].lower()+".csv"),dtype=float, delimiter=',')
                new_lCo2 = np.log10(new_Co2_df['eCO2'].values)
            elif (exp_attr[1]=='NorESM'):
                from . import gen_eCO2
                if exp_attr[3]=='Volc':
                    erf_data = pd.read_csv(os.path.expanduser(f"~/climate_data/SSP_inputdata/ERF_NorESM_rcp45VolcConst.csv"))
                elif exp_attr[3]=='VolcConst':
                    erf_data = pd.read_csv(os.path.expanduser(f"~/climate_data/SSP_inputdata/ERFanthro_NorESM_rcp45Volc.csv"))
                model_outputlCo2 = gen_eCO2.calculate_equivalent_co2(erf_data['ERF_anthro'].values)
                new_lCo2 = np.concatenate((np.log10(ekf.data[:(1980-1850),2]), np.log10(ekf.data[(1981-1850),2] - model_outputlCo2[0]  + model_outputlCo2)))            #compute_update(xhatf[k-1],0,volcs[k-1],case[startyr-2015+k-1],caseA[startyr-2015+k-1]+1)
            #np.log10(data3[:,1+rcp]),data3[:,6+rcp]
            new_anthro_clouds = np.concatenate((ekf.anthro_clouds,data3[handoffyr-2015: 1850+new_iter-2015,6+rcp]+1))
            #breakpoint()
            
            ekf.opt_depth=new_opt_depth
            ekf.tsi = new_tsi
            ekf.R_tvar = new_R_tvar
            ekf.Roc_tvar=new_Roc_tvar
            ekf.lCo2=new_lCo2
            ekf.anthro_clouds = new_anthro_clouds       
            ekf.n_iters=new_iter
        new_observ = np.transpose(np.array([temperature-given_preind_base+preind_base + ekf.offset,ohca_meas/ekf.zJ_from_W ]))
        sz = (new_iter,2) # size of array is now different
        sz2d=(new_iter,2,2)
        ekf.xhat=np.zeros(sz)      # a posteri estimate of x
        ekf.P=np.zeros(sz2d)         # a posteri error estimate
        ekf.F=np.zeros(sz2d)         # state transitions
        ekf.xhatminus=np.zeros(sz) # a priori estimate of x
        ekf.Pminus=np.zeros(sz2d)    # a priori error estimate
        ekf.K=np.zeros(sz2d)         # gain or blending factor
        ekf.xhathat=np.zeros(sz)   # smoothed a priori estimate of x
        ekf.Phat=np.zeros(sz2d)      # smoothed posteri error estimate 
        ekf.Khat=np.zeros(sz2d)      # smoothed gain or blending factor
        ekf.Shat=np.zeros(sz2d)
        ekf.S=np.zeros(sz2d)
        ekf.xblind=np.zeros(sz)
        ekf.lml=np.zeros(sz)
        ekf.y=np.zeros(sz)
        ekf.qqy=np.zeros(sz)
        
        #breakpoint()

        means , ses ,means2 ,ses2 = ekf.ekf_run(new_observ,new_iter,retPs=3)
        return means - ekf.offset-preind_base + given_preind_base , np.sqrt(np.abs(ses)),means2 -ekf.offset-preind_base + given_preind_base,np.sqrt(np.abs(ses2))


   
    #(trailfilter_xhat,P2s,S2s,trailfilter_xhatm)=ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)
