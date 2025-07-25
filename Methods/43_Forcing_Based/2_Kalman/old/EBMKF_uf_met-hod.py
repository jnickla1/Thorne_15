import numpy as np
import pandas as pd
import config
from . import EBMKF_Nicklas as ekf
#EBMKF_Nicklas as ekf
from netCDF4 import Dataset
import pdb

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])

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

    ekf.n_iters = 174 
    if experiment_type == "historical":
        means[0:ekf.n_iters], ses[0:ekf.n_iters],means2[0:ekf.n_iters],ses2[0:ekf.n_iters] = ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)
        return means - ekf.offset-preind_base, np.sqrt(np.abs(ses)),means2 -ekf.offset-preind_base,np.sqrt(np.abs(ses2))

    else:
        exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126_constVolc #
        new_iter=len(years)
        given_preind_base = np.mean(temperature[0:50])

        

        if (exp_attr[1]=='ESM1-2-LR'):

            new_opt_depth= np.full(new_iter, (1/ekf.involcavg-9.7279))
            new_opt_depth[0:ekf.n_iters]=ekf.opt_depth
            
            import xarray as xr
            ohca_later =Dataset(config.CLIMATE_DATA_PATH+"/ESM1-2-LR/opottempmint/"+exp_attr[2].lower()+"_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_l = ohca_later[:].__array__()
            ohca_ly = average_every_n(ohca_l[model_run,:], 12)
            ohca_earlier = Dataset(config.CLIMATE_DATA_PATH+"/climate_data/ESM1-2-LR/opottempmint/historical_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_e = ohca_earlier[:].__array__()
            ohca_ey = average_every_n(ohca_e[model_run,:], 12)
            ohca_meas = np.concatenate((ohca_ey,ohca_ly+ohca_ey[-1]))

        
        elif (exp_attr[1]=='NorESM'):
            #TODO fill this in later
            #opt_depth
            ohca_meas=1
            

        new_tsi = np.full(new_iter, ekf.sw_in)
        new_tsi[0:ekf.n_iters] = ekf.data[:,8]

        new_R_tvar =np.full(new_iter, np.mean(ekf.R_tvar[-11:-1]))
        new_R_tvar[0:ekf.n_iters]=ekf.R_tvar
        

        new_Roc_tvar =np.full(new_iter, np.mean(ekf.Roc_tvar[-11:-1]))
        new_Roc_tvar[0:ekf.n_iters]= ekf.Roc_tvar
        
        data3 =  np.genfromtxt(open(config.CLIMATE_DATA_PATH+"/SSP_inputdata/KF6projectionSSP.csv", "rb"),dtype=float, delimiter=',')
        SSPnames=[126,434,245,370,585]
        if exp_attr[2]=='RCP45':
            find_case = 245
        else:
            find_case = int(exp_attr[2][3:])
        rcp = SSPnames.index(find_case)
        handoffyr = 1850+ekf.n_iters
        new_lCo2 = np.concatenate((ekf.lCo2, np.log10(data3[handoffyr-2015: 1850+new_iter-2015,1+rcp])))
        new_anthro_clouds = np.concatenate((ekf.anthro_clouds,data3[handoffyr-2015: 1850+new_iter-2015,6+rcp]+1))
        #compute_update(xhatf[k-1],0,volcs[k-1],case[startyr-2015+k-1],caseA[startyr-2015+k-1]+1)
        #np.log10(data3[:,1+rcp]),data3[:,6+rcp]
        
        #breakpoint()
        new_observ = np.transpose(np.array([temperature-given_preind_base+preind_base + ekf.offset,ohca_meas/ekf.zJ_from_W ]))
        ekf.opt_depth=new_opt_depth
        ekf.tsi = new_tsi
        ekf.R_tvar = new_R_tvar
        ekf.Roc_tvar=new_Roc_tvar
        ekf.lCo2=new_lCo2
        ekf.anthro_clouds = new_anthro_clouds       
        ekf.n_iters=new_iter
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
