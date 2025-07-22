import numpy as np
import pandas as pd
import config
#EBMKF_Nicklas as ekf
from netCDF4 import Dataset
import pdb
from scipy.optimize import minimize
import matplotlib.pyplot as plt


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
    from . import EBMKF_Nicklas3 as ekf
    
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
    n_iters_nofit = 125
    
    if experiment_type == "historical":
        #only changed cloud feedback
        tot_iters=174
        ekf.n_iters = n_iters_nofit
        raw_opt_depth=ekf.opt_depth.copy()
        wt_opt_depths = 1/(ekf.opt_depth+9.7279)
        N = 35
        nwt_opt_depths=np.empty(len(ekf.opt_depth)); nwt_opt_depths[:]=ekf.involcavg
        cN=int(np.ceil(N/2))
        fN=int(np.floor(N/2))
        for i in range((fN),(len(nwt_opt_depths)-1)):
            lasta=i+cN;firsta=i-fN
            nwt_opt_depths[i] = (np.sum(wt_opt_depths[(firsta):i+1])+ ekf.involcavg*(cN-1))/N
        #computing half-average - future is assumed to be the average 
        nopt_depths=(1/nwt_opt_depths-9.7279)
        ekf.opt_depth=nopt_depths
                
        means[0:n_iters_nofit], ses[0: n_iters_nofit],_,_ = ekf.ekf_run(ekf.observ, n_iters_nofit,retPs=3)

        initial_guess = [ekf.gad_prior_mean, ekf.fdbkA_prior_mean]
        #_,Pfirst_pass = ekf.ekf_run(ekf.observ,tot_iters,retPs=True)
        for n_iter in range(n_iters_nofit,tot_iters+1):
            ekf.n_iters = n_iter

            sz = (n_iter,2)
            N = 30
            z_run_means =np.empty(sz); z_run_means[:] = np.nan #ensure that this starts as an all-nan array
            cN=int(np.ceil(N/2)); fN=int(np.floor(N/2));
            for i in range((fN),(n_iter-cN+1)):
                lasta=i+cN;firsta=i-fN
                z_run_means[i,:] = np.mean(ekf.observ[firsta:lasta,:],axis=0);
            if n_iter%5==0:
                result = minimize(ekf.efk_reeval_run_likeli3, initial_guess, 
                    args=(z_run_means, n_iter, ekf.observ,raw_opt_depth), #np.sqrt(Pfirst_pass[:,0,0]), np.sqrt(Pfirst_pass[:,1,1])),
                    method='L-BFGS-B' , tol=0.01 )#, options = {'disp' : True} )# Gradient-based optimization method
                # Extract the MAP estimates
                initial_guess = result.x
                print(result.x)
            
            means_trial, ses_trial, means2_trial, ses2_trial = ekf.ekf_run(ekf.observ, n_iter,retPs=3)
            means[n_iter-1] = means_trial[n_iter-1]
            ses[n_iter-1] = ses_trial[n_iter-1]
            if n_iter ==tot_iters:
                means2[0:n_iter]= means2_trial
                ses2[0:n_iter] = ses2_trial
                        
        return means - ekf.offset-preind_base, np.sqrt(np.abs(ses) ),means2 -ekf.offset-preind_base,np.sqrt(np.abs(ses2 ) )

    else:
        exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126_constVolc #
        new_iter=len(years)
        given_preind_base = np.mean(temperature[0:75]) #just a bit longer than 50 yrs

        if (exp_attr[1]=='ESM1-2-LR'):

            new_opt_depth= np.full(new_iter, (1/ekf.involcavg-9.7279))
            raw_opt_depth=new_opt_depth.copy()
            #new_opt_depth_inst= np.full(new_iter, (1/ekf.involcavg-9.7279))
            #new_opt_depth_inst[0:ekf.n_iters]=ekf.opt_depth

            wt_opt_depths = 1/(ekf.opt_depth+9.7279)
            N = 30
            nwt_opt_depths=np.empty(len(ekf.opt_depth)); nwt_opt_depths[:]=ekf.involcavg
            cN=int(np.ceil(N/2))
            fN=int(np.floor(N/2))
            for i in range((fN),(len(nwt_opt_depths)-1)):
                lasta=i+cN;firsta=i-fN
                nwt_opt_depths[i] = (np.sum(wt_opt_depths[(firsta):i+1])+ ekf.involcavg*(cN-1))/N
        #computing half-average - future is assumed to be the average 
            nopt_depths=(1/nwt_opt_depths-9.7279)
    
            new_opt_depth[0:ekf.n_iters]=nopt_depths
            raw_opt_depth[0:ekf.n_iters]=ekf.opt_depth
            
            import xarray as xr
            ohca_later =Dataset(config.CLIMATE_DATA_PATH+"/ESM1-2-LR/opottempmint/"+exp_attr[2].lower()+"_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_l = ohca_later[:].__array__()
            ohca_ly = average_every_n(ohca_l[model_run,:], 12)
            ohca_earlier = Dataset(config.CLIMATE_DATA_PATH+"/ESM1-2-LR/opottempmint/historical_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_e = ohca_earlier[:].__array__()
            ohca_ey = average_every_n(ohca_e[model_run,:], 12)
            ohca_meas = np.concatenate((ohca_ey,ohca_ly ))
            ohca_meas = ohca_meas - ohca_meas[0] #start at 0

        
        elif (exp_attr[1]=='NorESM'):
            #TODO fill this in later once we get a working method
            return empser,empser,empser,empser
            

        if ekf.n_iters != new_iter:
            new_tsi = np.full(new_iter, ekf.sw_in)
            new_tsi = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp119_1750-2500.csv")['solar'][100:(100+new_iter+1)].values
            new_tsi = new_tsi + ekf.sw_in-np.mean(new_tsi)

            new_R_tvar =np.full(new_iter, np.mean(ekf.R_tvar[-11:-1]))
            new_R_tvar[0:ekf.n_iters]=ekf.R_tvar
            

            new_Roc_tvar =np.full(new_iter, np.mean(ekf.Roc_tvar[-11:-1]))
            new_Roc_tvar[0:ekf.n_iters]= ekf.Roc_tvar
            
            data3 =  np.genfromtxt(open(oconfig.CLIMATE_DATA_PATH+"/SSP_inputdata/KF6projectionSSP.csv", "rb"),dtype=float, delimiter=',')
            SSPnames=[126,434,245,370,585]
            if exp_attr[2]=='RCP45':
                find_case = 245
            else:
                find_case = int(exp_attr[2][3:])
            rcp = SSPnames.index(find_case)
            handoffyr = 1850+ekf.n_iters
            new_Co2_df = pd.read_csv(open(config.CLIMATE_DATA_PATH+"/SSP_inputdata/eCO2_"+exp_attr[1]+"_"+exp_attr[2].lower()+".csv"),dtype=float, delimiter=',')
            new_lCo2 = np.log10(new_Co2_df['eCO2'].values)
            #new_lCo2 = np.concatenate((ekf.lCo2, np.log10(data3[handoffyr-2015: 1850+new_iter-2015,1+rcp])))
            new_anthro_clouds = np.concatenate((ekf.anthro_clouds,data3[handoffyr-2015: 1850+new_iter-2015,6+rcp]+1))
            #compute_update(xhatf[k-1],0,volcs[k-1],case[startyr-2015+k-1],caseA[startyr-2015+k-1]+1)
            #np.log10(data3[:,1+rcp]),data3[:,6+rcp]
            
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
        
        

        n_iters_nofit = 125
        tot_iters = new_iter #251
        #breakpoint()

        means[0:n_iters_nofit], ses[0: n_iters_nofit],_,_ = ekf.ekf_run(new_observ, n_iters_nofit,retPs=3)

        initial_guess = [ekf.gad_prior_mean, ekf.fdbkA_prior_mean]
        _,Pfirst_pass = ekf.ekf_run(new_observ,tot_iters,retPs=True)
        for n_iter in range(n_iters_nofit,tot_iters+1):
            ekf.n_iters = n_iter
            sz = (n_iter,2)
            N = 30
            z_run_means =np.empty(sz); z_run_means[:] = np.nan #ensure that this starts as an all-nan array
            cN=int(np.ceil(N/2)); fN=int(np.floor(N/2));
            for i in range((fN),(n_iter-cN+1)):
                lasta=i+cN;firsta=i-fN
                z_run_means[i,:] = np.mean(new_observ[firsta:lasta,:],axis=0);

            if (n_iter%5==0):
                result = minimize(ekf.efk_reeval_run_likeli3, initial_guess, 
                    args=(z_run_means, n_iter,new_observ,raw_opt_depth), #,new_opt_depth_inst
                    #args=(z_run_means, n_iter, np.sqrt(Pfirst_pass[:,0,0]), np.sqrt(Pfirst_pass[:,1,1])),  
                    method='L-BFGS-B' , tol=0.01 ,bounds = ((0.2, 1.75), (-.5, 1.5)))#, options = {'disp' : True} )# Gradient-based optimization method
                # Extract the MAP estimates
                initial_guess = result.x
                print(result.x,n_iter+1850)

            if (False):
                means_trial, ses_trial, means2_trial, ses2_trial = ekf.ekf_run(new_observ, n_iter,retPs=3,plottin=True)
                ekf.gad=1.2 ; ekf.fdbkA= 0.3# ekf.fdbkA_prior_mean
                ekf.precompute_coeffs(False)
                ekf.ekf_run(new_observ, n_iter,retPs=3,plottin=True)
                plt.show()

            else:
                means_trial, ses_trial, means2_trial, ses2_trial = ekf.ekf_run(new_observ, n_iter,retPs=3)
            
            means[n_iter-1] = means_trial[n_iter-1]
            ses[n_iter-1] = ses_trial[n_iter-1]
            if n_iter ==tot_iters:
                means2[0:n_iter]= means2_trial
                ses2[0:n_iter] = ses2_trial

      
        return means - ekf.offset-preind_base + given_preind_base , np.sqrt(np.abs(ses) + np.square(temps_1std )  -
            0*np.square(temps_1std[-1] )), means2 -ekf.offset-preind_base + given_preind_base,np.sqrt(np.abs(ses2) + np.square(temps_1std)- 0* np.square(temps_1std[-1] ))


   
    #(trailfilter_xhat,P2s,S2s,trailfilter_xhatm)=ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)
