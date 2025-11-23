import numpy as np
import pandas as pd
import config
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


def ta_smooth(orig_opt_depths,fill_value,optical=True, N=30):
    if(optical):
        wt_opt_depths = 1/(orig_opt_depths+9.7279)
    else:
        wt_opt_depths =orig_opt_depths
    nwt_opt_depths=np.full(len(orig_opt_depths),float(fill_value))
    cN=int(np.ceil(N/2))
    fN=int(np.floor(N/2))
    for i in range((fN),(len(nwt_opt_depths))):
        lasta=i+cN;firsta=i-fN
        nwt_opt_depths[i] = (np.sum(wt_opt_depths[(firsta):i+1])+ fill_value*(cN-1))/float(N)

    #computing half-average - future is assumed to be the average
    if(optical):
        return(1/nwt_opt_depths-9.7279)
    else:
        return(nwt_opt_depths)


from . import EBMKF_Nicklas4 as ekf #only change cloud feedback, but dynamically
if not hasattr(ekf, "_orig_state_saved"):
    ekf._orig_state_saved = True
    # Arrays
    ekf._orig_Q = ekf.Q.copy()
    ekf._orig_cM2 = ekf.cM2.copy()
    ekf._orig_opt_depth = ekf.opt_depth.copy()
    ekf._orig_tsi = ekf.tsi.copy()
    ekf._orig_R_tvar = ekf.R_tvar.copy()
    ekf._orig_Roc_tvar = ekf.Roc_tvar.copy()
    ekf._orig_TOA_meas_artif_var = ekf.TOA_meas_artif_var.copy()
    ekf._orig_lCo2 = ekf.lCo2.copy()
    ekf._orig_anthro_clouds = ekf.anthro_clouds.copy()
    # Scalars
    ekf._orig_dfrA_float_var = ekf.dfrA_float_var
    ekf._orig_gad = ekf.gad

def restore_baselines():
    ekf.Q = ekf._orig_Q.copy()
    ekf.cM2 = ekf._orig_cM2.copy()
    ekf.opt_depth = ekf._orig_opt_depth.copy()
    ekf.tsi = ekf._orig_tsi.copy()
    ekf.R_tvar = ekf._orig_R_tvar.copy()
    ekf.Roc_tvar = ekf._orig_Roc_tvar.copy()
    ekf.TOA_meas_artif_var = ekf._orig_TOA_meas_artif_var.copy()
    ekf.lCo2 = ekf._orig_lCo2.copy()
    ekf.anthro_clouds = ekf._orig_anthro_clouds.copy()
    ekf.dfrA_float_var = ekf._orig_dfrA_float_var
    ekf.gad = ekf._orig_gad
    

def run_method(years, temperature0, uncert, model_run, experiment_type):
    temperature = temperature0.astype(np.float64)
    #temperature has 2024, dont have all forcings for this yet - must update ekf.n_iters
    data_orig = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data_orig.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50])
    
    #ekf.temps = temperature[0:ekf.n_iters] +ekf.offset
    (temps_CIl, temps_CIu) = uncert
    temps_CIl = temps_CIl.astype(np.float64)
    temps_CIu = temps_CIu.astype(np.float64)
    temps_1std = (temps_CIu - temps_CIl) / 4
    
    #ekf.R_tvar=np.square(temps_1std[0:ekf.n_iters])
    #ekf.fdbkA = 0.35 #original paper value for this estimate
    #ekf.precompute_coeffs(False)
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    means2 = empser.copy()
    ses2 = empser.copy()


    exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126_constVolc #   
    
    if experiment_type == "historical":
        restore_baselines()
        ekf.n_iters = 174
        unf_new_opt_depth = ekf.data[:,3]*0.001
        ekf.opt_depth=ta_smooth(unf_new_opt_depth,(1/(unf_new_opt_depth[20]+9.7279)))
        given_preind_base = np.mean(temperature[0:50])
        temps = temperature-given_preind_base+preind_base + ekf.offset #ensure this matches expected starting temperature in 1850
        frac_blocked = 1-(5.5518/(unf_new_opt_depth+9.9735)) #(9.068/(comb_recAOD_cleaned+9.7279))
        erf_data= pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")
        erf_data_solar = erf_data['solar'].to_numpy()[(1850-1750):]
        solar_full = erf_data_solar[:ekf.n_iters]+ 340.4099428 - 0.108214
        tot_volc_erf = (- solar_full*frac_blocked + 151.22)
        contrails = erf_data['contrails'][(1850-1750):(2024-1750)].values
        TOA_correction = (tot_volc_erf - np.mean(tot_volc_erf[90:100])+ contrails)
        #TOA_correction = TOA_correction #- np.mean(TOA_correction)*2/3
        #TOA_correction[0:75] = TOA_correction[0:75] + 0.2
        
        TOA_correction[np.where(TOA_correction> -0.25)] = 0
        #breakpoint()
        #clamp down some of the uncertainties so that the KF pays most attention to the temperatures
        ekf.Q = ekf.Q * 0.25
        ekf.cM2 = ekf.cM2 * 0.5
        ekf.dfrA_float_var = ekf.dfrA_float_var *.25
        #now this is somewhat flat and only makes corrections when there is a major volcanic eruption
        erf_trapez=(1*TOA_correction[0:-2]+2*TOA_correction[1:-1]+0*TOA_correction[2:] )/3 #want to snap back as quickly as possible
        ohca_meas = ekf.oc_meas*ekf.zJ_from_W
        tsi_orig = ekf.tsi.copy()
        ekf.tsi = 2* ta_smooth(ekf.tsi,ekf.sw_in,optical=False,N=34) - ekf.sw_in  #np.full(np.shape(ekf.tsi),np.mean(ekf.tsi))
        #makes centered mean
        ekf.tsi[16:] = 2*ekf.tsi[16:]-ekf.tsi[:-16]
        #projecting forward into future
        
        #temps[2:-1] = temps[2:-1] - erf_trapez/(ekf.heatCp - ekf.Cs)
        #ohca_meas[2:] = ohca_meas[2:] - (erf_trapez/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*50 #not perfect but good for now, also unclear where the 50 comes from
        
        temps[9:ekf.n_iters] = temps[9:ekf.n_iters] - erf_trapez[:-7]/(ekf.heatCp - ekf.Cs)/2 #surface temp of 20 yr mean still depressed
        ohca_meas[9:ekf.n_iters] = ohca_meas[9:ekf.n_iters] - (erf_trapez[:-7]/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*50/2
        
        TOA_meas_artif1 =  ekf.TOA_meas_artif - TOA_correction * .75
        new_observ = np.array([[temps[:ekf.n_iters],ohca_meas/ekf.zJ_from_W ,TOA_meas_artif1]]).T
        #temps = temperature-given_preind_base+preind_base + ekf.offset
        #new_observ = np.array([[temps[:ekf.n_iters],ekf.oc_meas,ekf.TOA_meas_artif]]).T
        ekf.dfrA_float_var = ekf.dfrA_float_var/40
        means[0:ekf.n_iters], ses[0:ekf.n_iters],means2[0:ekf.n_iters],ses2[0:ekf.n_iters] = ekf.ekf_run(new_observ,ekf.n_iters,retPs=3)
        #means[0:ekf.n_iters] = means[0:ekf.n_iters] + (tsi_orig-ekf.tsi)/(ekf.heatCp) #adding TSI as a final correction doesnt work
        return means +given_preind_base - ekf.offset-preind_base, np.sqrt(np.abs(ses))*1.2,means2 -ekf.offset-preind_base,np.sqrt(np.abs(ses2))


    elif (exp_attr[0]=="histens"):
        restore_baselines()
        ekf.n_iters = len(years)
        new_iter=len(years)
        #keeping optical depth (volcanoes) and solar forcing the same
        unf_new_opt_depth0 = ekf.data[:,3]*0.001
        #add one more to the opt_depth - haven't yet worked out a new data source
        unf_new_opt_depth = np.append(unf_new_opt_depth0,[np.average(unf_new_opt_depth0[-10:])])

        ekf.opt_depth=ta_smooth(unf_new_opt_depth,ekf.involcavg)
        given_preind_base = np.mean(temperature[0:50])

        temps = temperature-given_preind_base+preind_base + ekf.offset #ensure this matches expected starting temperature in 1850, same baseline as HadCRUT
        frac_blocked = 1-(5.5518/(unf_new_opt_depth+9.9735)) #(9.068/(comb_recAOD_cleaned+9.7279))
        erf_data= pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")
        erf_data_solar = erf_data['solar'].to_numpy()[(1850-1750):]
        solar_full = erf_data_solar[:ekf.n_iters]+ 340.4099428 - 0.108214
        tot_volc_erf = (- solar_full*frac_blocked + 151.22)
        contrails = erf_data['contrails'][(1850-1750):(2024-1750)].values
        TOA_correction = (tot_volc_erf - np.mean(tot_volc_erf[90:100])+ contrails)
        TOA_correction[np.where(TOA_correction> -0.25)] = 0
        #clamp down some of the uncertainties so that the KF pays most attention to the temperatures
        ekf.Q = ekf.Q * 0.25
        ekf.cM2 = ekf.cM2 * 0.5
        ekf.dfrA_float_var = ekf.dfrA_float_var *.25
        #now this is somewhat flat and only makes corrections when there is a major volcanic eruption
        erf_trapez=(1*TOA_correction[0:-2]+2*TOA_correction[1:-1]+0*TOA_correction[2:] )/3 

        ohca_meas_all = np.load(config.CODEBASE_PATH +"/OHC_ensemble/ohca_sampled_ensemble.npy")
        ohca_meas = ohca_meas_all[model_run,0:ekf.n_iters]
        ohca_meas = ohca_meas - ohca_meas[0]
        ohca_uncert_all = np.load(config.CODEBASE_PATH +"/OHC_ensemble/ohca_uncertainty_ensemble.npy")
        ohca_uncert = np.abs(ohca_uncert_all[model_run,0:ekf.n_iters]) #remember stored as complex number, this gives us the RMS uncertainty
        ekf.Roc_tvar=np.square(ohca_uncert+5)/ekf.zJ_from_W/ekf.zJ_from_W #set minimum of 5ZJ uncertainty in 2010 standard year
        

        import os
        cur_path = os.path.dirname(os.path.realpath(__file__))
        #also reading this from a FaIR file
        TOA_meas_artif_all = np.load(cur_path + "/heads_tails_forcing/all-headstails_current_TOA_ensemble.npy")
        TOA_meas_artif = TOA_meas_artif_all[model_run,0,0:ekf.n_iters] #.to_numpy() if pandas array
        TOA_meas_artif1 =  TOA_meas_artif - TOA_correction *.75 #(tot_volc_erf ) * .4 #again removing the volcanic signal for the ta run
        
        data3 =  np.genfromtxt(open(config.CLIMATE_DATA_PATH+"/SSP_inputdata/KF6projectionSSP.csv", "rb"),dtype=float, delimiter=',')
        #ekf.anthro_clouds = np.append(ekf.anthro_clouds,[data3[2024-2015,6+2]+1])
        ekf.anthro_clouds = np.append(ekf.anthro_clouds,[2024*4.482E-5 -.018479 ])
        new_tsi = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp119_1750-2500.csv")['solar'][100:(100+ekf.n_iters+1)].values
        
        new2_tsi = np.append(ekf.tsi,[new_tsi[2024-1850]+ekf.sw_in-np.mean(new_tsi)])
        ekf.lCo2= np.append(ekf.lCo2,[np.log10(data3[2024-2015,1+2]-data3[2023-2015,1+2] + ekf.data[2023-1850,2])])
        ekf.tsi = 2* ta_smooth(new2_tsi,ekf.sw_in,optical=False,N=34) - ekf.sw_in  #np.full(np.shape(ekf.tsi),np.mean(ekf.tsi))
        #makes centered mean
        ekf.tsi[16:] = 2*ekf.tsi[16:]-ekf.tsi[:-16]

        temps[9:ekf.n_iters] = temps[9:ekf.n_iters] - erf_trapez[:-7]/(ekf.heatCp - ekf.Cs)/2 #surface temp of 20 yr mean still depressed
        ohca_meas[9:ekf.n_iters] = ohca_meas[9:ekf.n_iters] - (erf_trapez[:-7]/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*50/2
        new_observ = np.array([[temps[:ekf.n_iters],ohca_meas/ekf.zJ_from_W ,TOA_meas_artif1]]).T
        ekf.dfrA_float_var = ekf.dfrA_float_var/40
        #lots of trouble to add one more output point
        ndim=4
        sz = (ekf.n_iters,ndim,1) # size of array
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

        ekf.R_tvar =np.zeros(ekf.n_iters)
        ekf.Roc_tvar = np.zeros(ekf.n_iters)
        ekf.TOA_meas_artif_var = np.zeros(ekf.n_iters)
        means[0:ekf.n_iters], ses[0:ekf.n_iters],means2[0:ekf.n_iters],ses2[0:ekf.n_iters] = ekf.ekf_run(new_observ,ekf.n_iters,retPs=3)
        
        #import matplotlib.pyplot as plt
        #breakpoint()
        #plt.plot(years,means +given_preind_base - ekf.offset-preind_base - offset_satcal)
        #plt.show()
        #we're already putting back the given preind base
            
        return means +given_preind_base - ekf.offset-preind_base , np.sqrt(np.abs(ses))*1.2,means2 +given_preind_base -ekf.offset-preind_base ,np.sqrt(np.abs(ses2))

    else:
        restore_baselines()
        new_iter=len(years)
        given_preind_base = np.mean(temperature[0:50])
        erf_data = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")

        if (exp_attr[1]=='ESM1-2-LR'):
            ekf.n_iters = 174 
            unf_new_opt_depth= np.full(new_iter, (1/ekf.involcavg-9.7279))
            unf_new_opt_depth[0:ekf.n_iters]=ekf.data[:,3]*0.001 #infill old optical depths based on existing historical AOD record, later ta-smoothed
            
            ohca_later =Dataset(config.CLIMATE_DATA_PATH+"/ESM1-2-LR/opottempmint/"+exp_attr[2].lower()+"_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_l = ohca_later[:].__array__()
            ohca_ly = average_every_n(ohca_l[model_run,:], 12)
            ohca_earlier = Dataset(config.CLIMATE_DATA_PATH+"/ESM1-2-LR/opottempmint/historical_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_e = ohca_earlier[:].__array__()
            ohca_ey = average_every_n(ohca_e[model_run,:], 12)
            ohca_meas = np.concatenate((ohca_ey,ohca_ly)).astype(np.float64)
            ohca_meas = ohca_meas - ohca_meas[0] #start at 0 #supplies full record

            rt_later =Dataset(config.CLIMATE_DATA_PATH+"/ESM1-2-LR/rt/global_rt_"+exp_attr[2].lower()+".nc", 'r').variables['rt']
            rt_l = rt_later[:].__array__()
            rt_ly = average_every_n(rt_l[model_run,:], 12)
            rt_earlier = Dataset(config.CLIMATE_DATA_PATH+"/ESM1-2-LR/rt/global_rt_historical.nc", 'r').variables['rt']
            rt_e = rt_earlier[:].__array__()
            rt_ey = average_every_n(rt_e[model_run,:], 12)
            new_rtTOA = np.concatenate((rt_ey,rt_ly))

        
        elif (exp_attr[1]=='NorESM'):
            # need to change thermal conductivity ??
            ekf.n_iters = 174 
            unf_new_opt_depth= np.full(new_iter, (1/ekf.involcavg-9.7279))
            #unf_new_opt_depth[0:ekf.n_iters]=ekf.data[:,3]*0.001 #infill old optical depths based on existing historical AOD record, later ta-smoothed DONE - replace with NorESM simulation
            aod_spinup = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/NorESM1-M-historical/hist_aod.nc", 'r').variables['__xarray_dataarray_variable__']
            aod_s = aod_spinup[:].__array__()
            long_past_index = ((gen_orig_number(model_run,60) -1 )// 20)
            unf_new_opt_depth[(1850-1850):(2006-1850)]= average_every_n(aod_s[long_past_index,:], 12)/1000
            
            ohca_meas=np.zeros(new_iter).astype(np.float64)
            #ohca_meas[0:ekf.n_iters] = ekf.ocean_heat_measured/ekf.zJ_from_W   #DONE - replace with NorESM simulation
            ohca_spinup = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/NorESM1-M-historical/hist_ohca_mon.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_s = ohca_spinup[:].__array__()
            ohca_s_yr = average_every_n(ohca_s[long_past_index,:], 12)
            ohca_meas[(1850-1850):(2006-1850)]= ohca_s_yr - ohca_s_yr[0]
            
            #ohca_earlier = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/OHCA/historicalVolc_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_earlier = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/OHCA/historicalVolc_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
            ohca_e = ohca_earlier[:].__array__()
            #ohca_meas[(1980-1850):(2006-1850)] = average_every_n(ohca_e[model_run,:], 12) + ohca_meas[(1979-1850)]
            ohca_meas[(1980-1850):(2006-1850)] = ohca_e[(model_run if model_run< 37 else model_run-1) ,:] + 2*ohca_meas[(1979-1850)] - ohca_meas[(1978-1850)] #missing one model_run ohca, so shift them

            aod_earlier = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/BethkeEtAl2017/historicalVolc_aod.nc", 'r').variables['__xarray_dataarray_variable__']
            aod_e = aod_earlier[:].__array__()
            unf_new_opt_depth[(1980-1850):(2006-1850)]= average_every_n(aod_e[model_run,:], 12)/1000
            
            new_rtTOA = np.zeros(new_iter)
            rt_spinup = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/NorESM1-M-historical/rtmt/rtmt_NorESM1-M_historical.nc", 'r').variables['rtmt']
            rt_s = rt_spinup[:].__array__()
            new_rtTOA[(1850-1850):(2006-1850)]=average_every_n(rt_s[long_past_index,:], 12)
            rt_earlier = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/BethkeEtAl2017/historicalVolc_rtmt/rtmt_NorESM1-M_historicalVolc.nc", 'r').variables['rtmt']
            rt_e = rt_earlier[:].__array__()
            new_rtTOA[(1980-1850):(2006-1850)]=average_every_n(rt_e[model_run,:], 12)
            if exp_attr[3]=='Volc':
                #opt_depth        
                aod_later = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/BethkeEtAl2017/rcp45Volc_aod.nc", 'r').variables['__xarray_dataarray_variable__']
                aod_l = aod_later[:].__array__()
                unf_new_opt_depth[(2006-1850):(2100-1850)]= average_every_n(aod_l[model_run,:], 12)/1000
                #ohca
                ohca_later = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/OHCA/rcp45Volc_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
                ohca_l = ohca_later[:].__array__()
                ohca_meas[(2006-1850):(2100-1850)]= ohca_l[model_run,:] + 2*ohca_meas[(2005-1850)]-ohca_meas[(2004-1850)]
                rt_later =Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/BethkeEtAl2017/rcp45"+exp_attr[3]+"_rtmt/rtmt_NorESM1-M_rcp45"+exp_attr[3]+".nc", 'r').variables['rtmt']
                rt_l = rt_later[:].__array__()
                new_rtTOA[(2006-1850):(2100-1850)]= average_every_n(rt_l[model_run,:], 12)
            

            elif exp_attr[3]=='VolcConst':
                aod_later = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/BethkeEtAl2017/rcp45VolcConst_partial20_aod.nc", 'r').variables['__xarray_dataarray_variable__']
                aod_l = aod_later[:].__array__()
                unf_new_opt_depth[(2006-1850):(2100-1850)]= average_every_n(aod_l[model_run%20,:], 12)/1000
                ohca_later = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/OHCA/rcp45VolcConst_partial20_ohca.nc", 'r').variables['__xarray_dataarray_variable__']
                ohca_l = ohca_later[:].__array__()
                ohca_meas[(2006-1850):(2100-1850)]=ohca_l[model_run if (model_run>2 and model_run<14) else model_run%14,:]  + 2*ohca_meas[(2005-1850)]-ohca_meas[(2004-1850)]
                rt_later =Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/BethkeEtAl2017/rcp45"+exp_attr[3]+"_rtmt/rtmt_NorESM1-M_rcp45"+exp_attr[3]+"_20.nc", 'r').variables['rtmt']
                rt_l = rt_later[:].__array__()
                new_rtTOA[(2006-1850):(2100-1850)]= average_every_n(rt_l[model_run%20,:], 12)
                
            
            contrails = erf_data['contrails'][(1980-1750):(2100-1750)].values
            unf_new_opt_depth[(1980-1850):] = unf_new_opt_depth[(1980-1850):] + (contrails-0.015)/.18*.04  # get rid of gradual decline in the baseline over 21st century


                
        new_opt_depth = ta_smooth(unf_new_opt_depth,(1/(unf_new_opt_depth[20]+9.7279))) #ekf.involcavg)
        temps = temperature-given_preind_base+preind_base + ekf.offset #ensure this matches expected starting temperature in 1850

        frac_blocked = 1-(5.5518/(unf_new_opt_depth+9.9735)) #(9.068/(comb_recAOD_cleaned+9.7279))
        erf_data_solar = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")['solar'].to_numpy()[(1850-1750):]
        solar_full = erf_data_solar[:len(temps)]+ 340.4099428 - 0.108214
        tot_volc_erf = - solar_full*frac_blocked + 151.22
       # erf_trapez=(1*tot_volc_erf[0:-2]+2*tot_volc_erf[1:-1]+0*tot_volc_erf[2:] )/4 #want to snap back as quickly as possible
        #we are running the _ta method, so we want to feed in the temperature and ocean heat record that is already as close as possible to the 20yr mean

        contrails = erf_data['contrails'][(1850-1750):(len(temps)+100)].values
        TOA_correction = (tot_volc_erf - np.mean(tot_volc_erf[90:100])+ contrails)
        TOA_correction[np.where(TOA_correction> -0.25)] = 0 #-0.25 #correction to keep floating TOA adjustment from going wild
        #TOA_correction-= np.mean(TOA_correction) 
        #TOA_correction[0:75] = TOA_correction[0:75] + 0.2 don't need this fix
        erf_trapez=(1*TOA_correction[0:-2]+2*TOA_correction[1:-1]+0*TOA_correction[2:] )/3 #want to snap back as quickly as possible

        #erf_trapez-= np.mean(erf_trapez) #only going to surf temp and ohca
       # breakpoint()

        #temps[2:] = temps[2:] - erf_trapez/(ekf.heatCp - ekf.Cs)
        #ohca_meas[2:] = ohca_meas[2:] - (erf_trapez/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*50 #not perfect but good for now, also unclear where the 50 comes from
        temps[9:] = temps[9:] - erf_trapez[:-7]/(ekf.heatCp - ekf.Cs)/2 #surface temp of 20 yr mean still depressed
        ohca_meas[9:] = ohca_meas[9:] - (erf_trapez[:-7]/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*50/2
            
        if ekf.n_iters != new_iter:
            new_tsi = np.full(new_iter, ekf.sw_in)
            new_tsi = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp119_1750-2500.csv")['solar'][100:(100+new_iter+1)].values
            new_tsi = new_tsi + ekf.sw_in-np.mean(new_tsi)

            new_R_tvar =np.full(new_iter, np.mean(ekf.R_tvar[150:174]))
            new_R_tvar[0:ekf.n_iters]=ekf.R_tvar[0:ekf.n_iters]
            

            new_Roc_tvar =np.full(new_iter, np.mean(ekf.Roc_tvar[-50:])) #np.mean(ekf.Roc_tvar[-10:]))
            new_Roc_tvar[0:ekf.n_iters]= ekf.Roc_tvar[0:ekf.n_iters]
            
            data3 =  np.genfromtxt(open(config.CLIMATE_DATA_PATH+"/SSP_inputdata/KF6projectionSSP.csv", "rb"),dtype=float, delimiter=',')
            SSPnames=[126,434,245,370,585]
            if exp_attr[2]=='RCP45':
                find_case = 245
            else:
                find_case = int(exp_attr[2][3:])
            rcp = SSPnames.index(find_case)
            handoffyr = 1850+ekf.n_iters
            if (exp_attr[1]=='ESM1-2-LR'):
                new_Co2_df = pd.read_csv(open(config.CLIMATE_DATA_PATH+"/SSP_inputdata/eCO2_"+exp_attr[1]+"_"+exp_attr[2].lower()+".csv"),dtype=float, delimiter=',')
                new_lCo2 = np.log10(new_Co2_df['eCO2'].values)
                new_lCo2 -= new_lCo2[0]
                #warming_scale = ( np.average(ekf.temps[(1990-1850):(2010-1850)]) - np.average(ekf.temps[0:50])) / ( np.average(temperature[(1990-1850):(2010-1850)]) - np.average(temperature[0:50]))
                new_lCo2 *= (ekf.lCo2[(2013-1850)] - ekf.lCo2[0])/ new_lCo2[(2013-1850)]  #1 for 370, .9 for 245 and 126
                ekf.fdbkA = 0.25 #this going way down was helpful
                ekf.gad = 0.67 * 0.85
                new_lCo2 += ekf.lCo2[0]
                erf_data = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERF_"+exp_attr[1]+"_"+exp_attr[2].lower()+".csv")
                erf_natvolc = erf_data['ERF_natural']
            elif (exp_attr[1]=='NorESM'):
                from . import gen_eCO2
                if exp_attr[3]=='Volc':
                    erf_data = pd.read_csv(config.CODEBASE_PATH+"/Methods/43_Forcing_Based/3_Human_Induced/GWI_data/ERFanthro_NorESM_rcp45Volc_full.csv")
                    #erf_data = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERF_NorESM_rcp45Volc.csv")
                    erf_natvolc = pd.read_csv(config.CODEBASE_PATH+"/Methods/43_Forcing_Based/3_Human_Induced/GWI_data/ERFnatural_NorESM_rcp45Volc_full.csv").to_numpy()[:,model_run+1]
                elif exp_attr[3]=='VolcConst':
                    erf_data = pd.read_csv(config.CODEBASE_PATH+"/Methods/43_Forcing_Based/3_Human_Induced/GWI_data/ERF_NorESM_rcp45VolcConst_full.csv")
                    erf_natvolc = erf_data['ERF_natural']
                model_outputlCo2 = gen_eCO2.calculate_equivalent_co2(erf_data['ERF_wmghg'].values) #erf_data['ERF_anthro'].values)
                new_lCo2 = np.log10(model_outputlCo2)
                #rescaling because the Co2 scaling doesnt match somehow
                new_lCo2 -= new_lCo2[0]
                #warming_scale = ( np.average(ekf.temps[(1990-1850):(2010-1850)]) - np.average(ekf.temps[0:50])) / ( np.average(temperature[(1990-1850):(2010-1850)]) - np.average(temperature[0:50]))
                new_lCo2 *= (ekf.lCo2[(2013-1850)] - ekf.lCo2[0])/ new_lCo2[(2013-1850)] * 0.93 #allowed to adjust since it's wmghg *0.88 with standard gad #0.9 works well
                new_lCo2 += ekf.lCo2[0]

                ekf.gad = 0.67 * 1.35
                
                #np.concatenate((np.log10(ekf.data[:(1980-1850),2]), np.log10(ekf.data[(1981-1850),2] - model_outputlCo2[0]  + model_outputlCo2)))            #compute_update(xhatf[k-1],0,volcs[k-1],case[startyr-2015+k-1],caseA[startyr-2015+k-1]+1)
            #np.log10(data3[:,1+rcp]),data3[:,6+rcp]
            new_anthro_clouds = np.concatenate((ekf.anthro_clouds,data3[handoffyr-2015: 1850+new_iter-2015,6+rcp]+1))

            #clamp down some of the uncertainties so that the KF pays most attention to the temperatures

            ekf.Q = ekf.Q * 0.15
            ekf.cM2 = ekf.cM2 * 0.4
            ekf.dfrA_float_var = ekf.dfrA_float_var *.15
            
            ekf.opt_depth=new_opt_depth # * 2- 0.005
            ekf.tsi = new_tsi
            #same as in historical - shouldn't affect much, just smooths out this input
            ekf.tsi = 2* ta_smooth(ekf.tsi,ekf.sw_in,optical=False,N=34) - ekf.sw_in  #np.full(np.shape(ekf.tsi),np.mean(ekf.tsi))
            #makes centered mean
            ekf.tsi[16:] = 2*ekf.tsi[16:]-ekf.tsi[:-16]

            ekf.R_tvar = new_R_tvar  #/10 SR15 try
            ekf.Roc_tvar=new_Roc_tvar
            ekf.lCo2=new_lCo2
            ekf.anthro_clouds = new_anthro_clouds
            #TOA_means_artif = ta_smooth(new_rtTOA,0.00001,optical=False) + (new_rtTOA - erf_natvolc)/2 #   + np.mean(erf_volc) #half is the prior 30 yrs, half is what we expect the TOA balance to be removing this year's volcanism
            ekf.T02= np.mean(temps[(2000-1850):(2005-1850)]) #ensure critical T02 temperature is correct
            ekf.temps = temps
            ekf.precompute_coeffs(True) #calculates new CO2 parameters relative to observed temperature change

            
            #TOA_meas_artif0 =  new_rtTOA - (erf_natvolc - erf_data_solar[:len(new_rtTOA )]) * .75  #+ np.mean(erf_volc)  to be removing this year's volcanism
            TOA_meas_artif0 =  new_rtTOA.astype(np.float64) #TOA_meas_artif1*TOA_scale
            TOA_meas_artif = TOA_meas_artif0 - np.mean(TOA_meas_artif0[0:70])/2 + (np.mean(ekf.TOA_meas_artif[2001-1850:2023-1850]) -np.mean( TOA_meas_artif0[2001-1850:2023-1850]))/2
                      #this should be balanced at the start originally, and also match recent observations
            #TOA_scale = np.mean(ekf.TOA_crop)/np.mean(TOA_meas_artif1[2001-1850:2023-1850]) 


        

            #ekf.Q[0,0]=ekf.Q[0,0]/2;
            
            #Cd= 136.5*deepdepth/1000
            #np.sum(TOA_meas_artif)
            print("CONSISTENCY CHECK")
            
            surf_heat = (temps[-1]-np.mean(temps[0:50]))*(ekf.heatCp-ekf.Cs)
            print(f"OHCA + surf = {ohca_meas[-1]/ekf.zJ_from_W + surf_heat}")
            print(f"sum TOA = {np.sum(TOA_meas_artif)}")
            print("Error:")
            TOA_error = (np.sum(TOA_meas_artif)-ohca_meas[-1]/ekf.zJ_from_W - surf_heat)/np.sum(TOA_meas_artif)
            print(TOA_error)
            
            ohca_scale = (np.sum(TOA_meas_artif)- surf_heat)/ohca_meas[-1]*ekf.zJ_from_W
            print("OHCA_scale:")
            print(ohca_scale)

            #now fix the deep ocean depth / gad which means fitting values for blind model
            #preferably run just up to 2020s
            
            #ekf.precompute_coeffs(False)
            #ekf.n_iters=new_iter

        

        #new_observ = np.transpose(np.array([temperature-given_preind_base+preind_base + ekf.offset,ohca_meas/ekf.zJ_from_W ]))
        ekf.TOA_meas_artif_var = np.full(new_iter,ekf.TOA_crop_var/4) #just put constant variance on this
        ndim =4
        if True: #(abs(TOA_error < 0.2)):
            TOA_meas_artif1 =  TOA_meas_artif* (1-TOA_error)  - TOA_correction * .75 #TOA_meas_artif # already corrected
            #TOA_meas_artif1 =  TOA_meas_artif  - TOA_correction * .75
            #print("TOA correction")
           # print((TOA_correction))
            #ohca_meas = ohca_meas*ohca_scale
        ##else:
        #   TOA_meas_artif1 =  TOA_meas_artif  - TOA_correction * .75
        #    ohca_meas = ohca_meas*ohca_scale
        #    #ekf.Roc_tvar[170:]=new_Roc_tvar[170:]*4 #expanding uncertainty here since we had to correct the ocean
        #temps2=np.copy(temps) #repeating the logic of SR15 here
        #temps2[15:] = 15/8*ta_smooth(temps,0,optical=False, N=30)[15:]
        #temps2[7:] = 2*temps2[7:]-temps2[:-7]
        
        new_observ = np.array([[temps,ohca_meas/ekf.zJ_from_W,TOA_meas_artif1]]).T #*ohca_scale testing temperature offset
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

        
        return means - ekf.offset-preind_base + given_preind_base , np.sqrt(np.abs(ses))*2,means2 -ekf.offset-preind_base + given_preind_base,np.sqrt(np.abs(ses2))*2


   
    #(trailfilter_xhat,P2s,S2s,trailfilter_xhatm)=ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3)
