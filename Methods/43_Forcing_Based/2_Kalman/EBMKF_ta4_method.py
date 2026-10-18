from pathlib import Path as _V11Path
import sys as _V11Sys
_V11Sys.path.insert(0,str(_V11Path(__file__).resolve().parent))
# Allow this existing method file to be launched directly as a diagnostic tool.
if __name__ == '__main__' and not __package__:
    _V11Sys.path.insert(0,str(_V11Path(__file__).resolve().parents[3]))
    __package__ = 'Methods.43_Forcing_Based.2_Kalman'

from ebmkf_v11_support import (nonnegative_sqrt as _v11_sqrt,
    prepare_observed_two_state as _v11_prepare_observed,
    scenario_aerosols as _v11_scenario_aerosols, volcanic_erf as _v11_volcanic_erf,
    parameter_choice_density as _v11_choice_density,
    plot_four_state_components as _v11_plot_components,
    four_state_diagnostic_series as _v11_diagnostic_series,
    plot_four_state_correction as _v11_plot_correction)
import numpy as np
import pandas as pd
import config
#EBMKF_Nicklas as ekf
from netCDF4 import Dataset
import pdb
import xarray as xr
import os as _os
import pdb


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
        wt_opt_depths = 1/ekf.volcanic_denominator(orig_opt_depths)
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
        return ekf.raw_aod_from_inverse_denominator(nwt_opt_depths)
    else:
        return(nwt_opt_depths)


# Future starting parameters copied from the author-tuned ta2 cores on
# 2026-10-17: MPI uses Nicklas2b; NorESM uses Nicklas2c. These are explicit
# ta4 setup values for review, not imports that run another model or auto-fit.
# gad is the starting exchange coefficient; the existing 2020-OHC adjustment
# inside ekf_run remains active for long future records.
FUTURE_PARAMETER_SETTINGS = {
    'ESM1-2-LR': {
        'gad': 0.75, 'fdbkA': 0.20,
        'AEROSOL_FORCING_PERCENTILE': 85.0, 'eCO2_FORCING_PERCENTILE': 22.0,
    },
    'NorESM': {
        'gad': 1.05, 'fdbkA': 0.35,
        'AEROSOL_FORCING_PERCENTILE': 65.0, 'eCO2_FORCING_PERCENTILE': 45.0,
    },
}

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
    ekf._orig_TOA_meas_artif = ekf.TOA_meas_artif.copy() 
    ekf._orig_lCo2 = ekf.lCo2.copy()
    ekf._orig_anthro_clouds = ekf.anthro_clouds.copy()
    # Scalars
    ekf._orig_dfrA_float_var = ekf.dfrA_float_var
    ekf._orig_gad = ekf.gad
    ekf._orig_fdbkA = ekf.fdbkA
    ekf._orig_AEROSOL_FORCING_PERCENTILE = ekf.AEROSOL_FORCING_PERCENTILE
    ekf._orig_eCO2_FORCING_PERCENTILE = ekf.eCO2_FORCING_PERCENTILE
    ekf._orig_temps = ekf.temps.copy()
    ekf._orig_oc_meas = ekf.oc_meas.copy()

def restore_baselines():
    ekf.use_observational_forcing_percentiles = True
    ekf.Q = ekf._orig_Q.copy()
    ekf.cM2 = ekf._orig_cM2.copy()
    ekf.opt_depth = ekf._orig_opt_depth.copy()
    ekf.tsi = ekf._orig_tsi.copy()
    ekf.R_tvar = ekf._orig_R_tvar.copy()
    ekf.Roc_tvar = ekf._orig_Roc_tvar.copy()
    ekf.TOA_meas_artif = ekf._orig_TOA_meas_artif.copy()
    ekf.TOA_meas_artif_var = ekf._orig_TOA_meas_artif_var.copy()
    ekf.lCo2 = ekf._orig_lCo2.copy()
    ekf.anthro_clouds = ekf._orig_anthro_clouds.copy()
    ekf.dfrA_float_var = ekf._orig_dfrA_float_var
    ekf.gad = ekf._orig_gad
    ekf.fdbkA = ekf._orig_fdbkA
    ekf.AEROSOL_FORCING_PERCENTILE = ekf._orig_AEROSOL_FORCING_PERCENTILE
    ekf.eCO2_FORCING_PERCENTILE = ekf._orig_eCO2_FORCING_PERCENTILE
    ekf.temps = ekf._orig_temps.copy()
    ekf.oc_meas = ekf._orig_oc_meas.copy()
    ekf.precompute_coeffs(False)
    

def _apply_future_parameters(model_name):
    """Select model-specific starts before constructing future forcing arrays.

    Recompute the selected observational aerosol/eCO2 prefix first. The aerosol
    choice also selects the matching SSP curve; eCO2 keeps the existing
    model-derived future recipe and its observational normalization. Historical
    and heads/tails calls use restore_baselines() and never enter this helper.
    """
    for name, value in FUTURE_PARAMETER_SETTINGS[model_name].items():
        setattr(ekf, name, value)
    ekf.precompute_coeffs(False)


def run_method(years, temperature0, uncert, model_run, experiment_type, diagnostics=False):
    diagnostic_requested = diagnostics or _os.environ.get("EBMKF_DEBUG_PLOTS") == "1"
    temperature = temperature0.astype(np.float64)
    #temperature has 2024, dont have all forcings for this yet - must update ekf.n_iters
    data_orig = pd.read_csv("./Common_Data/HadCRUT.5.1.csv")
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
        _v11_choice_density(ekf)
        ekf.n_iters = len(years)
        unf_new_opt_depth = ekf.data[:ekf.n_iters,3]*0.001
        ekf.opt_depth=ta_smooth(unf_new_opt_depth,(1/ekf.volcanic_denominator(unf_new_opt_depth[20])),N=40)*1.4 #want a muted volcanic signal
        given_preind_base = np.mean(temperature[0:50])
        temps = temperature-given_preind_base+preind_base + ekf.offset #ensure this matches expected starting temperature in 1850

        erf_data= pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")
        erf_data_solar = erf_data['solar'].to_numpy()[(1850-1750):]
        solar_full = erf_data_solar[:ekf.n_iters]+ 340.4099428 - 0.108214
        tot_volc_erf = _v11_volcanic_erf(ekf, unf_new_opt_depth, solar_full)
        contrails = erf_data['contrails'][(1850-1750):(1850+ekf.n_iters-1750)].values
        TOA_correction = (tot_volc_erf - np.mean(tot_volc_erf[90:100])+ contrails)
        #TOA_correction = TOA_correction #- np.mean(TOA_correction)*2/3
        #TOA_correction[0:75] = TOA_correction[0:75] + 0.2
        
        TOA_correction[np.where(TOA_correction> -0.25)] = 0
        #breakpoint()
        #clamp down some of the uncertainties so that the KF pays most attention to the temperatures
        ekf.Q = ekf.Q * 0.2
        ekf.cM2 = ekf.cM2 * .6
        ekf.dfrA_float_var = ekf.dfrA_float_var *.35 #.5 .25
        #now this is somewhat flat and only makes corrections when there is a major volcanic eruption
        erf_trapez=(1*TOA_correction[0:-2]+2*TOA_correction[1:-1]+0*TOA_correction[2:] )/3 #want to snap back as quickly as possible
        ohca_meas = ekf.oc_meas*ekf.zJ_from_W
        tsi_orig = ekf.tsi.copy()
        ekf.tsi = 2* ta_smooth(ekf.tsi,ekf.sw_in,optical=False,N=34) - ekf.sw_in  #np.full(np.shape(ekf.tsi),np.mean(ekf.tsi))
        #makes centered mean
        ekf.tsi[:18] = np.mean(ekf.tsi[34:68])
        ekf.tsi[16:] = 2*ekf.tsi[16:]-ekf.tsi[:-16]
        #projecting this cycle forward into future
        
        #temps[2:-1] = temps[2:-1] - erf_trapez/(ekf.heatCp - ekf.Cs)
        #ohca_meas[2:] = ohca_meas[2:] - (erf_trapez/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*50 #not perfect but good for now, also unclear where the 50 comes from
        pulse = erf_trapez / (ekf.heatCp - ekf.Cs) / 6
        temps[2:] -= pulse          # Warming correction 1–2 years after eruption
        temps[9:] += pulse[:-7]     # Cooling correction 8–9 years after eruption #surface temp of 20 yr mean still depressed
        ohca_meas[9:ekf.n_iters] = ohca_meas[9:ekf.n_iters] - (erf_trapez[:-7]/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*25  #25
        
        TOA_meas_artif1 =  ekf.TOA_meas_artif[:ekf.n_iters] - TOA_correction  #0.75
        
##        import matplotlib.pyplot as plt
##        plt.figure()
##        plt.plot(years,TOA_meas_artif1)
##        plt.figure()
##        plt.plot(years,ekf.tsi)
##        plt.figure()
##        plt.plot(years,ekf.opt_depth)
##        plt.figure()
##        plt.plot(years,temps)
##        plt.show()
        new_observ = np.array([[temps[:ekf.n_iters],ohca_meas[:ekf.n_iters]/ekf.zJ_from_W ,TOA_meas_artif1]]).T
        #temps = temperature-given_preind_base+preind_base + ekf.offset
        #new_observ = np.array([[temps[:ekf.n_iters],ekf.oc_meas,ekf.TOA_meas_artif]]).T
        ekf.dfrA_float_var = ekf.dfrA_float_var/40
        if diagnostics:
            global _last_historical_diagnostics
            _last_historical_diagnostics = {
                'years': np.asarray(years).copy(),
                'adjusted_observations': np.asarray(new_observ[:,:,0]).copy(),
                'raw_ocean_ZJ': ekf._orig_oc_meas[:ekf.n_iters].copy()*ekf.zJ_from_W,
                'raw_TOA': ekf._orig_TOA_meas_artif[:ekf.n_iters].copy(),
                'TOA_correction': TOA_correction.copy(),
                'raw_AOD': unf_new_opt_depth.copy(),
                'raw_solar': tsi_orig.copy(),
            }
        means[0:ekf.n_iters], ses[0:ekf.n_iters],means2[0:ekf.n_iters],ses2[0:ekf.n_iters] = ekf.ekf_run(new_observ,ekf.n_iters,retPs=3,diagnostics=False)
        if _os.environ.get('EBMKF_DEBUG_PLOTS') == '1':
            figures, diagnostic_series = _show_run_diagnostics(years,new_observ,experiment_type,model_run,
                                                               temperature0-given_preind_base+preind_base)
            if _os.environ.get('EBMKF_DEBUG_BREAKPOINT') == '1':
                # Inspect ekf, years, new_observ and diagnostic_series; c continues.
                breakpoint()

        #means[0:ekf.n_iters] = means[0:ekf.n_iters] + (tsi_orig-ekf.tsi)/(ekf.heatCp) #adding TSI as a final correction doesnt work
        return means +given_preind_base - ekf.offset-preind_base, _v11_sqrt(ses)*1.2,means2 +given_preind_base -ekf.offset-preind_base,_v11_sqrt(ses2)


    elif (exp_attr[0]=="histens"):
        restore_baselines()
        ekf.n_iters = len(years)
        new_iter=len(years)
        #keeping optical depth (volcanoes) and solar forcing the same
        unf_new_opt_depth = ekf.data[:ekf.n_iters,3]*0.001
        ekf.opt_depth=ta_smooth(unf_new_opt_depth,(1/ekf.volcanic_denominator(unf_new_opt_depth[20])),N=40)*1.4 #want a muted volcanic signal
        given_preind_base = np.mean(temperature[0:50])

        temps = temperature-given_preind_base+preind_base + ekf.offset #ensure this matches expected starting temperature in 1850, same baseline as HadCRUT

        erf_data= pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")
        erf_data_solar = erf_data['solar'].to_numpy()[(1850-1750):]
        solar_full = erf_data_solar[:ekf.n_iters]+ 340.4099428 - 0.108214
        tot_volc_erf = _v11_volcanic_erf(ekf, unf_new_opt_depth, solar_full)
        contrails = erf_data['contrails'][(1850-1750):(new_iter+100)].values
        TOA_correction = (tot_volc_erf - np.mean(tot_volc_erf[90:100])+ contrails)
        TOA_correction[np.where(TOA_correction> -0.25)] = 0
        #clamp down some of the uncertainties so that the KF pays most attention to the temperatures
        ekf.Q = ekf.Q * 0.2
        ekf.cM2 = ekf.cM2 * 0.6
        ekf.dfrA_float_var = ekf.dfrA_float_var *.35
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
        if(len(years)==175):
            TOA_meas_artif_all = np.load(cur_path + "/heads_tails_forcing/all-headstails_current_TOA_ensemble_24.npy")
        elif(len(years)==176):
            TOA_meas_artif_all = np.load(cur_path + "/heads_tails_forcing/all-headstails_current_TOA_ensemble_25.npy")
        TOA_meas_artif0 = TOA_meas_artif_all[model_run,0,0:ekf.n_iters] #.to_numpy() if pandas array
        TOA_meas_artif = TOA_meas_artif0 - np.mean(TOA_meas_artif0[0:70])/2 + (np.mean(ekf.TOA_meas_artif[2001-1850:2023-1850]) -np.mean( TOA_meas_artif0[2001-1850:2023-1850]))/2
        # Balance the FaIR-derived TOA input before volcanic preprocessing.
        # Historical uses its assembled observational/CERES TOA input directly.
        surf_heat = (temps[-1]-np.mean(temps[0:50]))*(ekf.heatCp-ekf.Cs)
        TOA_error = (np.sum(TOA_meas_artif)-ohca_meas[-1]/ekf.zJ_from_W - surf_heat)/np.sum(TOA_meas_artif)
        #this should be balanced at the start originally, and also match recent observations and match the ocean observations
        TOA_meas_artif1 =  TOA_meas_artif* (1-TOA_error) - TOA_correction #*.75 #(tot_volc_erf ) * .4 #again removing the volcanic signal for the ta run
        
        data3 =  np.genfromtxt(open(config.CLIMATE_DATA_PATH+"/SSP_inputdata/KF6projectionSSP.csv", "rb"),dtype=float, delimiter=',')
        #ekf.anthro_clouds = np.append(ekf.anthro_clouds,[data3[2024-2015,6+2]+1])
        #Did we figure out why we're appending and if we neded to append more to anthro_clouds, tsi, lCo2
        # Frozen 9c already contains the complete historical aerosol series.
        new_tsi = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp119_1750-2500.csv")['solar'][100:(100+ekf.n_iters+1)].values
        if(len(years)==175 and len(ekf.lCo2)<175):
            new2_tsi = np.append(ekf.tsi,[new_tsi[2024-1850]+ekf.sw_in-np.mean(new_tsi)])
            ekf.lCo2= np.append(ekf.lCo2,[np.log10(data3[2024-2015,1+2]-data3[2023-2015,1+2] + ekf.data[2023-1850,2])])
        else:
             new2_tsi = ekf.tsi
        ekf.tsi = 2* ta_smooth(new2_tsi,ekf.sw_in,optical=False,N=34) - ekf.sw_in  #np.full(np.shape(ekf.tsi),np.mean(ekf.tsi))
        #makes centered mean
        ekf.tsi[:18] = np.mean(ekf.tsi[34:68])
        ekf.tsi[16:] = 2*ekf.tsi[16:]-ekf.tsi[:-16]
        pulse = erf_trapez / (ekf.heatCp - ekf.Cs) / 6
        temps[2:] -= pulse          # Warming correction 1–2 years after eruption
        temps[9:] += pulse[:-7]     # Cooling correction 8–9 years after eruption #surface temp of 20 yr mean still depressed
        ohca_meas[9:ekf.n_iters] = ohca_meas[9:ekf.n_iters] - (erf_trapez[:-7]/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*25
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
        if(len(years)==175 and len(ekf.R_tvar)<175):
            ekf.R_tvar =np.append(ekf.R_tvar, (len(years)-174)*[np.average(ekf.R_tvar[-10:])])
            ekf.TOA_meas_artif_var = np.append(ekf.TOA_meas_artif_var, (len(years)-174)*[np.average(ekf.TOA_meas_artif_var[-10:])])
        means[0:ekf.n_iters], ses[0:ekf.n_iters],means2[0:ekf.n_iters],ses2[0:ekf.n_iters] = ekf.ekf_run(new_observ,ekf.n_iters,retPs=3,diagnostics=False)
        #import matplotlib.pyplot as plt
        #breakpoint()
        #plt.plot(years,means +given_preind_base - ekf.offset-preind_base - offset_satcal)
        #plt.show()
        #we're already putting back the given preind base
            
        return means +given_preind_base - ekf.offset-preind_base , _v11_sqrt(ses)*1.2,means2 +given_preind_base -ekf.offset-preind_base ,_v11_sqrt(ses2)

    else:
        restore_baselines()
        _apply_future_parameters(exp_attr[1])
        new_iter=len(years)
        given_preind_base = np.mean(temperature[0:50])
        erf_data = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")

        if (exp_attr[1]=='ESM1-2-LR'):
            ekf.n_iters = 174 
            unf_new_opt_depth= np.full(new_iter, ekf.raw_aod_from_inverse_denominator(ekf.involcavg))
            unf_new_opt_depth[0:ekf.n_iters]=ekf.data[:ekf.n_iters,3]*0.001 #infill old optical depths based on existing historical AOD record, later ta-smoothed
            
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
            unf_new_opt_depth= np.full(new_iter, ekf.raw_aod_from_inverse_denominator(ekf.involcavg))
            #unf_new_opt_depth[0:ekf.n_iters]=ekf.data[:ekf.n_iters,3]*0.001 #infill old optical depths based on existing historical AOD record, later ta-smoothed DONE - replace with NorESM simulation
            aod_spinup = Dataset(config.CLIMATE_DATA_PATH+"/NorESM_volc/NorESM1-M-historical/hist_aod.nc", 'r').variables['__xarray_dataarray_variable__']
            aod_s = aod_spinup[:].__array__()
            long_past_index = ((gen_orig_number(model_run,60) -1 )// 20)
            unf_new_opt_depth[(1850-1850):(2006-1850)]= average_every_n(aod_s[long_past_index,:], 12)/1000
            
            ohca_meas=np.zeros(new_iter).astype(np.float64)
            #ohca_meas[0:ekf.n_iters] = ekf.ocean_heat_measured[:ekf.n_iters]/ekf.zJ_from_W   #DONE - replace with NorESM simulation
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


                
        new_opt_depth = ta_smooth(unf_new_opt_depth, (1/ekf.volcanic_denominator(unf_new_opt_depth[20])), N=40)*1.4  # same volcanic smoothing as historical/histens
        temps = temperature-given_preind_base+preind_base + ekf.offset #ensure this matches expected starting temperature in 1850

        erf_data_solar = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")['solar'].to_numpy()[(1850-1750):]
        solar_full = erf_data_solar[:len(temps)]+ 340.4099428 - 0.108214
        tot_volc_erf = _v11_volcanic_erf(ekf, unf_new_opt_depth, solar_full)
       # erf_trapez=(1*tot_volc_erf[0:-2]+2*tot_volc_erf[1:-1]+0*tot_volc_erf[2:] )/4 #want to snap back as quickly as possible
        #we are running the _ta method, so we want to feed in the temperature and ocean heat record that is already as close as possible to the 20yr mean

        contrails = erf_data['contrails'][(1850-1750):(len(temps)+100)].values
        TOA_correction = (tot_volc_erf - np.mean(tot_volc_erf[90:100])+ contrails)
        TOA_correction[np.where(TOA_correction> -0.25)] = 0 #-0.25 #correction to keep floating TOA adjustment from going wild
        #TOA_correction-= np.mean(TOA_correction) 
        #TOA_correction[0:75] = TOA_correction[0:75] + 0.2 don't need this fix
        erf_trapez=(1*TOA_correction[0:-2]+2*TOA_correction[1:-1]+0*TOA_correction[2:] )/3 #want to snap back as quickly as possible

        #erf_trapez-= np.mean(erf_trapez) #only going to surf temp and ohca
        #breakpoint()

        pulse = erf_trapez / (ekf.heatCp - ekf.Cs) / 6
        temps[2:] -= pulse          # Warming correction 1–2 years after eruption
        temps[9:] += pulse[:-7]     # Cooling correction 8–9 years after eruption #surface temp of 20 yr mean still depressed
        ohca_meas[9:] = ohca_meas[9:] - (erf_trapez[:-7]/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*50/2
            
        if ekf.n_iters != new_iter:
            new_tsi = np.full(new_iter, ekf.sw_in)
            new_tsi = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp119_1750-2500.csv")['solar'][100:(100+new_iter)].values
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
                # Feedback/exchange starts are selected above with both percentiles.
                
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

                
                #np.concatenate((np.log10(ekf.data[:(1980-1850),2]), np.log10(ekf.data[(1981-1850),2] - model_outputlCo2[0]  + model_outputlCo2)))            #compute_update(xhatf[k-1],0,volcs[k-1],case[startyr-2015+k-1],caseA[startyr-2015+k-1]+1)
            #np.log10(data3[:,1+rcp]),data3[:,6+rcp]
            new_anthro_clouds = np.concatenate((ekf.data[:ekf.n_iters,7]-ekf.aerosol_reference,_v11_scenario_aerosols(ekf, years[ekf.n_iters:], experiment_type)))

            #clamp down some of the uncertainties so that the KF pays most attention to the temperatures

            ekf.Q = ekf.Q * 0.15
            ekf.cM2 = ekf.cM2 * 0.4
            ekf.dfrA_float_var = ekf.dfrA_float_var *.15
            
            ekf.opt_depth=new_opt_depth # * 2- 0.005
            ekf.tsi = new_tsi
            #same as in historical - shouldn't affect much, just smooths out this input
            ekf.tsi = 2* ta_smooth(ekf.tsi,ekf.sw_in,optical=False,N=34) - ekf.sw_in  #np.full(np.shape(ekf.tsi),np.mean(ekf.tsi))
            #makes centered mean
            ekf.tsi[16:] = 2*ekf.tsi[16:]-ekf.tsi[:-16] #perhaps try with only 8 years, as that brings right to present - TODO LATER

            ekf.R_tvar = new_R_tvar  #/10 SR15 try
            ekf.Roc_tvar=new_Roc_tvar
            ekf.lCo2=new_lCo2
            ekf.anthro_clouds = new_anthro_clouds
            #TOA_means_artif = ta_smooth(new_rtTOA,0.00001,optical=False) + (new_rtTOA - erf_natvolc)/2 #   + np.mean(erf_volc) #half is the prior 30 yrs, half is what we expect the TOA balance to be removing this year's volcanism
            ekf.T02= np.mean(temps[(2000-1850):(2005-1850)]) #ensure critical T02 temperature is correct
            ekf.temps = temps
            # Keep full model/scenario forcing arrays during future recalibration.
            ekf.use_observational_forcing_percentiles = False
            ekf.precompute_coeffs(True) #calculates new CO2 parameters relative to observed temperature change

            
            #TOA_meas_artif0 =  new_rtTOA - (erf_natvolc - erf_data_solar[:len(new_rtTOA )]) * .75  #+ np.mean(erf_volc)  to be removing this year's volcanism
            TOA_meas_artif0 =  new_rtTOA.astype(np.float64) #TOA_meas_artif1*TOA_scale
            print(np.mean(TOA_meas_artif0[0:70]))
            ceres_offset = (np.mean(ekf.TOA_meas_artif[2001-1850:2023-1850]) -np.mean( TOA_meas_artif0[2001-1850:2023-1850])) 
            print(ceres_offset)

            
            TOA_meas_artif = TOA_meas_artif0 - np.mean(TOA_meas_artif0[0:70]) #+ (np.mean(ekf.TOA_meas_artif[2001-1850:2023-1850]) -np.mean( TOA_meas_artif0[2001-1850:2023-1850]))/2
                      #this should be balanced at the start originally, and also match recent observations
            #TOA_scale = np.mean(ekf.TOA_crop)/np.mean(TOA_meas_artif1[2001-1850:2023-1850]) 


        

            #ekf.Q[0,0]=ekf.Q[0,0]/2;
            
            #Cd= 136.5*deepdepth/1000
            #np.sum(TOA_meas_artif)
            print("CONSISTENCY CHECK")
            
            surf_heat = (temps[-1]-np.mean(temps[0:50]))*(ekf.heatCp-ekf.Cs)
            print(f"OHCA + surf = {ohca_meas[-1]/ekf.zJ_from_W + surf_heat}")
            print(f"sum TOA = {np.sum(TOA_meas_artif)}")
            ohcafull = ohca_meas[-1]
            difference = np.sum(TOA_meas_artif) * ekf.zJ_from_W - ohcafull

            extras_dir = _os.path.join(_os.path.dirname(__file__), "extras")
            _os.makedirs(extras_dir, exist_ok=True)
            out_file = _os.path.join(extras_dir, f"{exp_attr[2].lower()}.csv")
            if not diagnostic_requested:
                with open(out_file, "a") as f:
                    f.write(f"{model_run},{ohcafull},{difference}\n")

            out_file2 = _os.path.join(extras_dir, f"RT{exp_attr[2].lower()}.csv")
            print("Baseline TOA characteristics")
            print(np.mean(ekf.TOA_meas_artif[2001-1850:2023-1850]))
            print(np.std(ekf.TOA_meas_artif[2001-1850:2023-1850]))
            if not diagnostic_requested:
                with open(out_file2, "a") as f:
                    f.write(f"{model_run},{np.mean(TOA_meas_artif0[0:70])},{ceres_offset},{np.std(TOA_meas_artif0[2001-1850:2023-1850])}\n")

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
            TOA_meas_artif1 = TOA_meas_artif * (1-TOA_error) - TOA_correction  # full volcanic correction, as in historical/histens
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

        ekf.n_iters = new_iter  # diagnostics and state arrays must use the same simulation window
        means , ses ,means2 ,ses2 = ekf.ekf_run(new_observ,new_iter,retPs=3,diagnostics=False)
        if _os.environ.get('EBMKF_DEBUG_PLOTS') == '1':
            figures, diagnostic_series = _show_run_diagnostics(years,new_observ,experiment_type,model_run,
                                                               temperature0-given_preind_base+preind_base)
            if _os.environ.get('EBMKF_DEBUG_BREAKPOINT') == '1':
                # Inspect ekf, years, new_observ and diagnostic_series; c continues.
                breakpoint()

        if not diagnostic_requested:
            import matplotlib.pyplot as plt
            fig , axVolc = plt.subplots(figsize=(8,3.5))  #5)) #
            plt.subplots_adjust(bottom=0.15) ############
            
            plt.plot(years,unf_new_opt_depth, color='maroon')
            axVolc.set_ylabel('AOD (unitless)') #'Fraction Converted to
            axVolc.set_xlabel('Year')
            axVolc.set_title( "Corrected, Pre-Filtered Volcanic Signal (member #10)") #"SSP370")             #f"Hidden SW Radiative Flux Adjustment for EBMKF-{modelID}")
            extras_dir = _os.path.join(_os.path.dirname(__file__), "extras")
            plt.savefig(extras_dir+f"/EBMKFta4-NorESM-AODcorr1.png",dpi=300,format="png")
        
            
           #plt.show()
            
        return means - ekf.offset-preind_base + given_preind_base , _v11_sqrt(ses)*2,means2 -ekf.offset-preind_base + given_preind_base,_v11_sqrt(ses2)*2


   
    #(trailfilter_xhat,P2s,S2s,trailfilter_xhatm)=ekf.ekf_run(ekf.observ,ekf.n_iters,retPs=3,diagnostics=False)


def plot_historical_diagnostics(save_dir=None, show=True):
    """Run the actual historical ta4 wrapper and show its inputs/states/budget.

    This is an opt-in diagnostic run, not a calibration or cache refresh.
    All plotted observations/forcings are captured from run_method, and the
    blind replay uses those same adjusted forcings without assimilation.
    """
    import json
    from pathlib import Path
    import matplotlib.pyplot as plt
    root = Path(__file__).resolve().parents[3]
    frame = pd.read_csv(root/'Common_Data/HadCRUT.5.1.csv')
    years = frame['Time'].to_numpy()
    temperature = frame['Anomaly'].to_numpy()
    uncertainty = (frame['Lower'].to_numpy(), frame['Upper'].to_numpy())
    plt.close('all')  # discard incidental import-time figures
    result = run_method(years, temperature, uncertainty, 0, 'historical', diagnostics=True)
    d = _last_historical_diagnostics
    obs = d['adjusted_observations']
    n = len(years)
    filtered = ekf.xhat[:n,:,0].copy()
    retro = ekf.xhathat[:n,:,0].copy()
    # Replay precisely the wrapper's current forcing arrays, not import-time xblind.
    blind = np.zeros((n,3))
    blind[0] = [ekf.Teq1850, ekf.oc_meas[0], 0]
    sw = np.zeros(n); lw = np.zeros(n); exchange = np.zeros(n)
    factors = np.zeros((n,4))
    for k in range(n):
        state = blind[max(0,k-1)]
        if k:
            step = ekf.compute_update_hidden(np.array([[state[0],state[1],0.,1.]]).T,k)
            blind[k] = step[:3,0]
        T,H = state[:2]
        factors[k] = [ekf.tsi[k]/ekf.sw_in,
                      ekf.volcanic_denominator(ekf.volcanic_aod_reference)/ekf.volcanic_denominator(ekf.opt_depth[k]),
                      1+ekf.dfaA*(T-ekf.T02)+ekf.anthro_clouds[k]/ekf.reference_absorbed_sw,
                      1+ekf.dfaS*(T-ekf.T02)]
        sw[k] = ekf.reference_absorbed_sw*np.prod(factors[k])
        lw[k] = ekf.outbndf*T**(4-ekf.powp1)*(1-ekf.B0*ekf.lCo2[k])
        theta = (H-ekf.Cs*(T-ekf.Teq1850))/ekf.Cd+ekf.oc1850
        exchange[k] = ekf.gad*(T-ekf.Teq1850-theta+ekf.oc1850)
    storage = blind[:,1]+(ekf.heatCp-ekf.Cs)*(blind[:,0]-ekf.Teq1850)
    closure = np.diff(storage)-(sw-lw)[1:]
    assert np.max(np.abs(closure))<1e-9, 'Blind energy inventory fails to close'
    assert np.allclose(blind[1:,2],(sw-lw)[1:],atol=1e-10,rtol=0)
    # Show temperature anomalies on the source HadCRUT baseline, with absolute
    # reference temperatures printed separately. No mixed Celsius/Kelvin axes.
    anom = lambda t: np.asarray(t)-ekf.offset
    figs = []
    fig, axes = plt.subplots(4,2,figsize=(14,13),layout='constrained',sharex=True)
    figs.append(('ta4-inputs',fig)); a=axes.ravel()
    a[0].plot(years,temperature,color='0.65',label='HadCRUT5.1 annual')
    a[0].plot(years,anom(obs[:,0]),color='C3',label='Adjusted ta4 temperature input')
    a[0].set_ylabel('Temperature anomaly (°C)')
    a[1].plot(years,d['raw_ocean_ZJ'],color='0.65',label='Original OHC ensemble mean')
    a[1].plot(years,obs[:,1]*ekf.zJ_from_W,color='C3',label='Adjusted ta4 OHC input');a[1].set_ylabel('OHC (ZJ)')
    a[2].plot(years,d['raw_TOA'],color='0.65',label='Original synthetic/CERES TOA input')
    a[2].plot(years,obs[:,2],color='C3',label='Adjusted ta4 TOA input');a[2].plot(years,d['TOA_correction'],color='C2',label='Removed TOA correction');a[2].set_ylabel('Net TOA (W/m²)')
    a[3].plot(years,d['raw_AOD'],color='0.65',label='Frozen annual AOD')
    a[3].plot(years,ekf.opt_depth[:n],color='C3',label='ta4 smoothed/scaled AOD');a[3].set_ylabel('AOD (dimensionless)')
    a[4].plot(years,d['raw_solar'],color='0.65',label='Frozen annual TSI / 4')
    a[4].plot(years,ekf.tsi[:n],color='C3',label='ta4 solar input');a[4].set_ylabel('Solar input (W/m²)')
    a[5].plot(years,10**ekf.lCo2[:n],label='Equivalent CO₂ input');a[5].set_ylabel('Equivalent CO₂ (ppm)')
    a[6].plot(years,ekf.data[:n,7],color='0.65',label='Raw ARI + ACI')
    a[6].plot(years,ekf.anthro_clouds[:n],color='C3',label='Reference-centred ARI + ACI');a[6].set_ylabel('Aerosol ERF (W/m²)')
    a[7].axis('off')
    a[7].text(.02,.98,
        f"Teq1850: {ekf.Teq1850:.6f} K = {ekf.Teq1850-273.15:.6f} °C\n"
        f"Equilibrium reference years: {ekf.equilibrium_reference_years}\n"
        f"T02: {ekf.T02:.6f} K; reference years: {ekf.radiation_reference_years}\n"
        f"Absolute alignment: 1960–1989 mean = {ekf.JonesOffset:.2f} °C\n"
        f"Cs={ekf.Cs:.6f}; Cd={ekf.Cd:.6f}; heatCp={ekf.heatCp:.6f}\n"
        f"gamma={ekf.gad:g}; fdbkA={ekf.fdbkA:g}; fdbkS={ekf.fdbkS:g}; WV/LR={ekf.fdbkW:g}\n"
        f"beta2/dfaA={ekf.dfaA:.8f}; beta3/dfaS={ekf.dfaS:.8f}\n"
        f"Sref={ekf.reference_absorbed_sw:.6f} W/m²\n"
        f"Aerosol percentile={ekf.AEROSOL_FORCING_PERCENTILE:g}; eCO₂ percentile={ekf.eCO2_FORCING_PERCENTILE:g}\n"
        f"Volcanic sensitivity={ekf.volcanic_target_efficiency:g} ± {ekf.volcanic_target_efficiency_uncertainty:g} W/m²/AOD\n"
        "Ta4 multipliers: Q ×0.2; cM2 ×0.6; floating variance ×0.35/40\n"
        "These are input adjustments; they are not changes to the plotted baseline.",
        va='top',transform=a[7].transAxes,fontsize=10,linespacing=1.6)
    for ax in a[:7]:ax.legend(fontsize=8);ax.grid(alpha=.2)
    for ax in axes[-1]:ax.set_xlabel('Year')
    fig.suptitle('EBMKF-ta4: annual inputs and the actual adjusted inputs used by the filter')

    fig, axes = plt.subplots(3,2,figsize=(14,11),layout='constrained',sharex=True)
    figs.append(('ta4-states',fig));a=axes.ravel()
    centered=np.full(n,np.nan)
    for i in range(10,n-10):centered[i]=.5*(temperature[i-10:i+10].mean()+temperature[i-9:i+11].mean())
    a[0].plot(years,temperature,color='0.8',label='Annual temperature')
    a[0].plot(years,centered,'k',label='Centered 20-year benchmark')
    a[0].plot(years,anom(filtered[:,0]),label='ta4 filtered posterior')
    a[0].plot(years,anom(retro[:,0]),label='ta4 retrospective')
    a[0].plot(years,anom(blind[:,0]),'--',color='darkgoldenrod',label='Blind, same adjusted forcings')
    a[0].axhline(ekf.Teq1850-ekf.offset,color='0.4',ls=':',label='Initialization reference');a[0].set_ylabel('Temperature anomaly (°C)')
    a[1].plot(years,obs[:,1]*ekf.zJ_from_W,color='0.6',label='Adjusted OHC input')
    a[1].plot(years,filtered[:,1]*ekf.zJ_from_W,label='Filtered OHC')
    a[1].plot(years,blind[:,1]*ekf.zJ_from_W,'--',color='darkgoldenrod',label='Blind OHC');a[1].set_ylabel('OHC (ZJ)')
    for label,x,color in [('Filtered',filtered,'C0'),('Blind',blind,'darkgoldenrod')]:
        deep=(x[:,1]-ekf.Cs*(x[:,0]-ekf.Teq1850))/ekf.Cd+ekf.oc1850
        a[2].plot(years,deep-273.15,color=color,label=label+' deep temperature')
    a[2].set_ylabel('Deep temperature (°C)')
    a[3].plot(years,obs[:,2],color='0.6',label='Adjusted TOA input')
    a[3].plot(years,filtered[:,2],label='Posterior diagnosed TOA state')
    a[3].plot(years[1:],(sw-lw)[1:],'--',color='darkgoldenrod',label='Blind SW − LW');a[3].set_ylabel('Net TOA (W/m²)')
    a[4].plot(years,filtered[:,3],label='Inferred SW multiplier')
    a[4].axhline(1,color='0.5',ls=':');a[4].set_ylabel('SW multiplier (dimensionless)')
    a[5].plot(years,anom(obs[:,0])-anom(filtered[:,0]),label='Temperature input minus posterior')
    a[5].plot(years[1:],ekf.y[1:n,0,0],label='Temperature innovation (before update)');a[5].set_ylabel('Temperature residual (°C)')
    for ax in a:ax.legend(fontsize=8);ax.grid(alpha=.2)
    for ax in axes[-1]:ax.set_xlabel('Year')
    fig.suptitle('EBMKF-ta4: states, assimilation and blind replay on a common temperature baseline')

    fig, axes=plt.subplots(2,2,figsize=(13,8),layout='constrained',sharex=True)
    figs.append(('ta4-energy',fig));a=axes.ravel()
    a[0].plot(years,sw,label='Absorbed shortwave');a[0].plot(years,lw,label='Outgoing longwave');a[0].set_ylabel('Radiation (W/m²)')
    for j,label in enumerate(['Solar / S0','Volcanic transmission / reference','Cloud feedback + aerosol factor','Surface feedback factor']):a[1].plot(years,factors[:,j],label=label)
    a[1].set_ylabel('Multiplicative SW factors');a[1].set_title('Factors multiply; they are not additive forcing components')
    a[2].plot(years[1:],(sw-lw)[1:],label='SW − LW')
    a[2].plot(years[1:],exchange[1:],label='Surface-to-deep exchange')
    a[2].plot(years[1:],np.diff(storage),ls='--',label='Change in total storage');a[2].set_ylabel('W/m² (annual step)')
    a[3].plot(years[1:],closure,label='Δ[H + (C − Cs)ΔT] − (SW − LW)')
    a[3].set_ylabel('Blind energy residual (W/m²)');a[3].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    for ax in a:ax.legend(fontsize=8);ax.grid(alpha=.2)
    for ax in axes[-1]:ax.set_xlabel('Year')
    fig.suptitle('Blind physical budget: exchange cancels from total storage; assimilation is separate')
    correction_series = _v11_diagnostic_series(ekf,years,obs)
    figs.append(('ta4-correction',_v11_plot_correction(correction_series,'Historical EBMKF-ta4')))
    report = {
        'temperature_input':str(root/'Common_Data/toyKFmodelData9c.csv'),
        'observed_input':str(root/'Common_Data/HadCRUT.5.1.csv'),
        'helper_file':str(Path(__file__).with_name('ebmkf_v11_support.py')),
        'capacity_reference_file':str(root/'Common_Data/HadCRUT5_capacity_reference_5.0.2.0.csv'),
        'max_frozen_vs_observation_temperature_difference':float(np.max(np.abs(ekf.data[:n,1]-temperature))),
        'equilibrium_reference_years':list(ekf.equilibrium_reference_years),
        'Teq1850_K':float(ekf.Teq1850),'Teq1850_C':float(ekf.Teq1850-273.15),
        'T02_K':float(ekf.T02),'offset_K':float(ekf.offset),
        'AEROSOL_FORCING_PERCENTILE':float(ekf.AEROSOL_FORCING_PERCENTILE),
        'eCO2_FORCING_PERCENTILE':float(ekf.eCO2_FORCING_PERCENTILE),
        'volcanic_target_efficiency':float(ekf.volcanic_target_efficiency),
        'parameter_choice_density':_v11_choice_density(ekf,print_report=False),
        'fdbkA':float(ekf.fdbkA),'heatCp':float(ekf.heatCp),
        'max_blind_energy_residual':float(np.max(np.abs(closure))),
        'final_filtered_anomaly':float(anom(filtered[-1,0])),
        'final_adjusted_blind_anomaly':float(anom(blind[-1,0])),
        'final_original_blind_anomaly':float(anom(ekf.xblind[n-1,0,0])),
        'initial_actual_forcing_net_flux':float((sw-lw)[0]),
        'note':'Teq1850 is the calibration/initialization reference, not a solved equilibrium under actual 1850 forcing. Older capacity-reference HadCRUT does not set this temperature.'}
    print(json.dumps(report,indent=2))
    if save_dir:
        destination=Path(save_dir);destination.mkdir(parents=True,exist_ok=True)
        for name,fig in figs:fig.savefig(destination/(name+'.png'),dpi=160)
        (destination/'ta4-diagnostics.json').write_text(json.dumps(report,indent=2))
        pd.DataFrame({'year':years,'annual_temperature_anomaly':temperature,
                      'adjusted_temperature_anomaly':anom(obs[:,0]),
                      'filtered_temperature_anomaly':anom(filtered[:,0]),
                      'retrospective_temperature_anomaly':anom(retro[:,0]),
                      'blind_adjusted_temperature_anomaly':anom(blind[:,0]),
                      'blind_absorbed_SW':sw,'blind_outgoing_LW':lw,
                      'blind_ocean_exchange':exchange,'blind_total_storage_internal':storage}).to_csv(destination/'ta4-diagnostics.csv',index=False)
    if show:plt.show()
    return {'figures':figs,'report':report,'outputs':result}



def _show_run_diagnostics(years,new_observ,experiment_type,model_run,reference_temperature=None):
    import json
    import matplotlib.pyplot as plt
    from pathlib import Path
    # Import-time raw figures are separate from these active-wrapper diagnostics.
    plt.close('all')
    title=f'{experiment_type}, member {model_run} (EBMKF-ta4)'
    figures,series=_v11_plot_components(ekf,years,new_observ,title,reference_temperature=reference_temperature)
    destination=_os.environ.get('EBMKF_DEBUG_SAVE_DIR')
    if destination:
        destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
        prefix=f'{experiment_type}_{model_run}_ta4'
        for name,figure in zip(['components','innovations','TOA','correction'],figures):
            figure.savefig(destination/f'{prefix}_{name}.png',dpi=160)
        columns={'year':series['year'],'original_temperature_anomaly':series['reference_temperature_anomaly']}
        for label in ['filtered','retrospective']:
            for key in ['base_SW','SW_correction','absorbed_SW','outgoing_LW','net_flux']:
                columns[label+'_'+key]=series[label+'_'+key]
            columns[label+'_SW_multiplier']=series[label][:,3]
            columns[label+'_temperature_anomaly']=series[label][:,0]-ekf.offset
            columns[label+'_OHC_ZJ']=series[label][:,1]*ekf.zJ_from_W
        columns.update({'adjusted_temperature_anomaly':series['observations'][:,0]-ekf.offset,
                        'adjusted_OHC_ZJ':series['observations'][:,1]*ekf.zJ_from_W,
                        'adjusted_TOA':series['observations'][:,2],
                        'TOA_nominal_sigma':series['TOA_nominal_sigma'],
                        'TOA_effective_sigma':series['TOA_effective_sigma'],
                        'TOA_effective_lower95':series['observations'][:,2]-1.96*series['TOA_effective_sigma'],
                        'TOA_effective_upper95':series['observations'][:,2]+1.96*series['TOA_effective_sigma'],
                        'active_AOD':ekf.opt_depth[:len(years)],'active_solar':ekf.tsi[:len(years)],
                        'active_eCO2':10**ekf.lCo2[:len(years)],'active_aerosol_ERF':ekf.anthro_clouds[:len(years)]})
        pd.DataFrame(columns).to_csv(destination/f'{prefix}_diagnostics.csv',index=False)
        report={'experiment':experiment_type,'member':int(model_run),'year_start':int(years[0]),'year_end':int(years[-1]),
                'AEROSOL_FORCING_PERCENTILE':float(ekf.AEROSOL_FORCING_PERCENTILE),
                'eCO2_FORCING_PERCENTILE':float(ekf.eCO2_FORCING_PERCENTILE),'fdbkA':float(ekf.fdbkA),
                'gad_after_OHC_rescaling':float(ekf.gad),
                'future_starting_parameters':FUTURE_PARAMETER_SETTINGS.get(experiment_type.split('_')[1], {}) if experiment_type.startswith('fut_') else {},
                'volcanic_target_efficiency':float(ekf.volcanic_target_efficiency),
                'volcanic_target_efficiency_uncertainty':float(ekf.volcanic_target_efficiency_uncertainty),
                'correction_definition':series['correction_definition'],
                'first_correction_value': 'undefined: no prediction step at initialization',
                'correction_units':'W m-2; positive adds absorbed shortwave',
                'correction_note':'Same previous state and current-year forcings; not a change in total temperature and not the volcanic input correction.',
                'blind_note':'Replay at active post-wrapper parameters, with SW multiplier fixed to one.',
                'TOA_interval_definition':series['TOA_interval_definition'],
                'TOA_nominal_interval_definition':'Input +/- 1.96*sqrt(9*TOA_var), before added N*cM2 observation noise.'}
        (destination/f'{prefix}_diagnostics.json').write_text(json.dumps(report,indent=2))
    if _os.environ.get('EBMKF_DEBUG_NO_SHOW') != '1':
        plt.show(block=_os.environ.get('EBMKF_DEBUG_BREAKPOINT') != '1')
    return figures,series


def _diagnostic_main():
    import argparse
    from pathlib import Path
    parser=argparse.ArgumentParser(description='Inspect active EBMKF-ta4 inputs, states, innovations and the inferred extra SW correction; no evaluation caches saved.')
    parser.add_argument('experiment',nargs='?',default='historical',choices=['historical',
        'fut_ESM1-2-LR_SSP126_constVolc','fut_ESM1-2-LR_SSP245_constVolc','fut_ESM1-2-LR_SSP370_constVolc',
        'fut_NorESM_RCP45_Volc','fut_NorESM_RCP45_VolcConst'])
    parser.add_argument('member',nargs='?',type=int,default=0,help='Zero-based ensemble member')
    parser.add_argument('--breakpoint',action='store_true',help='Pause inside run_method after plotting')
    parser.add_argument('--save-dir',type=Path,help='Optional PNG/CSV/parameter report destination')
    parser.add_argument('--no-show',action='store_true',help='Run/save without GUI windows')
    args=parser.parse_args()
    if args.member<0:parser.error('Member must be non-negative')
    destination=args.save_dir.resolve() if args.save_dir else None
    _os.chdir(_V11Path(__file__).resolve().parents[3])
    if args.experiment=='historical':
        frame=pd.read_csv('Common_Data/HadCRUT.5.1.csv')
        years=frame.Time.to_numpy();temperature=frame.Anomaly.to_numpy()
        uncertainty=(frame.Lower.to_numpy(),frame.Upper.to_numpy())
    else:
        import evaluate_future_scenarios as future
        future.exp_attr=args.experiment.split('_');arrays=future.collect_data(future.exp_attr)
        if args.member>=len(arrays[0]):parser.error('Member outside available ensemble')
        temperature,start,years,sim,_=future.select_data(future.years_past[-1],args.member,*arrays)
        n=len(future.years_past)
        uncertainty=tuple(np.concatenate((temperature[:n]+ci-future.temps_obs_past,
                           sim[start:]-np.mean(future.temps_obs_past[-10:]-ci[-10:])))
                          for ci in [future.temps_CIl_past,future.temps_CIu_past])
    if not len(years)==len(temperature)==len(uncertainty[0])==len(uncertainty[1]):
        raise ValueError('Inconsistent diagnostic input lengths')
    _os.environ['EBMKF_DEBUG_PLOTS']='1'
    _os.environ['EBMKF_DEBUG_BREAKPOINT']='1' if args.breakpoint else '0'
    _os.environ['EBMKF_DEBUG_NO_SHOW']='1' if args.no_show else '0'
    if destination:_os.environ['EBMKF_DEBUG_SAVE_DIR']=str(destination)
    result=run_method(years,temperature,uncertainty,args.member,args.experiment,diagnostics=True)
    if not all(np.isfinite(a).all() for a in result):raise ValueError('Non-finite ta4 result')
    print(f'Checked ta4 {args.experiment}, member {args.member}, {years[0]}–{years[-1]}; no evaluation caches or legacy extras overwritten.')


if __name__=='__main__':
    _diagnostic_main()
