import numpy as np
from scipy import stats
import pandas as pd
import os
import pdb
from fut_evaluation_script import collect_data, gen_orig_number
import config
from netCDF4 import Dataset

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])

def max_every_ny(lst, n, y):
    """Calculates the average of every n elements in a list."""
    return np.array([np.max(lst[max(0,i-n*y):i + n]) for i in range(0, len(lst), n)])

def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor

def seasonal_cycle(arr):
    arr = np.asarray(arr)
    n_years = arr.size // 12
    return arr[:n_years*12].reshape(n_years, 12).mean(axis=0)

def run_method(years, temperature, uncert, model_run, experiment_type):

    
    len_trend=15 #taking 15 yr trend in data
    nmon_offset = 33/12 #from Cannon 2025 https://www.nature.com/articles/s41558-025-02247-8
    nmon_min_max = (76 - (-28))/12 #90% confidence range
    
    empser = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    cur_path = os.path.dirname(os.path.realpath(__file__))
    
    monthly  = pd.read_csv(cur_path+"/../../../Common_Data/HadCRUT.5.0.2.0.analysis.summary_series.global.monthly.csv")


    
    temperature_mon = monthly["Anomaly (deg C)"].to_numpy()
    offset50 = np.mean(temperature[0:50]) - np.mean(temperature_mon[0:(50*12)])
    temperature_mon_uncert = (monthly["Upper confidence limit (97.5%)"] - monthly["Lower confidence limit (2.5%)"]).to_numpy()

    
    if (experiment_type!="historical"):
        (sims_tas, stime_mon, sims_tas_hist,  stime_mon_hist) = collect_data(experiment_type.split("_"))

        this_hsim_yr = average_every_n(sims_tas_hist[model_run,:], 12)
        this_sim = sims_tas[model_run,:] #keep monthly
        this_hsim = sims_tas_hist[model_run,:]
        

        #offset_sim = np.mean( this_sim_yr[0:start_sim] - temps_obs_past[-start_sim:])
        #must decide how to offset the simulation - can do it so the first 50 yrs are 0 as we did for obs

        exp_attr = experiment_type.split("_")
        if (exp_attr[1]=='ESM1-2-LR'):
            offsetn = np.mean(temperature_mon[0:(50*12)])- np.mean(this_hsim_yr[0:50]) #preindustrial baseline
            temperature_mon_wcycle = np.concatenate((this_hsim,this_sim))+offsetn

            
        elif (exp_attr[1]=='NorESM'):
            #replacing temps_obs_past
            long_past_index = ((gen_orig_number(model_run,np.shape(sims_tas)[0])-1 )// 20) #either 1, 2, or 3, still in right order
            long_past_data_loc = config.CLIMATE_DATA_PATH+'/NorESM_volc/NorESM1-M-historical/hist_aave_tas.nc'
            variable = Dataset(long_past_data_loc, 'r').variables['tas']
            long_past_tas_array = variable[:].__array__()
            long_past_tas = long_past_tas_array[long_past_index,:]
            temps_obs_past = average_every_n(temperature_mon,12)
            offsetn =  -np.mean( this_hsim_yr[0:20]) + np.mean(temps_obs_past[(-1850+1980):(-1850+1980+20)])
            # baseline 1980-2000 match
            temperature_mon_wcycle = np.concatenate((long_past_tas[0:((-1850+1980)*12)],this_hsim,this_sim))+offsetn

        #HAVE TO REMOVE SEASONAL CYCLE
        cycle = seasonal_cycle(temperature_mon_wcycle[0:(12*100)])
        cycle_balanced = cycle - np.mean(cycle)
        temperature_mon = temperature_mon_wcycle - np.tile(cycle_balanced,len(years))
        
    months = np.arange(1850, 1850 + len(temperature_mon)/12, 1/12)
    monthly_mean =  np.full(np.shape(years)[0]*12,np.nan)
    monthly_var =  np.full(np.shape(years)[0]*12,np.nan)
    for i in range((len_trend) * 12 -1, (len(years)) * 12): #run for every month, don't yet have full 2025 but not provided
        # Fit the OLS trendline for each 15-year chunk
        regX = months[(i-(len_trend)*12) + 1:(i+1)]
        regY = temperature_mon[(i-(len_trend)*12) + 1:(i+1)]
        regres = stats.linregress(regX, regY) #can also give stdrr of slope, intercept
        slope = regres.slope
        intercept =regres.intercept
        sstdrr = regres.stderr


        monthly_mean[i] = np.min( temperature_mon[ i- 12 + 1:(i+1)])  + nmon_offset * slope

        monthly_var[i] = (sstdrr * nmon_offset)**2 + (nmon_min_max/3.6 * slope)**2 + (temperature_mon_uncert[min(i,len(temperature_mon_uncert)-1)])**2 /12

    return max_every_ny( monthly_mean, 12,4) + offset50, 2*np.sqrt( average_every_n( monthly_var, 12)), empser.copy(), empser.copy()
