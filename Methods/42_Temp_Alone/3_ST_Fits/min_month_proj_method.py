import numpy as np
from scipy import stats
import pandas as pd
import os
import pdb

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])

def max_every_ny(lst, n, y):
    """Calculates the average of every n elements in a list."""
    return np.array([np.max(lst[max(0,i-n*y):i + n]) for i in range(0, len(lst), n)])

def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor

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
    offset50 = np.mean(temperature[0:50]) - np.mean(temperature_mon[0:50*12])
    temperature_mon_uncert = (monthly["Upper confidence limit (97.5%)"] - monthly["Lower confidence limit (2.5%)"]).to_numpy()
    months = np.arange(1850, 1850 + len(temperature_mon)/12, 1/12)
    
    if (experiment_type!="historical"):
        return empser, empser, empser, empser #not valid for future runs, return blanks

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
        monthly_var[i] = (sstdrr * nmon_offset)**2 + (nmon_min_max/3.6 * slope)**2 + (temperature_mon_uncert[i])**2 /12

    return max_every_ny( monthly_mean, 12,4) + offset50, 2*np.sqrt( average_every_n( monthly_var, 12)), empser.copy(), empser.copy()
