import numpy as np
from scipy import stats

import matplotlib.pyplot as plt


def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor



def run_method(years, temperature, uncert, model_run, experiment_type):
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    regLinX = np.arange(1935,1970+1)
    regY = temperature[regLinX-1850]
    means[regLinX-1850] = np.mean(regY)
    nconst = len(regLinX)
    constV = np.var(regY)
    ses[regLinX-1850] = np.sqrt(constV /20)

    regLinXt = np.arange(1970,1975) - 1850
    means[regLinXt] = np.mean(regY) + (regLinXt+1850-1970)*0.11/5.5
    ses[regLinXt] = np.sqrt(constV /20)
    
    for endd in range(1975,years[-1]+1):
        #starting in 1901
        regX = np.arange(1970,endd+1)
        n=len(regX)
        regY = temperature[regX-1850]
        regres = stats.theilslopes(regY, regX, 0.3413*2) # +- 1 se
        slope = regres.slope
        intercept =regres.intercept
        sstdrr = (regres.high_slope - regres.low_slope)/2
        y_pred = slope * regX + intercept

        means[endd-1850] = y_pred[-1]
        ses[endd-1850] = np.sqrt( (n/2*sstdrr)**2 + constV /n)




    return means, ses

