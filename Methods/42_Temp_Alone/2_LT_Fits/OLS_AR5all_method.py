import numpy as np
from scipy import stats

#AR5 Chapt 2 pg 180 box 2.2
#1901–2012:
#0.075 ± 0.013 °C per decade (90%CI)

def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor



def run_method(years, temperature, uncert, model_run, experiment_type):

    slope = 0.075/10
    slope_assessed_err = 0.013 / 10
    
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs

    regX = np.arange(1901,2013)
    n=len(regX)
    regY = temperature[regX-1850]

    intercept = np.mean(regY - regX*slope)
    y_pred = slope * regX + intercept

    residuals = regY - y_pred
    
    # Calculate standard error of the estimate (s_e)
    s_e = np.sqrt(np.sum(residuals**2) / (n - 2))

    # Critical t-value for 90% confidence level
    t_crit = stats.t.ppf(0.90, df=n-2)
    SE_m = s_e / np.sqrt(np.sum((regX - np.mean(regX))**2))


    adjustment_factor = slope_assessed_err / (t_crit * SE_m)
    

    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    t_crit2 = stats.t.ppf(.5 + 0.3413, df=n-2) #one standard error, two sided
    regX2 = np.arange(1901,years[-1]+1)
    y_pred2=slope * regX2 + intercept
    means[regX2-1850] = y_pred2
    ses[regX2-1850] = confidence_interval(regX2, y_pred2, np.mean(regX), s_e, t_crit2, n, factor=adjustment_factor)



    return means, ses

