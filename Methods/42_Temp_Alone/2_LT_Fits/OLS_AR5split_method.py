import numpy as np
from scipy import stats

#AR5 Chapt 2 pg 180 box 2.2
#1901–1950:
#0.107 ± 0.026 °C per decade (90%CI)
#1951–2012
#0.106 ± 0.027

def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor



def run_method(years, temperature, uncert, model_run, experiment_type):
##FIRST HALF
    empser = np.full(np.shape(years),np.nan)
    slope = 0.107/10
    slope_assessed_err = 0.026 / 10

    regX = np.arange(1901,1950+1)
    n=len(regX)
    regY = temperature[regX-1850]
    intercept = np.mean(regY - regX*slope)
    y_pred = slope * regX + intercept

    residuals = regY - y_pred
    
    # Calculate standard error of the estimate (s_e)
    s_e = np.sqrt(np.sum(residuals**2) / (n - 2))

    # Critical t-value for 90% confidence level
    t_crit = stats.t.ppf(0.90, df=n-2)
    t_crit2 = stats.t.ppf(.5+0.3413, df=n-2)  #one standard error
    SE_m = s_e / np.sqrt(np.sum((regX - np.mean(regX))**2))


    adjustment_factor = slope_assessed_err / (t_crit * SE_m)
    
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)

    means[regX-1850] = y_pred
    ses[regX-1850] = confidence_interval(regX, y_pred, np.mean(regX), s_e, t_crit2, n, factor=adjustment_factor)

###SECOND HALF
    slope = 0.106/10
    slope_assessed_err = 0.027 / 10
    regX = np.arange(1951,2012+1)
    n=len(regX)
    regY = temperature[regX-1850]
    intercept = np.mean(regY - regX*slope)
    y_pred = slope * regX + intercept

    residuals = regY - y_pred
    
    # Calculate standard error of the estimate (s_e)
    s_e = np.sqrt(np.sum(residuals**2) / (n - 2))

    # Critical t-value for 90% confidence level
    t_crit = stats.t.ppf(0.90, df=n-2) #one standard error
    t_crit2 = stats.t.ppf(.5+0.3413, df=n-2)
    SE_m = s_e / np.sqrt(np.sum((regX - np.mean(regX))**2))


    adjustment_factor = slope_assessed_err / (t_crit * SE_m)

    regX2 = np.arange(1951,years[-1]+1)
    y_pred2=slope * regX2 + intercept
    means[regX2-1850] = y_pred2
    ses[regX2-1850] = confidence_interval(regX2, y_pred2, np.mean(regX), s_e, t_crit2, n, factor=adjustment_factor)


    return empser.copy(), empser.copy(), means, ses

