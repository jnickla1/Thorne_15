import numpy as np
from scipy import stats

def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor


def run_method(years, temperature, uncert, model_run, experiment_type):
    len_trend=15
    empser = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    
    for i in range(len_trend-1, len(years)):
        # Fit the OLS trendline for each 15-year chunk : 14yrs behind and this current year
        regX = years[i-len_trend+1:i+1]
        regY = temperature[i-len_trend+1:i+1]
        regres = stats.linregress(regX, regY) #can also give stdrr of slope, intercept
        slope = regres.slope
        intercept =regres.intercept
        sstdrr = regres.stderr
        y_pred = slope * regX + intercept
        residuals = regY - y_pred
        # Calculate standard error of the estimate (s_e)
        s_e = np.sqrt(np.sum(residuals**2) / (len_trend - 2))

        # Critical t-value for 1se confidence level
        t_crit = stats.t.ppf(.5 + 0.3413, df=len_trend-2) #one standard error, two sided

        means[i] = y_pred[-1]
        cis = confidence_interval(regX, y_pred, np.mean(regX), s_e, t_crit, len_trend )
        ses[i] = cis[-1]
    means2 = empser.copy()
    ses2 = empser.copy()
    for i in range(len_trend-1, len(years)-len_trend):
        # Fit the OLS trendline for each 15-year chunk : 14yrs behind and this current year and 15yrs after
        regX = years[i-len_trend+1:i+1+len_trend]
        regY = temperature[i-len_trend+1:i+1+len_trend]
        regres = stats.linregress(regX, regY) #can also give stdrr of slope, intercept
        slope = regres.slope
        intercept =regres.intercept
        sstdrr = regres.stderr
        y_pred = slope * regX + intercept
        residuals = regY - y_pred
        # Calculate standard error of the estimate (s_e)
        s_e = np.sqrt(np.sum(residuals**2) / (len_trend - 2))

        # Critical t-value for 1se confidence level
        t_crit = stats.t.ppf(.5 + 0.3413, df=len_trend-2) #one standard error, two sided

        means2[i] = y_pred[len_trend-1]
        cis = confidence_interval(regX, y_pred, np.mean(regX), s_e, t_crit, len_trend )
        ses2[i] = cis[len_trend-1]
        
    return means, ses, means2, ses2
