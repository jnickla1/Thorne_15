import numpy as np
from scipy import stats


def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor



def run_method(years, temperature, uncert, model_run, experiment_type):
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    for endd in range(1950,years[-1]+1):
        #starting in 1901
        regX = np.arange(1901,endd+1)
        n=len(regX)
        regY = temperature[regX-1850]
        slope, intercept, _, _, _ = stats.linregress(regX, regY) #can also give stdrr of slope, intercept
        y_pred = slope * regX + intercept
        residuals = regY - y_pred
        # Calculate standard error of the estimate (s_e)
        s_e = np.sqrt(np.sum(residuals**2) / (n - 2))

        # Critical t-value for 1 se confidence level
        t_crit = stats.t.ppf(.5+0.3413, df=n-2) #one standard error

        means[endd-1850] = y_pred[-1]
        cis = confidence_interval(regX, y_pred, np.mean(regX), s_e, t_crit, n )
        ses[endd-1850] = cis[-1]

    means2= means.copy()
    ses2 = ses.copy()
    means2[1901-1850:] = y_pred
    ses2[1901-1850:] = cis

    return means, ses, means2,ses2

