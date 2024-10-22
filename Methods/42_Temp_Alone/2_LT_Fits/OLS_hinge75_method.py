import numpy as np
from scipy import stats

#import matplotlib.pyplot as plt


def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor



def run_method(years, temperature, uncert, model_run, experiment_type):
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    regLinX = np.arange(1935,1975+1)
    regY = temperature[regLinX-1850]
    means[regLinX-1850] = np.mean(regY)
    ses[regLinX-1850] = np.sqrt(np.var(regY)/20)
    
    for endd in range(1976,years[-1]+1):
        #starting in 1901
        regX = np.arange(1974,endd+1)
        n=len(regX)
        regY = temperature[regX-1850]
        regres = stats.linregress(regX, regY) #can also give stdrr of slope, intercept
        slope = regres.slope
        intercept =regres.intercept
        sstdrr = regres.stderr
        y_pred = slope * regX + intercept
        residuals = regY - y_pred
        # Calculate standard error of the estimate (s_e)
        s_e = np.sqrt(np.sum(residuals**2) / (n - 2))

        # Critical t-value for 1se confidence level
        t_crit = stats.t.ppf(.5 + 0.3413, df=n-2) #one standard error, two sided

        means[endd-1850] = y_pred[-1]
        cis = confidence_interval(regX, y_pred, np.mean(regX), s_e, t_crit, n )
        ses[endd-1850] = cis[-1]

##        if(endd ==years[-1]):
##            plt.plot(regX, y_pred)
##            plt.fill_between(regX, y_pred+cis, y_pred-cis)
##            plt.plot(regX, (slope) * regX + intercept + (sstdrr) * (regX-np.mean(regX)), 'r')
##            plt.plot(regX, (slope) * regX + intercept- (sstdrr) * (regX-np.mean(regX)), 'r')
##            plt.figure()

    means2= means.copy()
    ses2 = ses.copy()
    means2[1974-1850:] = y_pred
    ses2[1974-1850:] = cis


    return means, ses, means2,ses2
