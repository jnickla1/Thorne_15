import numpy as np
from scipy import stats




def run_method(years, temperature, uncert, model_run, experiment_type):
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    regLinX = np.arange(1935,1975+1)
    regY = temperature[regLinX-1850]
    intercept =  np.mean(regY)
    means[regLinX-1850] = intercept
    
    ses[regLinX-1850] = np.sqrt(np.var(regY)/20)
    
    for endd in range(1976,years[-1]+1):
        #starting in 1901
        regX = np.arange(1974,endd+1)
        n=len(regX)
        regY = temperature[regX-1850]

        a, res, _, _ = np.linalg.lstsq(np.reshape((regX - regX[0]), (-1, 1)), regY-intercept, rcond=None) #make first arg 2d
        
        y_pred = a * (regX- regX[0]) + intercept
        # Calculate standard error of the estimate (s_e)
        s_e = np.sqrt(res / (n - 1))
        SE_m = s_e / np.sqrt(np.sum((regX - np.mean(regX))**2))
        
        # Critical t-value for 1 se confidence level
        #t_crit = stats.t.ppf(.5+0.3413, df=n-1) #one standard error

        means[endd-1850] = y_pred[-1]
        ses[endd-1850] = n*SE_m


    means2= means.copy()
    ses2 = ses.copy()
    ns2 = np.arange(1,years[-1]+2-1974)
    means2[1974-1850:] = y_pred
    ses2[1974-1850:] = ns2*SE_m


    return means, ses, means2,ses2

