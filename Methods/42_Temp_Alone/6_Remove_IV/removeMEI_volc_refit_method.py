import numpy as np
import statsmodels.api as sm
import pandas as pd

enso = pd.read_csv("./Common_Data/meiv_shift.csv")
omega = 2 * np.pi / 1 
data = np.genfromtxt(open("./Common_Data/toyKFmodelData8c.csv", "rb"),dtype=float, delimiter=',')


def run_method(years, temperature, uncert, model_run, experiment_type):
    regLinX = np.arange(1935,1970+1)
    regY = temperature[regLinX-1850]
    constV = np.var(regY)
    
    cos_2yr = np.cos(0.5 * omega * years)
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    means2 = empser.copy()
    ses2 = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    avg_len_l=10
    avg_len_u=11 #actually 10 yrs below, that year, 10 yrs after
    means21yr = empser.copy()
    for i in range(avg_len_l, len(years) - avg_len_u+1):
        chunk=temperature[i-avg_len_l:i+avg_len_u]
        means21yr[i] = np.mean(chunk)


    for endi in range(50, len(years) ):
        MEI_shift = 22
        exstart = -avg_len_l +MEI_shift 
        X = pd.DataFrame({
            'MEI': enso['AVG'][avg_len_l -MEI_shift +exstart +1 : (endi-avg_len_u+2 - MEI_shift)], #1971 start
            'AOD': data[avg_len_l +exstart        : (endi-avg_len_u+1), 3], 
            'TSI': data[avg_len_l +exstart        : (endi-avg_len_u+1) , 8],
            'cos2': cos_2yr[avg_len_l +exstart    : (endi-avg_len_u+1)], #linear fit
        })
        # Add a constant term for the intercept
        X = sm.add_constant(X)
        y = temperature[avg_len_l +exstart : (endi-avg_len_u+1)] - means21yr[avg_len_l +exstart : (endi-avg_len_u+1)]
        #here we are predicting the differences from the 21yr running mean
        model = sm.OLS(y, X).fit()

        X2 = pd.DataFrame({
            'MEI': enso['AVG'][1:(endi-MEI_shift+1)],
            'AOD': data[MEI_shift:endi, 3],
            'TSI': data[MEI_shift:endi, 8],
            'cos2': cos_2yr[MEI_shift:endi],
        })

        X2 = sm.add_constant(X2)
        pred2=model.get_prediction(X2)
        #print(pred2[endi-MEI_shift])
        means[endi]= temperature[endi]- pred2.predicted_mean[endi-MEI_shift -1 ] 
        ses[endi] = np.sqrt( temps_1std[endi]**2 + (pred2.se[endi-MEI_shift -1])**2+ constV)
    #print(pred2.predicted_mean)
    means2[MEI_shift:174]= temperature[MEI_shift:174] - pred2.predicted_mean
    ses2[MEI_shift:174] = np.sqrt( temps_1std[MEI_shift:174]**2 + (pred2.se)**2 + constV/4)
    return means, ses,means2, ses2

