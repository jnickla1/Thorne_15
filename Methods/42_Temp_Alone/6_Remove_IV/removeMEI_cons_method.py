import numpy as np
import statsmodels.api as sm
import pandas as pd
enso = pd.read_csv("./Common_Data/meiv_shift.csv")



def run_method(years, temperature, uncert, model_run, experiment_type):
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
##    X = pd.DataFrame({
##        'MEI': enso['AVG'][0:],
##        'AOD': data[21:, 3],
##        'TSI': data[21:, 8],
##        'cos2': cos_2yr[21:],
##        'time': years[21:]-1971,
##    })
    coefMEI= 0.1 #from Foster Rahmstorf 2011 http://dx.doi.org/10.1088/1748-9326/6/4/044022
    uncertMEI = 0.25 #also from above paper
    means[21:174]= temperature[21:174]-enso['AVG'][0:174]*coefMEI
    ses[21:174] = np.sqrt( temps_1std[21:174]**2 + (enso['AVG'][0:174]*coefMEI)**2)
    
    return means, ses, empser.copy(), empser.copy()

