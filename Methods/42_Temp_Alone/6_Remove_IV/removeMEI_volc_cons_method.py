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
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
#here we're just fitting a linear line
    X = pd.DataFrame({
        'MEI': enso['AVG'][108:], #1971 start
        'AOD': data[129:, 3],
        'TSI': data[129:, 8],
        'cos2': cos_2yr[129:174],
        'time': years[129:174]-1971,
    })
    # Add a constant term for the intercept
    X = sm.add_constant(X)
    y = temperature[129:174]
    #here we are trying to predict temperature
    model = sm.OLS(y, X).fit()
    st_idx = 1971- 1871
    X2 = pd.DataFrame({
        'MEI': enso['AVG'][1+st_idx:],
        'AOD': data[22+st_idx:, 3],
        'TSI': data[22+st_idx:, 8],
        'cos2': cos_2yr[22+st_idx:174],
        'time': np.zeros(174- (22+st_idx)),
    })

    X2 = sm.add_constant(X2)
    pred2=model.get_prediction(X2) #natural variablity component
   # print(pred2.predicted_mean)
    offset = np.mean(pred2.predicted_mean)
    means[22+st_idx:174]= temperature[22+st_idx:174] - pred2.predicted_mean +offset
    ses[22+st_idx:174] = np.sqrt( temps_1std[22+st_idx:174]**2 + (pred2.se)**2 + constV/4)
    
    return means, ses, empser.copy(), empser.copy()

