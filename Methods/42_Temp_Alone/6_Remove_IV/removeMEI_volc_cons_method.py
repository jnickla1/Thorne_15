import numpy as np
import statsmodels.api as sm
import pandas as pd
from netCDF4 import Dataset
import os
import pdb

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])


def run_method(years, temperature, uncert, model_run, experiment_type):
    regLinX = np.arange(1935,1970+1)
    regY = temperature[regLinX-1850] #all inputs start in 1850?
    constV = np.var(regY)
    omega = 2 * np.pi / 1
    cos_2yr = np.cos(0.5 * omega * years)
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs

    enso = pd.read_csv("./Common_Data/meiv_shift.csv")
    data = np.genfromtxt(open("./Common_Data/toyKFmodelData8c.csv", "rb"),dtype=float, delimiter=',')
    AODdata = data[:, 3]  #these two may get overwritten
    ensoA = enso['AVG']   #these two may get overwritten

    if experiment_type == "historical":
        TSIdata = data[:, 8] #these two may get overwritten
        start_shift = 1
        s_yr = 21 # could use this with NaN from file, +1 to account for month shift to start in August/Sept
        start_yr = start_shift+s_yr  #22 start in 1871
        end_yr = 174 #others are left as NaN


    else:
        exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126_constVolc #fut_NorESM_RCP45_Volc
        new_iter=len(years)
        given_preind_base = np.mean(temperature[0:50])

        dataOpenTSI = pd.read_csv(os.path.expanduser('~/')+"climate_data/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp126_1750-2500.csv")
        TSIdata = dataOpenTSI['solar'][100:(100+len(years))]
        

        if (exp_attr[1]=='ESM1-2-LR'):
            enso_data =Dataset(os.path.expanduser('~/')+"climate_data/ESM1-2-LR/combined/"+exp_attr[2].lower()+"_nino34_aave_tas.nc", 'r').variables['__xarray_dataarray_variable__']
            enso_arr = enso_data[:].__array__()
            ensoA = average_every_n(enso_arr[model_run,6:(-12+6)], 12) #start in August
            start_yr=1 #everyone else starts in 1851
            start_shift = 0 #ensoA now starts in 1850
            end_yr = len(ensoA)+start_yr #should be whole dataset
            AODdata0 = AODdata
            involcavg = np.mean(1/(AODdata0+9.7279))
            AODdata= np.full(len(years), (1/involcavg-9.7279))
            AODdata[0:len(AODdata0)] =AODdata0 
                                       
            
            
        elif (exp_attr[1]=='NorESM'):
            enso_data =Dataset(os.path.expanduser('~/')+"climate_data/NorESM_volc/BethkeEtAl2017/"+exp_attr[2].lower()+exp_attr[3]+"_nino34_tas.nc", 'r').variables['__xarray_dataarray_variable__']
            enso_arr = enso_data[:].__array__()
            ensoNew =np.zeros(len(TSIdata))
            ensoNew[22:(22+len(ensoA))]=ensoA
            ensofromNoresm = average_every_n(enso_arr[model_run,6:], 12) #start in August
            ensoNew[(1980-1850) :(1980-1850 + len(ensofromNoresm))]  = ensofromNoresm
            ensoA = ensoNew #create new hybrid record to allow for fitting window
            start_yr=0
            start_shift = 0
            end_yr = len(ensoA)+start_yr #should be whole dataset
            AOD_simD =Dataset(os.path.expanduser('~/')+"climate_data/NorESM_volc/BethkeEtAl2017/"+exp_attr[2].lower()+exp_attr[3]+"_aod.nc", 'r').variables['__xarray_dataarray_variable__']
            AOD_simA = AOD_simD[:].__array__()
            AOD_sim = average_every_n(AOD_simA[model_run,:], 12)
            AODdata0 = AODdata
            AODdata= np.zeros(len(years))
            AODdata[0:len(AODdata0)] =AODdata0
            AODdata[(1980-1850):(1980-1850+len(AOD_sim))] = AOD_sim
    
#here we're just fitting a linear line
    st_idx0 = 1979- 1850-start_yr #107

    X = pd.DataFrame({
        'MEI': ensoA[start_shift+st_idx0:], #1871 , start in 1979
        'AOD': AODdata[start_yr+st_idx0:],
        'TSI': TSIdata[start_yr+st_idx0:],
        'cos2': cos_2yr[(start_yr+st_idx0):end_yr],
        'time': years[(start_yr+st_idx0):end_yr]-1979, #time fit doesnt really matter
    })
    # Add a constant term for the intercept
    X = sm.add_constant(X)
    y = temperature[(start_yr+st_idx0):end_yr]
    #here we are trying to predict temperature 1971-2023
    
    model = sm.OLS(y, X).fit()
    st_idx = 1971- 1850-start_yr #100
    X2 = pd.DataFrame({
        'MEI': ensoA[start_shift+st_idx:],
        'AOD': AODdata[start_yr+st_idx:],
        'TSI': TSIdata[start_yr+st_idx:],
        'cos2': cos_2yr[start_yr+st_idx:end_yr],
        'time': np.zeros(end_yr- (start_yr+st_idx)),
    })

    X2 = sm.add_constant(X2)
    pred2=model.get_prediction(X2) #natural variablity component
   # print(pred2.predicted_mean)
    offset = np.mean(pred2.predicted_mean)
    means[start_yr+st_idx:end_yr]= temperature[start_yr+st_idx:end_yr] - pred2.predicted_mean +offset
    ses[start_yr+st_idx:end_yr] = np.sqrt( temps_1std[start_yr+st_idx:end_yr]**2 + (pred2.se)**2 + constV/4)
    #only making predictions from 1971 onwards.
    
    return means, ses, empser.copy(), empser.copy()

