import numpy as np
import statsmodels.api as sm
import pandas as pd
from sys import path
import os
from netCDF4 import Dataset

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])


def run_method(years, temperature, uncert, model_run, experiment_type):

    enso_mei = pd.read_csv("./Common_Data/meiv_shift.csv")
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    if experiment_type == "historical":
        start_shift = 1
        ensoA = enso_mei['AVG'][start_shift:]
        s_yr = 21 # could use this with NaN from file, +1 to account for month shift to start in August/Sept
        start_yr = start_shift+s_yr 
        end_yr = len(ensoA)+start_yr #others are left as NaN

    else:
        exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126_constVolc #
        new_iter=len(years)
        given_preind_base = np.mean(temperature[0:50])

        

        if (exp_attr[1]=='ESM1-2-LR'):
            enso_data =Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/ESM1-2-LR/combined/"+exp_attr[2].lower()+"_nino34_aave_tas.nc", 'r').variables['__xarray_dataarray_variable__']
            enso_arr = enso_data[:].__array__()
            ensoA = average_every_n(enso_arr[model_run,6:(-12+6)], 12) #start in August
            start_yr=1
            end_yr = len(ensoA)+start_yr #should be whole dataset
            
            
        elif (exp_attr[1]=='NorESM'):
            enso_data =Dataset(os.path.expanduser('~/')+"data/jnickla1/climate_data/NorESM_volc/BethkeEtAl2017/"+exp_attr[2].lower()+exp_attr[3]+"_nino34_tas.nc", 'r').variables['__xarray_dataarray_variable__']
            enso_arr = enso_data[:].__array__()
            ensoA = average_every_n(enso_arr[model_run,6:], 12) #start in August
            start_yr=(1980-1850)
            end_yr = len(ensoA)+start_yr #should be whole dataset
            
        
##        exp_attr = experiment_type.split("_")
##        end_yr=len(years)
##        dir_path = os.path.dirname(os.path.realpath(__file__))
##        path.insert(1, dir_path+'../3_ST_Fits/')
##        from etrend30y_method import run_method as run_method_30y_etrend
##        smoothed_temp,_,_,_ = run_method_30y_etrend(years, temperature, uncert, model_run, experiment_type)
##        
##        
##        if (exp_attr[1]=='ESM1-2-LR'):
##            start_yr = 1
##            data_loc = os.path.expanduser('~/')+'data/jnickla1/climate_data/'+exp_attr[1]+'/combined/'+exp_attr[2].lower()+'_nino34_aave_tas.nc'
##            dataset = Dataset(data_loc, 'r')
##            variable = dataset.variables['__xarray_dataarray_variable__']
##            sims_tas = variable[:].__array__()
##            ensoAuncor = average_every_n(sims_tas[model_run,8:-(12-8)], 12) #start in August/Sept
##            ensoA = ensoAuncor - smoothed_temp
##            
##        elif (exp_attr[1]=='NorESM'):
##            start_yr = 22
##            data_loc = os.path.expanduser('~/')+'data/jnickla1/climate_data/'+exp_attr[1]+'_volc/BethkeEtAl2017/'+exp_attr[2].lower()+exp_attr[3]+'_nino34_aave_tas.nc'
##            dataset = Dataset(data_loc, 'r')
##            variable = dataset.variables['__xarray_dataarray_variable__']
##            sims_tas = variable[:].__array__()
##            ensoAuncor = average_every_n(sims_tas[model_run,8:-(12-8)], 12) #start in August/Sept
##            ensoAlast = ensoAuncor - smoothed_temp[1980-1850:]
##            ensoA = np.concatenate((enso_mei['AVG'][1:],ensoAlast))
    
        

    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs

    coefMEI= 0.1 #from Foster Rahmstorf 2011 http://dx.doi.org/10.1088/1748-9326/6/4/044022
    uncertMEI = 0.25 #also from above paper
    means[start_yr:end_yr]= temperature[start_yr:end_yr]-ensoA*coefMEI
    ses[start_yr:end_yr] = np.sqrt( temps_1std[start_yr:end_yr]**2 + (ensoA*coefMEI)**2)
    
    return means, ses, empser.copy(), empser.copy()

