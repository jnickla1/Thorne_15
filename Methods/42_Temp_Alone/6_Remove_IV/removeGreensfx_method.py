import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
import os
import netCDF4

def run_method(years, temperature, uncert, model_run, experiment_type):

    empser = np.full(np.shape(years),np.nan)
    
    cur_path = os.path.dirname(os.path.realpath(__file__))
    try:
        filein = netCDF4.Dataset(cur_path+"/FilterAna_obs_HadCRUT5_"+experiment_type+".nc",'r')
        #preind_base = np.mean(temps_obs[0:50])
        comput_temps = np.array(filein.variables['tas_fbr_aa'])
        comput_uncert = np.full(np.shape(temperature),0.06)
        
        comput_temps2 = empser.copy()
        comput_temps2[0:174] = comput_temps

        #print(means[-1]+preind_base)
        return comput_temps2, comput_uncert, empser, empser
    except:
        return empser, empser, empser

