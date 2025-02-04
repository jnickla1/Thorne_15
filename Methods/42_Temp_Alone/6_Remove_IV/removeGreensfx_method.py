import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
import os
import netCDF4
import pdb

def run_method(years, temperature, uncert, model_run, experiment_type):

    empser = np.full(np.shape(years),np.nan)
    
    cur_path = os.path.dirname(os.path.realpath(__file__))

    #fut_ESM1-2-LR_SSP126_constVolc
    exp_attr = experiment_type.split("_")
    temp_mean = np.mean(temperature)
    comput_uncert = np.full(np.shape(temperature),0.06)

    try: #unclear why I'm using a try..except here
        if experiment_type=="historical":
            filein = netCDF4.Dataset(cur_path+"/filterAna_obs_HadCRUT5_"+experiment_type+".nc",'r')
            #preind_base = np.mean(temps_obs[0:50])
            comput_temps = np.array(filein.variables['tas_fbr_aa'])
            comput_temps2 = empser.copy()
            comput_temps2[0:174] = comput_temps
        
        elif exp_attr[1]=="ESM1-2-LR":
            filein = netCDF4.Dataset(cur_path+"/filterAna_results_MPI"+exp_attr[2][3:]+"_1850_2100_v1.nc",'r')
            comput_temps_full = np.array(filein.variables['tas_fbr_aa'])
            comput_temps = comput_temps_full[:,model_run]
            comput_mean = np.mean(comput_temps)
            comput_temps2 = empser.copy()
            comput_temps2[(years[0]-1850):(years[0]-1850+len(comput_temps))] = comput_temps + temp_mean - comput_mean

        elif exp_attr[1]=="NorESM":
            filein = netCDF4.Dataset(cur_path+"/filterAna_results_v"+exp_attr[3][1:]+"_1980_2099_v4.nc",'r')
            comput_temps_full = np.array(filein.variables['tas_fbr_aa'])
            comput_temps = comput_temps_full[:,model_run]
            comput_mean = np.mean(comput_temps)
            temp_mean = np.mean(temperature[(years[0]-1980):(years[0]-1980+len(comput_temps))])
            comput_temps2 = empser.copy()
            comput_temps2[(years[0]-1980):(years[0]-1980+len(comput_temps))] = comput_temps  + temp_mean - comput_mean
            

        #print(means[-1]+preind_base)
        return comput_temps2, comput_uncert, empser, empser
    except:
        return empser, empser, empser, empser

