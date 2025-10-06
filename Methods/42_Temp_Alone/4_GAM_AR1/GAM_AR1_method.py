import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
import os
import sys
import subprocess
def run_method(years, temperature, uncert, model_run, experiment_type):

    preind_base =  0 #np.mean(temps_obs[0:50])
    cur_path = os.path.dirname(os.path.realpath(__file__))
    load_file = cur_path+"/GAM_AR_Stephenson/output/gamAR1_fits_"+experiment_type+str(model_run)+".csv"

    if not(os.path.exists(load_file)):
        (temps_CIl, temps_CIu) = uncert
        #save to input_data
        loc_string = experiment_type+str(model_run)
        inp_save_path =  cur_path+"/GAM_AR_Stephenson/input_data/"+loc_string+".csv"

        # Define the headers in the dataframe
        df = pd.DataFrame({"Time" : years, "Anomaly" : temperature, "Lower" : temps_CIl, "Upper" : temps_CIu})
        # Save the array as a CSV file with headers
        df.to_csv(inp_save_path, index=False)
                   
        #call Rfunction, which will save output 
    
        R_command = ("Rscript --vanilla "+cur_path+"/GAM_AR_Stephenson/trend_gamAR1_iter.R "+loc_string)

        # Run the command using subprocess
        try:
            subprocess.run(R_command, shell=True, check=True)
            print("R command executed successfully.")
        except subprocess.CalledProcessError as e:
            sys.exit(f"Error while executing R command: {e}")

            
    comput_temps = pd.read_csv(load_file, sep=',') #has header
    comput_uncert = pd.read_csv(cur_path+"/GAM_AR_Stephenson/output/gamAR1_se_fits_"+experiment_type+str(model_run)+".csv", sep=',')
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)

    stidx = 29
    for endi in range(stidx, len(years)):

        means[endi] = comput_temps.values[endi-stidx, endi] -preind_base 
        ses[endi] = comput_uncert.values[endi-stidx, endi]
    #print(means[-1]+preind_base)
    return means, ses, comput_temps.values[endi-stidx, :] -preind_base, comput_uncert.values[endi-stidx, :]

