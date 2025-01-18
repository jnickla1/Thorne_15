import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
import subprocess
import pdb;
import os
import sys 

def run_method(years, temperature, uncert, model_run, experiment_type):

    dir_path = os.path.dirname(os.path.realpath(__file__))
    load_file = dir_path+"/Bayes_Sequential_Change_Point/output/"+experiment_type+str(model_run)+"BSCtemperatures.csv"
    
    if not(os.path.exists(load_file)):
        (temps_CIl, temps_CIu) = uncert
        #save to input_data
        loc_string = experiment_type+str(model_run)
        inp_save_path = dir_path+"/Bayes_Sequential_Change_Point/input_data/"+loc_string+".csv"

        # Define the headers in the dataframe
        df = pd.DataFrame({"Time" : years, "Anomaly" : temperature, "Lower" : temps_CIl, "Upper" : temps_CIu})
        # Save the array as a CSV file with headers
        df.to_csv(inp_save_path, index=False)
                   
        #call matlab function, which will save output 
    
        matlab_command = ( ###MIGHT NEED TO CHANGE YOUR VERSION OF MATLAB
            "/Applications/MATLAB_R2019a.app/bin/matlab -nodisplay -nosplash -nodesktop "
            "-r \"cd "+dir_path+"/Bayes_Sequential_Change_Point/; "
            "Bayes_changepoint_seq_iter('"+loc_string+"'); exit;\""
        )

        # Run the command using subprocess
        try:
            subprocess.run(matlab_command, shell=True, check=True)
            print("MATLAB command executed successfully.")
        except subprocess.CalledProcessError as e:
            sys.exit(f"Error while executing MATLAB command: {e}")
            

    comput_temps = pd.read_csv(load_file, sep=',', header=None)
    comput_uncert = pd.read_csv(dir_path+"/Bayes_Sequential_Change_Point/output/"+
                                experiment_type+str(model_run)+"BSCuncertainty.csv", sep=',', header=None)
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)

    stidx = 49 #same as in the Matlab File
    for endi in range(stidx, len(years)):

        means[endi] = comput_temps.values[endi-stidx, endi]  #used to subtract preind base
        ses[endi] = comput_uncert.values[endi-stidx, endi]
    #print(means[-1]+preind_base)
    return means, ses, comput_temps.values[endi-stidx, :] , comput_uncert.values[endi-stidx, :]

