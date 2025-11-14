import numpy as np
from scipy.stats import linregress
import os
import pandas as pd
import subprocess
import sys
from fut_evaluation_script import gen_orig_number

def run_method(years, temperature, uncert, model_run0, experiment_type):

    means_c = np.full(np.shape(years), np.nan)
    ses_c = np.full(np.shape(years), np.nan)
    means_r = np.full(np.shape(years), np.nan)
    ses_r = np.full(np.shape(years), np.nan)

    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4  # pm2 stdevs
    dir_path = os.path.dirname(os.path.realpath(__file__))

    if experiment_type == "historical":
        # Original historical processing
        lfcs_path = dir_path + "/lfca/ERSST_firstHadCRUT5_s1975.npy"
        #run /lfca/Python/run_lfca_hadcrut5.py if doesn't exist
        data_monthly = pd.read_csv(dir_path + "/../../../Common_Data/HadCRUT.5.0.2.0.analysis.summary_series.global.monthly.csv")
        temps_obs_monthly = data_monthly.loc[:, "Anomaly (deg C)"].to_numpy()
        preind_base = np.mean(temps_obs_monthly[0:50*12])
        
        st_idx = 1970 - 1850
        lfcs = np.real(np.load(lfcs_path))

        for endi in range(st_idx + 5, (2023-1850)):
            month_end_index = (endi-107)*12
            table_index = endi - st_idx
            temps_cropped = temps_obs_monthly[(1957-1850)*12 : (endi)*12]
            
            slope, intercept, r_value, p_value, std_err = linregress(lfcs[table_index, 0:month_end_index], temps_cropped)
            predicted_temps = np.real(slope) * lfcs[table_index, 0:month_end_index] + intercept
            residuals = temps_cropped - predicted_temps
            residual_std_dev = np.std(residuals)

            means_c[endi] = np.mean(predicted_temps[-13:-1])/2 - preind_base/2 + temperature[endi]/2 
            ses_c[endi] = residual_std_dev/np.sqrt(3)

        means_r[(1957-1850):(2022-1850)] = predicted_temps.reshape(-1, 12).mean(axis=1)/2 - preind_base/2 + temperature[(1957-1850):(2022-1850)]/2 
        ses_r[(1957-1850):(2022-1850)] = residual_std_dev/np.sqrt(3)

    else:
        # Future experiment processing
        #Determine baseline period based on experiment type
        exp_attr = experiment_type.split("_")
        if exp_attr[1] == 'ESM1-2-LR':
            model_run = gen_orig_number(model_run0,50)
        elif exp_attr[1] == 'NorESM':
            model_run = gen_orig_number(model_run0,60)
                
        cache_dir = os.path.join(dir_path, "/lfca/Python/lfca_cache")
        cache_file = os.path.join(cache_dir, f"lfca_{experiment_type}_run{model_run}.npy")
        
        # Check if LFCA results exist, if not, compute them
        if not os.path.exists(cache_file):
            print(f"LFCA cache not found for {experiment_type} run {model_run}. Computing...")
            try:
                # Call the LFCA computation script
                lfca_script = dir_path+  "/lfca/Python/run_lfca_climate_models.py"
                result = subprocess.run([sys.executable, lfca_script, experiment_type, str(model_run)], 
                                      capture_output=True, text=True, check=True)
                print("LFCA computation completed.")
            except subprocess.CalledProcessError as e:
                print(f"Error computing LFCA: {e}")
                print(f"STDOUT: {e.stdout}")
                print(f"STDERR: {e.stderr}")
                return means_c, ses_c, means_r, ses_r
        
        # Load the cached LFCA results
        try:
            lfcs_data = np.load(cache_file)
        except:
            print(f"Failed to load LFCA cache file: {cache_file}")
            return means_c, ses_c, means_r, ses_r
        

        
        if exp_attr[1] == 'ESM1-2-LR':
            # Use 1850-1900 as preindustrial baseline for ESM1-2-LR
            preind_start_idx = 0
            preind_end_idx = 50
            analysis_start_year = 1970  # Start analysis from 1970
        elif exp_attr[1] == 'NorESM':
            # Use 1980-2000 baseline for NorESM (matching original offset approach)
            preind_start_idx = 1980 - 1850
            preind_end_idx = 2000 - 1850
            analysis_start_year = 1990  # Start analysis from 1990
        else:
            print(f"Unknown model type: {exp_attr[1]}")
            return means_c, ses_c, means_r, ses_r
        
        # Calculate baseline offset
        preind_temp_base = np.mean(temperature[preind_start_idx:preind_end_idx])
        
        analysis_start_idx = analysis_start_year - 1850
        
        # Process each year from analysis start to end
        for endi in range(analysis_start_idx, len(years)):
            current_year = years[endi]
            
            # Find corresponding LFCA data row (based on years available in cache)
            year_offset = max(0, endi - analysis_start_idx)
            
            if year_offset >= lfcs_data.shape[0]:
                break
                
            # Get LFCA time series up to current year
            month_end_index = (endi - analysis_start_idx + 1) * 12
            if month_end_index > lfcs_data.shape[1]:
                month_end_index = lfcs_data.shape[1]
            
            lfcs_current = lfcs_data[year_offset, :month_end_index]
            
            # Remove NaN values
            valid_indices = ~np.isnan(lfcs_current)
            if np.sum(valid_indices) < 24:  # Need at least 2 years of data
                continue
            
            lfcs_valid = lfcs_current[valid_indices]
            
            # Get corresponding temperature data (monthly)
            temp_monthly = np.repeat(temperature[analysis_start_idx:endi+1], 12)
            temp_monthly = temp_monthly[:month_end_index][valid_indices]
            
            if len(temp_monthly) != len(lfcs_valid):
                # Handle length mismatch
                min_len = min(len(temp_monthly), len(lfcs_valid))
                temp_monthly = temp_monthly[:min_len]
                lfcs_valid = lfcs_valid[:min_len]
            
            if len(lfcs_valid) < 24:  # Still need sufficient data
                continue
            
            try:
                # Perform regression
                slope, intercept, r_value, p_value, std_err = linregress(lfcs_valid, temp_monthly)
                
                # Calculate predicted temperatures
                predicted_temps = slope * lfcs_valid + intercept
                residuals = temp_monthly - predicted_temps
                residual_std_dev = np.std(residuals)
                
                # Calculate current estimate (blend of prediction and observation)
                if len(predicted_temps) >= 12:
                    recent_prediction = np.mean(predicted_temps[-12:])  # Last year average
                else:
                    recent_prediction = predicted_temps[-1]
                
                means_c[endi] = (recent_prediction + temperature[endi]) / 2 - preind_temp_base
                ses_c[endi] = residual_std_dev / np.sqrt(3)
                
                # Retrospective estimates (annual averages)
                if len(predicted_temps) >= 12:
                    pred_annual = predicted_temps.reshape(-1, 12).mean(axis=1)
                    retro_start = analysis_start_idx
                    retro_end = min(retro_start + len(pred_annual), len(means_r))
                    
                    for j, pred_val in enumerate(pred_annual):
                        if retro_start + j < retro_end:
                            means_r[retro_start + j] = (pred_val + temperature[retro_start + j]) / 2 - preind_temp_base
                            ses_r[retro_start + j] = residual_std_dev / np.sqrt(3)
                            
            except Exception as e:
                print(f"Regression failed for year {current_year}: {e}")
                continue

    return means_c, ses_c, means_r, ses_r
