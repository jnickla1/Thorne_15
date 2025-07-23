import xarray as xr
import config
import numpy as np
import pandas as pd
import pdb
import os

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])

def convert_csv(ssp, cmip_model):
    if cmip_model == "ESM1-2-LR":
        cstr="/combined"
    elif cmip_model =="NorESM":
        cstr="_volc/BethkeEtAl2017"
    ts_file = config.CLIMATE_DATA_PATH+f"/{cmip_model}{cstr}/{ssp}_aave_tas.nc"
    if cmip_model =="NorESM" and ssp == "hist":
        ts_file = config.CLIMATE_DATA_PATH+f"/{cmip_model}_volc/NorESM1-M-historical/{ssp}_aave_tas.nc"
    with xr.open_dataset(ts_file) as ds:
        vararray = ds['tas'].values

    return [average_every_n(vararray[i,:], 12) for i in range(np.shape(vararray)[0])]


def combine_write(ssp0, cmip_model,volc=""):
    if ssp0=="hist":
        comb_rec = np.array(convert_csv(ssp0, cmip_model)).T
    else:
        hist_rec = np.array(convert_csv("historical"+volc, cmip_model))
        sim_rec = np.array(convert_csv(ssp0, cmip_model))
        comb_rec = np.concatenate((hist_rec.T, sim_rec.T))
    if cmip_model == "ESM1-2-LR":
        years = np.arange(1850, 2101)  # Generate years from 1850 to 2100
    elif cmip_model =="NorESM" and ssp0!="hist":
        years = np.arange(1980, 2100)
    elif cmip_model =="NorESM" and ssp0=="hist":
        years = np.arange(1850, 2006)

    df = pd.DataFrame(comb_rec, index=years)
    output_dir = config.CODEBASE_PATH+"/Methods/43_Forcing_Based/3_Human_Induced/GWI_data"
    output_path2 = os.path.join(output_dir, f"ts_{cmip_model}_{ssp0}.csv")
    df.to_csv(output_path2, index=True)

if __name__ == "__main__":
    #combine_write("ssp126","ESM1-2-LR")
    #combine_write("ssp245","ESM1-2-LR")
    #combine_write("ssp370","ESM1-2-LR")
    #combine_write("rcp45Volc","NorESM","Volc")
    #combine_write("rcp45VolcConst","NorESM","Volc")
    combine_write("hist","NorESM","Volc")
