import xarray as xr
import os
import numpy as np
import pandas as pd
import pdb


def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])

def convert_csv(ssp, cmip_model):
    if cmip_model == "ESM1-2-LR":
        cstr="/combined"
    elif cmip_model =="NorESM":
        cstr="_volc/BethkeEtAl2017"
    ts_file = os.path.expanduser(f"~/climate_data/{cmip_model}{cstr}/{ssp}_aod.nc")
    if cmip_model =="NorESM" and ssp == "hist":
        ts_file = os.path.expanduser(f"~/climate_data/{cmip_model}_volc/NorESM1-M-historical/{ssp}_aod.nc")
    with xr.open_dataset(ts_file) as ds:
        vararray = ds['__xarray_dataarray_variable__'].values

    return [average_every_n(vararray[i,:], 12) for i in range(np.shape(vararray)[0])]


def combine_write(ssp0, cmip_model,volc=""):
    if ssp0=="hist":
        sim_rec = np.array(convert_csv(ssp0, cmip_model))
        nens=3
        comb_recAOD = sim_rec.T/1000
    else:
        hist_rec = np.array(convert_csv("historical"+volc, cmip_model))
        sim_rec = np.array(convert_csv(ssp0, cmip_model))
        nens=20
        comb_recAOD = np.concatenate((hist_rec[0:nens,:].T, sim_rec.T))/1000

    erf_data = pd.read_csv(os.path.expanduser(f"~/climate_data/SSP_inputdata/ERFs-Smith-ar6/ERF_SSP245_1750-2500.csv"))
    
    
    
    if cmip_model == "ESM1-2-LR":
        years = np.arange(1850, 2101)  # Generate years from 1850 to 2100
        solar = erf_data['solar'][(1850-1750):(2101-1750)]
        contrails = erf_data['contrails'][(1850-1750):(2101-1750)]
    elif cmip_model =="NorESM" and ssp0!="hist":
        years = np.arange(1980, 2100)
        solar = erf_data['solar'][(1980-1750):(2100-1750)]
        contrails = erf_data['contrails'][(1980-1750):(2100-1750)]
    elif cmip_model =="NorESM" and ssp0=="hist":
        years = np.arange(1850, 2006)
        solar = erf_data['solar'][(1850-1750):(2006-1750)]
        contrails = erf_data['contrails'][(1850-1750):(2006-1750)]

    
    #df_solar = pd.read_csv(os.path.expanduser("~/Downloads/Thorne_15_codefigurestats/Common_Data/toyKFmodelData8c.csv"),header=None,index_col=0)
    #solar = df_solar[8]

    #aerosol_rad = erf_data['aerosol-radiation_interactions'][(1980-1750):(2100-1750)]
    #blackC = erf_data['bc_on_snow'][(1980-1750):(2100-1750)]
    #import statsmodels.api as sm
    #smoothed_aerosol_data = sm.nonparametric.lowess(aerosol_rad.values,range(1980,2100),frac=0.3)[:,1]

    comb_recAOD_cleaned = comb_recAOD + (np.tile(contrails,(nens,1)).T-0.015)/.18*.04  # get rid of gradual decline in the baseline over 21st century
    

    frac_blocked = 1-(5.5518/(comb_recAOD_cleaned+9.9735)) #(9.068/(comb_recAOD_cleaned+9.7279))
    #newly fitted values from excel spreadsheet
    solar_full = solar+ 340.4099428 - 0.108214
    tot_volc_erf = - np.tile(solar_full,(nens,1)).T*frac_blocked + 151.22 #23.7
    tot_nat_erf = tot_volc_erf + np.tile(solar,(nens,1)).T
    df = pd.DataFrame(tot_nat_erf, index=years)
    output_dir = os.path.expanduser("~/Downloads/Thorne_15_codefigurestats/Methods/43_Forcing_Based/3_Human_Induced/GWI_data")
    output_path2 = os.path.join(output_dir, f"ERFnatural_{cmip_model}_{ssp0}.csv")
    df.to_csv(output_path2, index=True)

if __name__ == "__main__":
    #combine_write("ssp126","ESM1-2-LR")
    #combine_write("ssp245","ESM1-2-LR")
    #combine_write("ssp370","ESM1-2-LR")
    #combine_write("rcp45Volc","NorESM","Volc")
    #combine_write("rcp45VolcConst_partial20","NorESM","Volc")
    combine_write("hist","NorESM")
