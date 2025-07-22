import numpy as np
import pandas as pd
import config
#EBMKF_Nicklas as ekf
from netCDF4 import Dataset
import pdb
import xarray as xr
from sklearn.linear_model import LinearRegression
from gen_eCO2 import calculate_equivalent_co2

def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])


def ta_smooth(orig_opt_depths,fill_value,optical=True):
    if(optical):
        wt_opt_depths = 1/(orig_opt_depths+9.7279)
    else:
        wt_opt_depths =orig_opt_depths
    N = 30
    nwt_opt_depths=np.full(len(orig_opt_depths),fill_value)
    cN=int(np.ceil(N/2))
    fN=int(np.floor(N/2))
    for i in range((fN),(len(nwt_opt_depths)-1)):
        lasta=i+cN;firsta=i-fN
        nwt_opt_depths[i] = (np.sum(wt_opt_depths[(firsta):i+1])+ fill_value*(cN-1))/N
    #computing half-average - future is assumed to be the average
    if(optical):
        return(1/nwt_opt_depths-9.7279)
    else:
        return(nwt_opt_depths)


def ext_timeser(timeseries, start_year=1850, end_year=2014, 
                                       fit_start=1980, project_end=2023, random_seed=320):
    """
    Extend a historical timeseries to project 10 more years using a linear trend + Gaussian noise.

    Parameters:
        timeseries (array-like or pd.Series): Timeseries data from start_year to end_year.
        start_year (int): First year of the input timeseries (default 1850).
        end_year (int): Last year of the input timeseries (default 2014).
        fit_start (int): Start year of the trend fit (default 1980).
        project_end (int): Final year of the projection (default 2024).
        random_seed (int or None): Seed for reproducibility.

    Returns:
        pd.Series: Extended timeseries from start_year to project_end.
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    # Ensure Series with year index
    if not isinstance(timeseries, pd.Series):
        years = np.arange(start_year, end_year + 1)
        timeseries = pd.Series(timeseries, index=years)

    # Subset for regression
    fit_years = np.arange(fit_start, end_year + 1).reshape(-1, 1)
    fit_values = timeseries.loc[fit_start:end_year].values

    # Fit linear regression
    model = LinearRegression()
    model.fit(fit_years, fit_values)
    trend = model.predict(fit_years)

    # Compute RMSE
    rmse = np.sqrt(np.mean((fit_values - trend) ** 2))

    # Project future values
    future_years = np.arange(end_year + 1, project_end + 1)
    future_X = future_years.reshape(-1, 1)
    future_trend = model.predict(future_X)
    future_noise = np.random.normal(loc=0, scale=rmse, size=len(future_years))
    future_values = future_trend + future_noise

    # Combine original and projected data
    extended_series = pd.concat([
        timeseries,
        pd.Series(future_values, index=future_years)
    ])
    return extended_series




#temperature has 2024, dont have all forcings for this yet - must update ekf.n_iters
data_orig = pd.read_csv("../../../Common_Data/HadCRUT5.csv")
temps_obs = data_orig.loc[:,"Anomaly"].to_numpy()
preind_base = np.mean(temps_obs[0:50])


#ekf.R_tvar=np.square(temps_1std[0:ekf.n_iters])
#ekf.fdbkA = 0.35 #original paper value for this estimate
#ekf.precompute_coeffs(False)
##empser  = np.full(np.shape(years),np.nan)
##means = empser.copy()
##ses = empser.copy()
##means2 = empser.copy()
##ses2 = empser.copy()
import EBMKF_Nicklas4 as ekf #only change cloud feedback, but dynamically


###FIRST RUN WITHOUT TA_SMOOTH is the total
total_warming = ekf.xh1s.copy()
total_warming_se  = ekf.stdP.copy()

#if experiment_type == "historical":
ekf.n_iters = 174
unf_new_opt_depth = ekf.data[:,3]*0.001
ekf.opt_depth=ta_smooth(unf_new_opt_depth,ekf.involcavg)
given_preind_base = np.mean(temps_obs[0:50])
temps = temps_obs-given_preind_base+preind_base + ekf.offset #ensure this matches expected starting temperature in 1850
frac_blocked = 1-(5.5518/(unf_new_opt_depth+9.9735)) #(9.068/(comb_recAOD_cleaned+9.7279))
erf_data_solar = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")['solar'].to_numpy()[(1850-1750):]
solar_full = erf_data_solar[:ekf.n_iters]+ 340.4099428 - 0.108214
tot_volc_erf = (- solar_full*frac_blocked + 151.22)
erf_trapez=(1*tot_volc_erf[0:-2]+2*tot_volc_erf[1:-1]+0*tot_volc_erf[2:] )/4 #want to snap back as quickly as possible, no volc signal
temps[2:-1] = temps[2:-1] - erf_trapez/(ekf.heatCp - ekf.Cs)
ohca_meas = ekf.oc_meas*ekf.zJ_from_W
#temps[7:] = temps[7:] + erf_trapez[:-5]/(ekf.heatCp - ekf.Cs)
ohca_meas[2:] = ohca_meas[2:] - (erf_trapez/(ekf.Cs + ekf.Cd))*ekf.zJ_from_W*200 #not perfect but good for now, also unclear where the 200 comes from
TOA_meas_artif1 =  ekf.TOA_meas_artif - (tot_volc_erf ) * .75 
new_observ = np.array([[temps[:ekf.n_iters],ohca_meas/ekf.zJ_from_W ,TOA_meas_artif1]]).T
ekf.dfrA_float_var = ekf.dfrA_float_var/40

###SECOND RUN WITH TA_SMOOTH is the anthro
warming,warming_se,_,_ = ekf.ekf_run(new_observ,ekf.n_iters,retPs=3)
anthro_warming = warming.copy()
anthro_warming_se = warming_se.copy()


###THRID RUN IS hist-aer
#remember to extend to future linearly with huge uncertainty to cover 2014-2024
#solar set to const
orig_tsi = ekf.tsi.copy()
ekf.tsi = np.full(ekf.n_iters, ekf.sw_in)
#volc set to const
ekf.opt_depth=np.full(ekf.n_iters,(1/ekf.involcavg-9.7279))
#aer varies - can keep the same
#co2 set to const preind
orig_lCo2 = ekf.lCo2.copy()
ekf.lCo2= np.full(ekf.n_iters,np.log10(ekf.data[0,2]))
#new_obs - be sure to extend to 2024. Uncertanties unchanged
aer_temps_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/ts/aave_ts_Amon_MPI-ESM1-2-LR_hist-aer_r4i1p1f1_gn_185001-201412.nc", 'r').variables['ts']
aer_temps = average_every_n(aer_temps_raw[:].__array__(), 12)
given_preind_base = np.mean(aer_temps[0:50])
aer_temps = aer_temps - given_preind_base + preind_base + ekf.offset
aer_ohca_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/ohca_aave/ohca_thetaot_Emon_MPI-ESM1-2-LR_hist-aer_r4i1p1f1_gn_185001-201412.nc", 'r').variables['__xarray_dataarray_variable__']
aer_ohca = average_every_n(aer_ohca_raw[:].__array__(), 12)
aer_ohca = aer_ohca - aer_ohca[0]
aer_toa_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/rt/combined_global_rt/hist-aer/global_rt_Amon_MPI-ESM1-2-LR_hist-aer_r4i1p1f1_gn_185001-186912.nc-201412.nc", 'r').variables['rt']
aer_toa = average_every_n(aer_toa_raw[:].__array__(), 12)
new_observ_aer = np.array([[ext_timeser(aer_temps),ext_timeser(aer_ohca/ekf.zJ_from_W) ,ext_timeser(aer_toa)]]).T
warming,warming_se,_,_ = ekf.ekf_run(new_observ_aer,ekf.n_iters,retPs=3)
aer_warming = warming.copy()
aer_warming_se = warming_se.copy()

###FOURTH RUN is hist-totalO3
#counts as a type of CO2
#recalculate CO2
forcingsdf = pd.read_csv(config.CLIMATE_DATA_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")
forcing_O3= forcingsdf['solar'].to_numpy()[(1850-1750):(2024-1750)]
ekf.lCo2= np.log10(calculate_equivalent_co2(forcing_O3))
#aer set to 0
orig_anthro_clouds = (ekf.data[:,7]+1)
ekf.anthro_clouds = np.full(ekf.n_iters, orig_anthro_clouds[0])
#solar set to const - don't change
#volc set to const - don't change
totalO3_temps_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/ts/aave_ts_Amon_MPI-ESM1-2-LR_hist-totalO3_r4i1p1f1_gn_185001-201412.nc", 'r').variables['ts']
totalO3_temps = average_every_n(totalO3_temps_raw[:].__array__(), 12)
given_preind_base = np.mean(totalO3_temps[0:50])
totalO3_temps = totalO3_temps - given_preind_base + preind_base + ekf.offset
totalO3_ohca_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/ohca_aave/ohca_thetaot_Emon_MPI-ESM1-2-LR_hist-totalO3_r4i1p1f1_gn_185001-201412.nc", 'r').variables['__xarray_dataarray_variable__']
totalO3_ohca = average_every_n(totalO3_ohca_raw[:].__array__(), 12)
totalO3_ohca = totalO3_ohca - totalO3_ohca[0]
totalO3_toa_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/rt/combined_global_rt/hist-totalO3/global_rt_Amon_MPI-ESM1-2-LR_hist-totalO3_r4i1p1f1_gn_185001-186912.nc-201412.nc", 'r').variables['rt']
totalO3_toa = average_every_n(totalO3_toa_raw[:].__array__(), 12)
new_observ_totalO3 = np.array([[ext_timeser(totalO3_temps),ext_timeser(totalO3_ohca/ekf.zJ_from_W) ,ext_timeser(totalO3_toa)]]).T
warming,warming_se,_,_= ekf.ekf_run(new_observ_totalO3,ekf.n_iters,retPs=3)
totalO3_warming = warming.copy()
totalO3_warming_se = warming_se.copy()

#averaging hist-aer and hist-totalO3 to write OHF

ohf_warming = np.mean(np.array([totalO3_warming, aer_warming]),axis=0) #check axis
ohf_warming_se = np.sqrt( np.mean(np.array([totalO3_warming_se**2, aer_warming_se**2]),axis=0))

###FIFTH RUN IS hist-GHG
#recalculate CO2 with only GHGs
forcing_GHG_allyrs= forcingsdf['co2'].to_numpy() + forcingsdf['ch4'].to_numpy() + forcingsdf['n2o'].to_numpy() + forcingsdf['other_wmghg'].to_numpy()
ekf.lCo2= np.log10(calculate_equivalent_co2(forcing_GHG_allyrs[(1850-1750):(2024-1750)]))
#solar set to const - don't change
#volc set to const - don't change
#aer set to 0 - dont change
GHG_temps_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/ts/aave_ts_Amon_MPI-ESM1-2-LR_hist-GHG_r4i1p1f1_gn_185001-201412.nc", 'r').variables['ts']
GHG_temps = average_every_n(GHG_temps_raw[:].__array__(), 12)
given_preind_base = np.mean(GHG_temps[0:50])
GHG_temps = GHG_temps - given_preind_base + preind_base+ ekf.offset
#GHG_temps - ekf.Teq1850 only about 1.1 for run 4, close to 0.95 of HadCrut
GHG_ohca_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/ohca_aave/ohca_thetaot_Emon_MPI-ESM1-2-LR_hist-GHG_r4i1p1f1_gn_185001-201412.nc", 'r').variables['__xarray_dataarray_variable__']
GHG_ohca = average_every_n(GHG_ohca_raw[:].__array__(), 12)
GHG_ohca = GHG_ohca - GHG_ohca[0]
GHG_toa_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/rt/combined_global_rt/hist-GHG/global_rt_Amon_MPI-ESM1-2-LR_hist-GHG_r4i1p1f1_gn_185001-186912.nc-201412.nc", 'r').variables['rt']
GHG_toa = average_every_n(GHG_toa_raw[:].__array__(), 12)
new_observ_GHG = np.array([[ext_timeser(GHG_temps),ext_timeser(GHG_ohca/ekf.zJ_from_W) ,ext_timeser(GHG_toa)]]).T
warming,warming_se,_,_ = ekf.ekf_run(new_observ_GHG,ekf.n_iters,retPs=3)
GHG_warming = warming.copy()
GHG_warming_se = warming_se.copy()
breakpoint()

###SIXTH RUN IS hist-nat
#solar set to original
ekf.tsi  = orig_tsi 
#volc set to original
ekf.opt_depth = unf_new_opt_depth
#co2 set to const preind
ekf.lCo2= np.full(ekf.n_iters,np.log10(ekf.data[0,2]))
#aer set to 0 - dont change
nat_temps_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/ts/aave_ts_Amon_MPI-ESM1-2-LR_hist-nat_r4i1p1f1_gn_185001-201412.nc", 'r').variables['ts']
nat_temps = average_every_n(nat_temps_raw[:].__array__(), 12)
given_preind_base = np.mean(nat_temps[0:50])
nat_temps = nat_temps - given_preind_base + preind_base + ekf.offset
nat_ohca_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/ohca_aave/ohca_thetaot_Emon_MPI-ESM1-2-LR_hist-nat_r4i1p1f1_gn_185001-201412.nc", 'r').variables['__xarray_dataarray_variable__']
nat_ohca = average_every_n(nat_ohca_raw[:].__array__(), 12)
nat_ohca = nat_ohca - nat_ohca[0]
nat_toa_raw = Dataset(config.CLIMATE_DATA_PATH+"/attributable-warming_ESM1-2-LR/rt/combined_global_rt/hist-nat/global_rt_Amon_MPI-ESM1-2-LR_hist-nat_r4i1p1f1_gn_185001-186912.nc-201412.nc", 'r').variables['rt']
nat_toa = average_every_n(nat_toa_raw[:].__array__(), 12)
new_observ_nat = np.array([[ext_timeser(nat_temps),ext_timeser(nat_ohca/ekf.zJ_from_W) ,ext_timeser(nat_toa)]]).T
warming,warming_se,_,_ = ekf.ekf_run(new_observ_nat,ekf.n_iters,retPs=3)
nat_warming = warming.copy()
nat_warming_se = warming_se.copy()

ekf.plt.close()



from scipy.stats import norm

def save_analytic_percentiles_from_mean_se(output_csv, year_range, data_dict, percentiles=[5, 17, 83, 95, 50]):
    """
    For each variable in data_dict, compute specified percentiles from N(mean, se) using z-scores.
    
    Parameters:
        output_csv (str): Filename to write output CSV.
        year_range (array-like): List of years, e.g., np.arange(1850, 2025).
        data_dict (dict): Keys are variable names, values are tuples (mean_vector, se_vector).
        percentiles (list): List of percentiles to compute (default: [5, 17, 50, 83, 95]).
    """
    df = pd.DataFrame({'year': year_range})
    
    z_scores = {p: norm.ppf(p / 100) for p in percentiles}

    for varname, (mean_vec, se_vec) in data_dict.items():
        for p in percentiles:
            z = z_scores[p]
            df[f'{varname}_p{p}'] = mean_vec + z * se_vec - ekf.Teq1850

    df.to_csv(output_csv, index=False)


years = np.arange(1850, 2024)
n = len(years)


data = {
    'aer': (aer_warming, aer_warming_se),
    'ghg': (GHG_warming, GHG_warming_se),
    'ohf': (ohf_warming, ohf_warming_se),
    'ant': (anthro_warming, anthro_warming_se),
    'nat': (nat_warming, nat_warming_se),
    'tot': (total_warming, total_warming_se)
}

save_analytic_percentiles_from_mean_se("Nicklas_GMST_timeseries.csv", years, data)
