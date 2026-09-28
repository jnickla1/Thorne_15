#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os

import sys
from pathlib import Path

# add parent directory to PYTHONPATH
sys.path.append(str(Path(__file__).resolve().parents[1]))

import EBMKF_Nicklas4 as ekf

# ---- you fill these in ----
CLIM_PATH = os.path.expanduser('~/')+"climate_data"
model_run = 0

colors = {"Volc": "orangered", "VolcConst": "dodgerblue"}
colors2 = {"Volc": "maroon", "VolcConst": "navy"}

def average_every_n(lst, n):
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])

def gen_orig_number(new_member_number, sz_ens):
    nums = np.arange(1, sz_ens + 1)
    reshuffled = sorted([f"{n}|" for n in nums])
    recovered_order = [int(s.rstrip("|")) for s in reshuffled]
    return recovered_order if new_member_number == -1 else recovered_order[new_member_number]

def load_tas_annual( rcp: str, volc_variant: str, model_run: int,
                         temps_obs_past: np.ndarray):
    """
    NorESM temperature series construction matching your select_data(..., NorESM) branch.

    Parameters
    ----------
    rcp : str
        e.g. "rcp45"
    volc_variant : str
        "Volc" or "VolcConst" or "NoVolc" (you'll ignore NoVolc later)
    model_run : int
        ensemble member index (0-based)
    temps_obs_past : np.ndarray
        HadCRUT5 anomalies with 1850-1899 baseline removed (like your code)

    Returns
    -------
    years_all : np.ndarray
        The concatenated year axis built exactly like your code.
    simall : np.ndarray
        The concatenated annual temperature anomaly series ("simall").
    """
    cdataprefix = CLIM_PATH + '/'
    # --- file locations exactly like collect_data(..., NorESM) ---
    fut_data_loc  = f"{cdataprefix}NorESM_volc/BethkeEtAl2017/{rcp}{volc_variant}_aave_tas.nc"
    hist_data_loc = f"{cdataprefix}NorESM_volc/BethkeEtAl2017/historicalVolc_aave_tas.nc"  # for Volc + VolcConst

    # --- read FUT (monthly) + time ---
    with Dataset(fut_data_loc, "r") as ds:
        sims_tas_row = ds.variables["tas"][model_run, :].__array__()
        stime_mon = ds.variables["time"][:].__array__() / 365.0 + 1850.0

    # --- read HIST (monthly) + time ---
    with Dataset(hist_data_loc, "r") as ds:
        sims_tas_hist_row = ds.variables["tas"][model_run, :].__array__()
        stime_mon_hist = ds.variables["time"][:].__array__() / 365.0 + 1850.0

    # monthly -> annual, exactly like select_data
    this_sim_yr  = average_every_n(sims_tas_row, 12)
    this_hsim_yr = average_every_n(sims_tas_hist_row, 12)

    stime_yrs  = np.floor(average_every_n(stime_mon, 12)).astype(int)
    shtime_yrs = np.floor(average_every_n(stime_mon_hist, 12)).astype(int)

    # ---- long past replacement, exactly like select_data ----
    # NOTE: your original uses np.shape(sims_tas)[0] for sz_ens, which for NorESM is 60.
    # Use 60 explicitly to match the intended grouping of NorESM1-M historical members.
    long_past_index = ((gen_orig_number(model_run, 60) - 1) // 20)

    long_past_data_loc = f"{cdataprefix}NorESM_volc/NorESM1-M-historical/hist_aave_tas.nc"
    with Dataset(long_past_data_loc, "r") as ds:
        long_past_row = ds.variables["tas"][long_past_index, :].__array__()
    long_past_tas = average_every_n(long_past_row, 12)

    # ---- offset_sim, exactly like your code ----
    # baseline 1980-2000 match (20 years)
    offset_sim = np.mean(this_hsim_yr[0:20] - temps_obs_past[(-1850 + 1980):(-1850 + 1980 + 20)])

    sim_corrected = this_sim_yr - offset_sim
    simh_corrected = this_hsim_yr - offset_sim
    long_past_tas_corrected = long_past_tas - offset_sim

    # ---- concatenation exactly like your code ----
    simall = np.concatenate((
        long_past_tas_corrected[0:(-1850 + 1980)],
        simh_corrected,
        sim_corrected
    ))

    fixyrs = 1980 - 3829
    years_all = np.arange(1850,len(simall)+1850)

    return years_all, simall


def load_all_inputs(rcp: str, volc_variant: str, new_iter: int, ekf, model_run: int = 0):
    
    """
    NorESM inputs loader, mirroring EBMKF_ta4_method.py (NorESM branch).

    Parameters
    ----------
    rcp : str
        e.g. "rcp45"
    volc_variant : str
        "Volc" or "VolcConst" (ignore "NoVolc" per your use-case)
    new_iter : int
        length of annual axis (typically 2101-1850 = 251)
    ekf : module/object
        Your EBMKF_Nicklas4 module (needed for involcavg / baseline constant)
    model_run : int
        ensemble member index (0-based)
    CLIM_PATH : str
        e.g. os.path.expanduser("~/") + "climate_data"
    GWI_DATA_PATH : str
        path to the GWI_data folder that contains the ERF_*_full.csv files

    Returns
    -------
    unf_new_opt_depth : (new_iter,) float
    ohca_meas         : (new_iter,) float
    eco2              : (N,) float   (equivalent CO2 series; length depends on CSV)
    new_rtTOA         : (new_iter,) float
    """
    if CLIM_PATH is None:
        raise ValueError("Please pass CLIM_PATH (e.g. ~/climate_data).")
    GWI_DATA_PATH = str(Path(__file__).resolve().parents[2]) + "/3_Human_Induced/GWI_data/"
    print(GWI_DATA_PATH)

    if volc_variant not in ["Volc", "VolcConst"]:
        raise ValueError("volc_variant must be 'Volc' or 'VolcConst'")

    # ---- AOD/optical depth ----
    unf_new_opt_depth = np.full(new_iter, (1/ekf.involcavg - 9.7279), dtype=np.float64)

    long_past_index = ((gen_orig_number(model_run, 60) - 1) // 20)

    # spinup: 1850-2005 (monthly -> annual)
    aod_spinup = Dataset(f"{CLIM_PATH}/NorESM_volc/NorESM1-M-historical/hist_aod.nc", "r") \
        .variables["__xarray_dataarray_variable__"][:].__array__()
    unf_new_opt_depth[(1850-1850):(2006-1850)] = average_every_n(aod_spinup[long_past_index, :], 12) / 1000.0

    # historicalVolc overwrites: 1980-2005 (monthly -> annual)
    aod_earlier = Dataset(f"{CLIM_PATH}/NorESM_volc/BethkeEtAl2017/historicalVolc_aod.nc", "r") \
        .variables["__xarray_dataarray_variable__"][:].__array__()
    unf_new_opt_depth[(1980-1850):(2006-1850)] = average_every_n(aod_earlier[model_run, :], 12) / 1000.0

    # future: 2006-2099 (monthly -> annual)
    if volc_variant == "Volc":
        aod_later = Dataset(f"{CLIM_PATH}/NorESM_volc/BethkeEtAl2017/{rcp}Volc_aod.nc", "r") \
            .variables["__xarray_dataarray_variable__"][:].__array__()
        unf_new_opt_depth[(2006-1850):(2100-1850)] = average_every_n(aod_later[model_run, :], 12) / 1000.0
    else:  # VolcConst
        aod_later = Dataset(f"{CLIM_PATH}/NorESM_volc/BethkeEtAl2017/{rcp}VolcConst_partial20_aod.nc", "r") \
            .variables["__xarray_dataarray_variable__"][:].__array__()
        unf_new_opt_depth[(2006-1850):(2100-1850)] = average_every_n(aod_later[model_run % 20, :], 12) / 1000.0

    # ---- OHCA ----
    ohca_meas = np.zeros(new_iter, dtype=np.float64)

    # spinup: 1850-2005 (monthly -> annual), rebased to start at 0
    ohca_spinup = Dataset(f"{CLIM_PATH}/NorESM_volc/NorESM1-M-historical/hist_ohca_mon.nc", "r") \
        .variables["__xarray_dataarray_variable__"][:].__array__()
    ohca_s_yr = average_every_n(ohca_spinup[long_past_index, :], 12)
    ohca_meas[(1850-1850):(2006-1850)] = ohca_s_yr - ohca_s_yr[0]

    # historicalVolc OHCA: 1980-2005 (ANNUAL already in your method)
    ohca_earlier = Dataset(f"{CLIM_PATH}/NorESM_volc/OHCA/historicalVolc_ohca.nc", "r") \
        .variables["__xarray_dataarray_variable__"][:].__array__()
    # "missing one model_run ohca, so shift them"
    mr_ohca = (model_run if model_run < 37 else model_run - 1)
    ohca_meas[(1980-1850):(2006-1850)] = (
        ohca_earlier[mr_ohca, :]
        + 2*ohca_meas[(1979-1850)] - ohca_meas[(1978-1850)]
    )

    # future OHCA: 2006-2099 (ANNUAL)
    if volc_variant == "Volc":
        ohca_later = Dataset(f"{CLIM_PATH}/NorESM_volc/OHCA/{rcp}Volc_ohca.nc", "r") \
            .variables["__xarray_dataarray_variable__"][:].__array__()
        ohca_meas[(2006-1850):(2100-1850)] = (
            ohca_later[model_run, :]
            + 2*ohca_meas[(2005-1850)] - ohca_meas[(2004-1850)]
        )
    else:
        ohca_later = Dataset(f"{CLIM_PATH}/NorESM_volc/OHCA/{rcp}VolcConst_partial20_ohca.nc", "r") \
            .variables["__xarray_dataarray_variable__"][:].__array__()
        idx = (model_run if (model_run > 2 and model_run < 14) else model_run % 14)
        ohca_meas[(2006-1850):(2100-1850)] = (
            ohca_later[idx, :]
            + 2*ohca_meas[(2005-1850)] - ohca_meas[(2004-1850)]
        )

    # ---- rtmt -> new_rtTOA ----
    new_rtTOA = np.zeros(new_iter, dtype=np.float64)

    # spinup: 1850-2005 (monthly -> annual)
    rt_spinup = Dataset(f"{CLIM_PATH}/NorESM_volc/NorESM1-M-historical/rtmt/rtmt_NorESM1-M_historical.nc", "r") \
        .variables["rtmt"][:].__array__()
    new_rtTOA[(1850-1850):(2006-1850)] = average_every_n(rt_spinup[long_past_index, :], 12)

    # historicalVolc overwrite: 1980-2005 (monthly -> annual)
    rt_earlier = Dataset(f"{CLIM_PATH}/NorESM_volc/BethkeEtAl2017/historicalVolc_rtmt/rtmt_NorESM1-M_historicalVolc.nc", "r") \
        .variables["rtmt"][:].__array__()
    new_rtTOA[(1980-1850):(2006-1850)] = average_every_n(rt_earlier[model_run, :], 12)

    # future: 2006-2099 (monthly -> annual)
    if volc_variant == "Volc":
        rt_later = Dataset(f"{CLIM_PATH}/NorESM_volc/BethkeEtAl2017/{rcp}Volc_rtmt/rtmt_NorESM1-M_{rcp}Volc.nc", "r") \
            .variables["rtmt"][:].__array__()
        new_rtTOA[(2006-1850):(2100-1850)] = average_every_n(rt_later[model_run, :], 12)
    else:
        rt_later = Dataset(f"{CLIM_PATH}/NorESM_volc/BethkeEtAl2017/{rcp}VolcConst_rtmt/rtmt_NorESM1-M_{rcp}VolcConst_20.nc", "r") \
            .variables["rtmt"][:].__array__()
        new_rtTOA[(2006-1850):(2100-1850)] = average_every_n(rt_later[model_run % 20, :], 12)

    # ---- CO2: equivalent CO2 from ERF_wmghg via gen_eCO2 (same as EKF NorESM branch) ----
    # NOTE: for Volc, your method uses ERFanthro_NorESM_rcp45Volc_full.csv (not the SSP_inputdata CSV).
    #       for VolcConst, it uses ERF_NorESM_rcp45VolcConst_full.csv.
    # We return the "model_outputlCo2" (equivalent CO2), which is the most consistent thing to plot.
    from gen_eCO2 import calculate_equivalent_co2  # same module used in your method file

    if volc_variant == "Volc":
        erf_data = pd.read_csv(f"{GWI_DATA_PATH}/ERFanthro_NorESM_{rcp}Volc_full.csv")
    else:
        erf_data = pd.read_csv(f"{GWI_DATA_PATH}/ERF_NorESM_{rcp}VolcConst_full.csv")

    eco2 = calculate_equivalent_co2(erf_data["ERF_wmghg"].values)

    # ---- contrails tweak to unf_new_opt_depth (copied behavior) ----
    # Your method uses contrails from the same erf_data for NorESM and adds this after 1980.
    if "contrails" in erf_data.columns:
        contrails = erf_data["contrails"][(1980-1750):(2100-1750)].values
        unf_new_opt_depth[(1980-1850):] = (
            unf_new_opt_depth[(1980-1850):] + (contrails - 0.015)/0.18*0.04
        )

    return unf_new_opt_depth, ohca_meas, eco2, new_rtTOA

def main(ekf):
    # Define x-axes (edit to match your actual time bases)
    # optical depth is annual-like length new_iter (years)
    # OHCA/RT are annual after average_every_n
    # CO2 length depends on your CSV
    years_annual = np.arange(1850, 2100)  # <-- replace with your actual years array

    fig, axes = plt.subplots(5, 1, figsize=(7, 10), sharex=True)
    axT, axOH, axCO2, axAOD, axRF = axes
    fig2, ax2= plt.subplots(figsize=(9,4))
    erf_data = erf_data = pd.read_csv(CLIM_PATH+"/SSP_inputdata/ERFs-Smith-ar6/ERF_ssp245_1750-2500.csv")
    contrails = erf_data['contrails'][(1980-1750):(2100-1750)].values

    # ----- OBSERVED series (black): wire these in -----
    # obs_temp, obs_ohca, obs_co2, obs_aod, obs_rf should be aligned to same x-axis you plot
    # axT.plot(years_annual, obs_temp, color="black")
    # etc.
    hist_record = pd.read_csv("../../../../Common_Data/HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv")
    HadCRUT_ts =  hist_record["Anomaly (deg C)"].values
    HadCRUT_ts = HadCRUT_ts - np.mean(HadCRUT_ts[0:50])
    for case in ["Volc", "VolcConst"]:
        unf_aod, ohca, eco2, rf = load_all_inputs("rcp45",case, len(years_annual), ekf,0)

        years, tas = load_tas_annual("rcp45",case, 0,HadCRUT_ts)

        
        axT.plot(years_annual, tas, color=colors[case], label=case)

        axOH.plot(years_annual[:len(ohca)], ohca, color=colors[case])
        axCO2.plot(years_annual[:len(eco2)], eco2, color=colors[case])
        axAOD.plot(years_annual[:len(unf_aod)], unf_aod, color=colors[case])
        ax2.plot(years_annual[:len(unf_aod)], unf_aod, color=colors[case],label = case + " corr0")
        unf_aod2 = unf_aod.copy()
        unf_aod2[(1980-1850):] = unf_aod2[(1980-1850):] + (contrails-0.015)/.18*.04
        ax2.plot(years_annual[:len(unf_aod)], unf_aod2, color=colors2[case],label= case + " corr1")
        axRF.plot(years_annual[:len(rf)], rf, color=colors[case])

    axT.plot(years_annual[:174], ekf.temps - np.mean(ekf.temps[0:50]), color='k',label="observed")
    axOH.plot(years_annual[:174], ekf.ocean_heat_measured, color='k')
    meas_CO2 = pd.read_csv("../../../../Common_Data/co2_etheridge_keeling.csv")
    axCO2.plot(years_annual[:174], meas_CO2['mean'][:174], color='k')
    axAOD.plot(years_annual[:174], ekf.data[:, 3] * 0.001, color='k')
    ax2.plot(years_annual[:174], ekf.data[:, 3] * 0.001, color='k',ls='--',label = "observed")
    ax2.axhline(0,color='grey', ls=':', lw=1)
    ax2.set_ylabel("AOD (unitless)")
    ax2.set_xlabel("year")
    ax2.set_title("Further Correction to NorESM AOD")
    ax2.legend(ncol=3)
    ax2.set_xlim([1980,2100])
    
    axRF.plot(years_annual[(2001-1850):173], ekf.TOA_crop, color='k')

    
    
    axT.set_title("Temperature (relative to 1850-1900)")
    axT.set_ylabel("°C")
    axOH.set_title("Ocean Heat Content (OHCA)")
    axOH.set_ylabel("ZJ")
    axCO2.set_title("CO₂")
    axCO2.set_ylabel("ppm")
    axAOD.set_title("Optical depth (unfiltered)")
    axAOD.set_ylabel("unitless")
    axRF.set_title("Net radiative forcing (annual)")
    axRF.set_ylabel("W/m2")

    axT.legend(loc="upper left", frameon=False)
    axes[-1].set_xlabel("Year")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    fig.suptitle("NorESM (Volc) Uncorrected Inputs")
    plt.tight_layout()
    plt.show()

    

    
    #fig.savefig('NorESM_inputs.png', dpi=600, bbox_inches='tight')
    fig2.savefig('NorESM_AO2_contrail_correction.png', dpi=600, bbox_inches='tight')

# You’ll call main(ekf) from wherever ekf is available in your project.

if __name__ == "__main__":
    # Example:
    # plot_all_inputs("mpi_inputs_panels.png")
    # ---- initialize EKF ----

    # ---- plot MPI inputs ----
    main(ekf)
