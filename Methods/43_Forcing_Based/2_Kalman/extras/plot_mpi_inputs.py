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

colors = {"ssp126": "green", "ssp245": "purple", "ssp370": "red"}

def average_every_n(arr, n):
    return np.array([np.mean(arr[i:i+n]) for i in range(0, len(arr), n)])

def load_mpi_tas_annual(cdataprefix, ssp: str, model_run: int):
    """
    Returns:
      years (int array)
      tas_annual_anom (float array): concatenated historical+future, baseline 1850-1899 = 0
    """
    fut_path  = f"{cdataprefix}ESM1-2-LR/combined/{ssp.lower()}_aave_tas.nc"
    hist_path = f"{cdataprefix}ESM1-2-LR/combined/historical_aave_tas.nc"

    # --- future ---
    with Dataset(fut_path, "r") as ds:
        tas_f_m = np.asarray(ds.variables["tas"][model_run, :], dtype=np.float64)
        t_f = np.asarray(ds.variables["time"][:], dtype=np.float64) / 365.0 + 1850.0

    # --- historical ---
    with Dataset(hist_path, "r") as ds:
        tas_h_m = np.asarray(ds.variables["tas"][model_run, :], dtype=np.float64)
        t_h = np.asarray(ds.variables["time"][:], dtype=np.float64) / 365.0 + 1850.0

    # monthly -> annual
    tas_f_y = average_every_n(tas_f_m, 12)
    tas_h_y = average_every_n(tas_h_m, 12)

    years_f = np.floor(average_every_n(t_f, 12)).astype(int)
    years_h = np.floor(average_every_n(t_h, 12)).astype(int)

    # baseline 1850–1899 using *historical* first 50 annual values
    baseline = tas_h_y[:50].mean()
    tas_all = np.concatenate([tas_h_y, tas_f_y]) - baseline
    years_all = np.concatenate([years_h, years_f])

    return years_all, tas_all


def load_mpi_esm1_2_lr_inputs(ssp: str, new_iter: int, ekf):
    # Optical depth (unfiltered)
    unf_new_opt_depth = np.full(new_iter, (1/ekf.involcavg - 9.7279))
    unf_new_opt_depth[:ekf.n_iters] = ekf.data[:, 3] * 0.001

    # OHCA (annual)
    ohca_later = Dataset(f"{CLIM_PATH}/ESM1-2-LR/opottempmint/{ssp}_ohca.nc", "r") \
        .variables["__xarray_dataarray_variable__"][:].__array__()
    ohca_earlier = Dataset(f"{CLIM_PATH}/ESM1-2-LR/opottempmint/historical_ohca.nc", "r") \
        .variables["__xarray_dataarray_variable__"][:].__array__()

    ohca_ly = average_every_n(ohca_later[model_run, :], 12)
    ohca_ey = average_every_n(ohca_earlier[model_run, :], 12)
    ohca_meas = np.concatenate((ohca_ey, ohca_ly)).astype(np.float64)
    ohca_meas -= ohca_meas[0]

    # Net RT/TOA (annual)
    rt_later = Dataset(f"{CLIM_PATH}/ESM1-2-LR/rt/global_rt_{ssp}.nc", "r")["rt"][:].__array__()
    rt_earlier = Dataset(f"{CLIM_PATH}/ESM1-2-LR/rt/global_rt_historical.nc", "r")["rt"][:].__array__()
    rt_ly = average_every_n(rt_later[model_run, :], 12)
    rt_ey = average_every_n(rt_earlier[model_run, :], 12)
    new_rtTOA = np.concatenate((rt_ey, rt_ly))

    # CO2 
    new_Co2_df = pd.read_csv(f"{CLIM_PATH}/SSP_inputdata/ERF_ESM1-2-LR_{ssp}.csv")
    eco2 = new_Co2_df["CO2"].values

    return unf_new_opt_depth, ohca_meas, eco2, new_rtTOA

def main(ekf):
    # Define x-axes (edit to match your actual time bases)
    # optical depth is annual-like length new_iter (years)
    # OHCA/RT are annual after average_every_n
    # CO2 length depends on your CSV
    years_annual = np.arange(1850, 2101)  # <-- replace with your actual years array

    fig, axes = plt.subplots(5, 1, figsize=(7, 10), sharex=True)
    axT, axOH, axCO2, axAOD, axRF = axes

    # ----- OBSERVED series (black): wire these in -----
    # obs_temp, obs_ohca, obs_co2, obs_aod, obs_rf should be aligned to same x-axis you plot
    # axT.plot(years_annual, obs_temp, color="black")
    # etc.

    for ssp in ["ssp126", "ssp245", "ssp370"]:
        unf_aod, ohca, eco2, rf = load_mpi_esm1_2_lr_inputs(ssp, new_iter=len(years_annual), ekf=ekf)

        years, tas = load_mpi_tas_annual(CLIM_PATH+"/", ssp, 0)

        
        axT.plot(years_annual, tas, color=colors[ssp], label=ssp.upper())

        axOH.plot(years_annual[:len(ohca)], ohca, color=colors[ssp])
        axCO2.plot(years_annual[:len(eco2)], eco2, color=colors[ssp])
        axAOD.plot(years_annual[:len(unf_aod)], unf_aod, color=colors[ssp])
        axRF.plot(years_annual[:len(rf)], rf, color=colors[ssp])

    axT.plot(years_annual[:174], ekf.temps - np.mean(ekf.temps[0:50]), color='k',label="observed")
    axOH.plot(years_annual[:174], ekf.ocean_heat_measured, color='k')
    meas_CO2 = pd.read_csv("../../../../Common_Data/co2_etheridge_keeling.csv")
    axCO2.plot(years_annual[:174], meas_CO2['mean'][:174], color='k')
    axAOD.plot(years_annual[:174], ekf.data[:, 3] * 0.001, color='k')
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
    fig.suptitle("MPI-ESM1-2-LR Inputs")
    plt.tight_layout()
    #plt.show()
    fig.savefig('MPIESM_inputs.png', dpi=600, bbox_inches='tight')

# You’ll call main(ekf) from wherever ekf is available in your project.

if __name__ == "__main__":
    # Example:
    # plot_all_inputs("mpi_inputs_panels.png")
    # ---- initialize EKF ----

    # ---- plot MPI inputs ----
    main(ekf)
    
