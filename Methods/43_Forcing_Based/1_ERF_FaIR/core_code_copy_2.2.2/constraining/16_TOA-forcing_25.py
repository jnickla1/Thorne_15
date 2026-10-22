#!/usr/bin/env python
# coding: utf-8
"""
16_TOA-forcing.py
For each of the 100 ensemble members (one per observational GMST realisation),
pick one constrained parameter set at random (reproducibly), run FaIR 2.2.2 with the
same setup as 04_dump-calibration.py, then drive the EBM with total ERF only and save
a 2 x 176 array for timebounds 1850-2025:
  row 0: TOA energy imbalance, W m^-2
  row 1: surface temperature anomaly, K
"""

import os

import numpy as np
import pandas as pd
import xarray as xr
from dotenv import load_dotenv
from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties

load_dotenv()

cal_v = os.getenv("CALIBRATION_VERSION")
fair_v = os.getenv("FAIR_VERSION")
orig_const = os.getenv("CONSTRAINT_SET")
if not all([cal_v, fair_v, orig_const]):
    raise RuntimeError("Missing env vars: CALIBRATION_VERSION, FAIR_VERSION, CONSTRAINT_SET")

cutoff_year = 2025
n_members = 100
t0, t1 = 1750, 2027                  # same time axis as 04_dump-calibration.py
export_start, export_end = 1850, 2026 #fine to have another year
n_export = export_end - export_start + 1   # 176
GLOBAL_SEED = 726371

base = "../../../../../"
outdir = f"{base}output/{orig_const}/toa_imbalance"
os.makedirs(outdir, exist_ok=True)

# -----------------------
# Inputs shared by all members (identical to 04)
# -----------------------
calib_dir = f"{base}output/fair-{fair_v}/v{cal_v}/{orig_const}/calibrations/"
df_methane = pd.read_csv(calib_dir + "CH4_lifetime.csv", index_col=0)
df_landuse = pd.read_csv(calib_dir + "landuse_scale_factor.csv", index_col=0)
df_lapsi = pd.read_csv(calib_dir + "lapsi_scale_factor.csv", index_col=0)

da_emissions = xr.load_dataarray(
    f"{base}output/fair-{fair_v}/v{cal_v}/{orig_const}/emissions/ssps_harmonized_1750-2499.nc"
)
da = da_emissions.loc[dict(config="unspecified", scenario="ssp245")].isel(
    timepoints=slice(0, t1 - t0)
)
emis_1750 = da.isel(timepoints=0)
emitted_species = [
    s for s in da.specie.values if not np.all(np.isnan(da.sel(specie=s).values))
]

trend_shape = np.ones(t1 - t0 + 1)
trend_shape[:271] = np.linspace(0, 1, 271)

species, properties = read_properties()
for s in ["Halon-1202", "NOx aviation", "Contrails"]:
    species.remove(s)

BASELINE_SPECIES = [
    "CH4", "N2O", "NOx", "Sulfur", "CO", "VOC", "BC", "OC", "NH3",
    "CCl4", "CHCl3", "CH2Cl2", "CH3Cl", "CH3Br", "Halon-1211", "SO2F2", "CF4",
]


def apply_fixed_calibrations(fair_obj):
    """Non-sampled calibrations, identical to parallel.py and 04_dump-calibration.py."""
    sc = fair_obj.species_configs
    fill(sc["unperturbed_lifetime"], df_methane.loc["historical_best", "base"], specie="CH4")
    for col, sp in [
        ("CH4", "CH4"), ("N2O", "N2O"), ("VOC", "VOC"), ("NOx", "NOx"),
        ("HC", "Equivalent effective stratospheric chlorine"),
    ]:
        fill(sc["ch4_lifetime_chemical_sensitivity"], df_methane.loc["historical_best", col], specie=sp)
    fill(sc["lifetime_temperature_sensitivity"], df_methane.loc["historical_best", "temp"])
    fill(
        sc["land_use_cumulative_emissions_to_forcing"],
        df_landuse.loc["historical_best", "CO2_AFOLU"],
        specie="CO2 AFOLU",
    )
    fill(sc["lapsi_radiative_efficiency"], df_lapsi.loc["historical_best", "BC"], specie="BC")
    fill(sc["forcing_efficacy"], 0.6, specie="Volcanic")
    for sp in BASELINE_SPECIES:
        fill(sc["baseline_emissions"], float(emis_1750.sel(specie=sp)), specie=sp)


bulk_properties = {
    "bulk": {
        "type": "unspecified",
        "input_mode": "forcing",
        "greenhouse_gas": False,
        "aerosol_chemistry_from_emissions": False,
        "aerosol_chemistry_from_concentration": False,
    }
}

# -----------------------
# Loop over ensemble members
# -----------------------
for member_idx in range(1, n_members + 1):
    constraint_set = f"{orig_const}r{member_idx}cut{cutoff_year}"
    post_dir = f"{base}output/fair-{fair_v}/v{cal_v}/{constraint_set}/posteriors/"

    df_params = pd.read_csv(post_dir + "calibrated_constrained_parameters222.csv", index_col=0)
    # solar trend is not carried into the 222 file (see 04), so read it from the raw dump
    df_raw = pd.read_csv(post_dir + "calibrated_constrained_parameters.csv", index_col=0)

    # reproducible random pick of one constrained parameter set
    rng = np.random.default_rng(GLOBAL_SEED + member_idx)
    pick = int(rng.integers(0, len(df_params)))
    run_id = df_params.index[pick]
    row = df_params.iloc[[pick]].copy()
    solar_trend = float(df_raw.loc[run_id, "fscale_solar_trend"])
    volc_scale = float(row["forcing_scale[Volcanic]"].iloc[0])
    solar_scale = float(row["forcing_scale[Solar]"].iloc[0])

    tmp_csv = f"/tmp/fair_selected_params_{member_idx:03d}.csv"
    row.index = [0]
    row.to_csv(tmp_csv)

    # ---- full-species run to get total ERF ----
    f = FAIR(ch4_method="Thornhill2021")
    f.define_time(t0, t1, 1)
    f.define_scenarios(["all"])
    f.define_species(species, properties)
    f.define_configs([0])
    f.allocate()

    for sp in f.species:
        if sp in emitted_species:
            fill(f.emissions, da.sel(specie=sp).values[:, None], scenario="all", specie=sp)

    f.fill_from_csv(forcing_file=f"{base}data2/forcing/volcanic_solar_2025.csv")
    fill(f.forcing, f.forcing.sel(specie="Volcanic") * volc_scale, specie="Volcanic")
    fill(
        f.forcing,
        f.forcing.sel(specie="Solar") * solar_scale + trend_shape[:, None, None] * solar_trend,
        specie="Solar",
    )

    f.fill_species_configs()
    f.override_defaults(tmp_csv)
    apply_fixed_calibrations(f)

    initialise(f.concentration, f.species_configs["baseline_concentration"])
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)
    f.run(progress=False)

    # ---- bulk run driven by total ERF only ----
    ff = FAIR()
    ff.define_time(t0, t1, 1)
    ff.define_scenarios(["all"])
    ff.define_configs([0])
    ff.define_species(["bulk"], bulk_properties)
    ff.allocate()

    fill(ff.forcing, f.forcing_sum.sel(scenario="all"), specie="bulk", scenario="all")
    ff.override_defaults(tmp_csv)

    initialise(ff.forcing, 0)
    initialise(ff.temperature, 0)
    initialise(ff.ocean_heat_content_change, 0)
    ff.run(progress=False)

    # ---- export timebounds 1850-2025 ----
    toa = ff.toa_imbalance.sel(scenario="all", config=0)
    years = toa["timebounds"].values
    mask = (years >= export_start) & (years <= export_end)
    toa_series = toa.values[mask]
    surf_series = ff.temperature.sel(scenario="all", config=0, layer=0).values[mask]

    if toa_series.shape[0] != n_export or surf_series.shape[0] != n_export:
        raise RuntimeError(
            f"Expected {n_export} years ({export_start}-{export_end}). "
            f"Got N={toa_series.shape[0]}, T={surf_series.shape[0]}"
        )

    outfile = os.path.join(outdir, f"{constraint_set}_toa_imbalance_member{member_idx:03d}.npy")
    np.save(outfile, np.vstack([toa_series, surf_series]))
    os.remove(tmp_csv)
    print(f"Saved {outfile}  (run id {run_id})")
