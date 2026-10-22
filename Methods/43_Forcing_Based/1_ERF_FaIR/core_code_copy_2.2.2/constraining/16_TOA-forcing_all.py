#!/usr/bin/env python++++
# coding: utf-8
"""
16_TOA-forcing.py  (per-member version)
--------------------------------------
Run once per ensemble member (100 total). For the given ENSEMBLE_RUN index,
select one calibrated parameter set at random (reproducibly), drive FaIR 2.2.2
with *total* ERF only, and save a 2x175 .npy array for years 1850–2024:
  row 0: TOA energy imbalance, W m^-2
  row 1: Surface temperature anomaly, K
"""

import os
import sys
import numpy as np
import pandas as pd
from dotenv import load_dotenv

load_dotenv()

# -----------------------
# Env + inputs
# -----------------------
cal_v      = os.getenv("CALIBRATION_VERSION")
fair_v     = os.getenv("FAIR_VERSION")
orig_const = os.getenv("CONSTRAINT_SET")
#ensemble_run = os.getenv("ENSEMBLE_RUN")
cutoff_year  = 2025 #os.getenv("CUTOFF_YEAR")

for ensemble_run in range(1,101):

    if not all([cal_v, fair_v, orig_const, ensemble_run, cutoff_year]):
        raise RuntimeError("Missing required env vars: CALIBRATION_VERSION, FAIR_VERSION, CONSTRAINT_SET, ENSEMBLE_RUN, CUTOFF_YEAR")

    try:
        member_idx = int(ensemble_run)  # 0..99 expected
    except Exception as e:
        raise RuntimeError(f"ENSEMBLE_RUN must be an integer (0..99). Got {ensemble_run!r}") from e

    if not (0 <= member_idx <= 100):
        raise RuntimeError(f"ENSEMBLE_RUN={member_idx} out of range [0,100]")

    constraint_set = f"{orig_const}r{member_idx}cut{cutoff_year}"  # mirrors earlier naming

    # Calibrated parameters CSV created by 04_dump-calibration.py
    calib_csv = f"../../../../../output/fair-{fair_v}/v{cal_v}/{constraint_set}/posteriors/calibrated_constrained_parameters222.csv"
    if not os.path.exists(calib_csv):
        raise FileNotFoundError(f"Calibrated parameters file not found: {calib_csv}")

    # Output directory and file
    outdir = f"../../../../../output/{orig_const}/toa_imbalance"
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{constraint_set}_toa_imbalance_member{member_idx:03d}.npy")

    # -----------------------
    # Load calibrated draws, pick one deterministically
    # -----------------------
    df_params = pd.read_csv(calib_csv, index_col=0)
    n_draws = len(df_params)
    if n_draws < 1:
        raise RuntimeError("No calibrated parameter rows found in CSV.")

    # One temp-history per CSV in this per-member workflow, so all rows belong to this member.
    start, stop = 0, n_draws

    # Reproducible selection: base seed + member index
    GLOBAL_SEED = 726371
    rng = np.random.default_rng(GLOBAL_SEED + member_idx)
    pick = int(rng.integers(low=start, high=stop))
    row = df_params.iloc[[pick]].copy()  # keep as DataFrame with single row
    config_label = row.index[0]  # use the CSV index label as the FaIR config id

    # Write a temp CSV containing ONLY the selected row, reindexed to 0
    tmp_csv = os.path.join('/tmp', f'fair_selected_params_{member_idx:03d}.csv')
    row_zero = row.copy()
    row_zero.index = [0]
    row_zero.to_csv(tmp_csv)

    # -----------------------
    # Build total ERF ("all") and run FaIR
    # -----------------------
    from fair import FAIR
    from fair.interface import fill, initialise
    from fair.io import read_properties

    t0, t1 = 1750, 2026
    export_start = 1850

    # Full-species to construct total forcing
    f = FAIR()
    f.define_time(t0, t1, 1)
    f.define_scenarios(["all"])

    species, properties = read_properties("../../../../../data2/calibration/v1.4.1/species_configs_properties_1.4.1.csv")
    f.define_species(species, properties)

    # Single config (index 0)
    f.define_configs([0])
    f.allocate()

    # Fill historical emissions and volcanic/solar forcing
    f.fill_from_csv(forcing_file="../../../../../data2/forcing/volcanic_solar_2025.csv")
    fill(
        f.forcing,
        f.forcing.sel(specie="Volcanic") * df_configs["forcing_scale[Volcanic]"].values.squeeze(),
        specie="Volcanic",
    )
    fill(
        f.forcing,
        f.forcing.sel(specie="Solar") * df_configs["forcing_scale[Solar]"].values.squeeze()
        + trend_shape[:, None, None] * solar_trend_values,
        specie="Solar",
    )

    # Apply per-config forcing scale overrides if present
    volc_key = "forcing_scale[Volcanic]"
    sol_key  = "forcing_scale[Solar]"
    trend_shape = np.ones((t1 - t0 + 1))
    trend_shape[:271] = np.linspace(0, 1, 271)

    def maybe_series(df, key, default):
        return df[key].values.squeeze() if key in df.columns else default

    solar_trend_col = next((c for c in row.columns if "solar" in c.lower() and "trend" in c.lower()), None)
    solar_trend_val = maybe_series(row, solar_trend_col, 0.0) if solar_trend_col else 0.0

    if volc_key in row.columns:
        fill(f.forcing, f.forcing.sel(specie="Volcanic") * row[volc_key].values.squeeze(), specie="Volcanic")
    if sol_key in row.columns:
        fill(
            f.forcing,
            f.forcing.sel(specie="Solar") * row[sol_key].values.squeeze() + trend_shape[:, None, None] * solar_trend_val,
            specie="Solar",
        )

    # Load per-species config defaults and override with selected calibration row
    f.fill_species_configs("../../../../../data2/calibration/v1.4.1/species_configs_properties_1.4.1.csv")
    f.override_defaults(tmp_csv)  # map by row index from CSV

    # Initial conditions
    initialise(f.concentration, f.species_configs["baseline_concentration"])
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)

    # Run to get per-species and summed forcing
    f.run()

    # Bulk model: drive with total forcing only
    ff = FAIR()
    ff.define_time(t0, t1, 1)
    ff.define_scenarios(["all"])
    ff.define_configs([0])

    bulk_species = ["bulk"]
    bulk_properties = {
        "bulk": {
            "type": "unspecified",
            "input_mode": "forcing",
            "greenhouse_gas": False,
            "aerosol_chemistry_from_emissions": False,
            "aerosol_chemistry_from_concentration": False,
        }
    }
    ff.define_species(bulk_species, bulk_properties)
    ff.allocate()

    # Total forcing from full run
    fill(ff.forcing, f.forcing_sum.sel(scenario="all"), specie="bulk", scenario="all")

    # Use the same selected calibration params for EBM3 climate response
    ff.override_defaults(tmp_csv)

    # Initial conditions
    initialise(ff.forcing, 0)
    initialise(ff.temperature, 0)
    initialise(ff.ocean_heat_content_change, 0)

    # Run to populate ff.toa_imbalance and temperatures
    ff.run()

    # -----------------------
    # Extract TOA imbalance and surface temperature for 1850–2024
    # -----------------------
    N_all = ff.toa_imbalance.sel(scenario="all")  # [timebounds, config]
    years = N_all["timebounds"].values
    mask = (years >= export_start) & (years <= t1)

    # Helper to extract surface temperature (layer 0 or box 0 depending on FaIR build)
    temp = ff.temperature.sel(scenario="all")
    # Handle possible dimension names
    layer_dim = None
    for cand in ["layer", "box"]:
        if cand in temp.dims:
            layer_dim = cand
            break

    if layer_dim is None:
        # assume leading dimension is layer
        surf = temp.isel({temp.dims[0]: 0})
    else:
        surf = temp.sel({layer_dim: 0})

    # Drop config dimension (single) and select years
    toa_series = N_all.sel(config=0).values.squeeze()[mask]
    surf_series = surf.sel(config=0).values.squeeze()[mask]

    if toa_series.shape[0] != 176 or surf_series.shape[0] != 176:
        raise RuntimeError(f"Expected 175 years (1850–2025). Got N={toa_series.shape[0]}, T={surf_series.shape[0]}")

    # Stack into 2x175 array
    out = np.vstack([toa_series, surf_series])

    # Save
    np.save(outfile, out)

    print(f"Saved {outfile}")
