#!/usr/bin/env python
# coding: utf-8

"""Global annual biomass burning emissions from GFED5.1 (1997-2022)."""

# GFED5.1 monthly files already contain per-species emissions in g per grid cell per
# month, so we only sum over months and grid cells. This replaces the GFED4.1s
# reconstruction (DM x partitioning x emission factors).
# NOx is reported "as NO" (same as GFED4), so the NO -> NO2 conversion in
# 07_make-slcf-emissions.py remains correct.
# NMVOC is taken as GFED5's total gaseous NMOC (NMOC_g).

import os

import pandas as pd
import xarray as xr
from dotenv import load_dotenv
from tqdm.auto import tqdm

load_dotenv()

print("Getting GFED5.1 data...")

cal_v = os.getenv("CALIBRATION_VERSION")
fair_v = os.getenv("FAIR_VERSION")
constraint_set = os.getenv("CONSTRAINT_SET")
progress = os.getenv("PROGRESS", "False").lower() in ("true", "1", "t")
gfed5_dir = os.getenv("GFED5_DIR")  # folder containing GFED5.1_monthly_YYYY.nc

start_year = 1997
end_year = 2022

# GFED5 variable name -> column name expected by 07 and 10
species_map = {
    "BC": "BC",
    "OC": "OC",
    "SO2": "SO2",
    "NOx": "NOx",
    "CO": "CO",
    "NMOC_g": "NMVOC",
    "NH3": "NH3",
    "CH4": "CH4",
    "N2O": "N2O",
}

rows = {}
for year in tqdm(range(start_year, end_year + 1), disable=1 - progress):
    fname = os.path.join(gfed5_dir, f"GFED5.1_monthly_{year}.nc")
    with xr.open_dataset(fname) as ds:
        n_months = ds.sizes["time"]
        if n_months != 12:
            raise ValueError(f"{fname} has {n_months} months, expected 12")
        # g per grid cell per month -> Tg per year
        rows[year] = {
            out_name: float(ds[gfed_name].sum(skipna=True)) / 1e12
            for gfed_name, out_name in species_map.items()
        }

gfed_df = pd.DataFrame.from_dict(rows, orient="index")
gfed_df.index.name = "year"
print(gfed_df)

os.makedirs(
    f"../../../../../output/fair-{fair_v}/v{cal_v}/{constraint_set}/emissions",
    exist_ok=True,
)
# filename kept as before so 07 and 10 need no changes
gfed_df.to_csv(
    f"../../../../../output/fair-{fair_v}/v{cal_v}/{constraint_set}/emissions/"
    "gfed4.1s_1997-2022.csv"
)
