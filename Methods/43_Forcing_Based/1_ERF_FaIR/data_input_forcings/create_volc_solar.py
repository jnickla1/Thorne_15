import numpy as np
import pandas as pd
from pathlib import Path

# ---------- 1. VOLCANIC: splice in the real recent-years product ----------
volc_new = pd.read_csv("volcanic_ERF_annual_timebounds_1750-2025.csv", index_col=0)["volcanic_ERF"]

# GloSSAC sanity check only (not used to compute values) -- confirms no missed eruption signal
glossac = pd.read_csv("glossac_aod550.csv", index_col=0).squeeze() / 1000.0  # save your pasted data here
print("GloSSAC 2023-2025 (sanity check, not used for values):")
print(glossac.loc[2023:2025])

# ---------- 2. SOLAR: recover TSI->ERF scale from the overlap, apply to real TSI ----------
tsi_obs = pd.read_csv("tsi_cdr_yearly.csv", index_col=0)["tsi"]  # from the NOAA CDR download

solar_existing = pd.read_csv("solar_erf_timebounds.csv", index_col=0)["erf"]

# Overlap window: years both series cover AND where solar_existing is still real
# (verify this cutoff against your own inspection -- 2014/2015 looked like the likely splice point)
overlap = range(2003, 2015)
tsi_overlap = tsi_obs.loc[overlap].values
erf_overlap = solar_existing.loc[overlap].values

# ERF is a linear function of TSI (ERF = k*(TSI - TSI_ref)), so ordinary least squares recovers k, offset
k, intercept = np.polyfit(tsi_overlap, erf_overlap, 1)
print(f"Recovered TSI->ERF: erf = {k:.5f} * tsi + {intercept:.5f}")

fitted = k * tsi_overlap + intercept
resid = erf_overlap - fitted
print("Overlap fit R2:", 1 - np.sum(resid**2) / np.sum((erf_overlap - erf_overlap.mean())**2))
# expect this near 1.0 -- it's a physical linear relationship, not an empirical power law like volcanic

solar_new = k * tsi_obs.loc[2023:2025] + intercept
print("\nUpdated solar ERF, 2023-2025:")
print(solar_new)
print("vs. old template values:")
print(solar_existing.loc[2023:2025])

# ---------- 3. Write both back into the wide volcanic_solar.csv ----------
df = pd.read_csv("volcanic_solar.csv", encoding="utf-8-sig")
year_cols = [c for c in df.columns if c not in ["Scenario", "Variable", "Region", "Unit"]]

for y in [2023, 2024, 2025]:
    df.loc[df["Variable"] == "Volcanic", str(y)] = volc_new.loc[y]
    df.loc[df["Variable"] == "Solar", str(y)] = solar_new.loc[y]

last_real_year, last_real_val = 2025, volc_new.loc[2025]   # 0.237954, from the found file
taper_years = 9                                             # same span the original file used (2025->2033)
target_year = last_real_year + taper_years                  # 2034 -- reaches 0 here
rate = last_real_val / taper_years

for i, y in enumerate(range(last_real_year + 1, target_year + 1)):
    df.loc[df["Variable"] == "Volcanic", str(y)] = max(last_real_val - rate * (i + 1), 0)
# years after target_year: leave whatever's already there (0), matching the original convention

import numpy as np

k_physical = 0.175  # (1 - albedo)/4, albedo ~0.3; adjust if your source cites a specific efficacy-adjusted value
tsi_ref = np.mean(tsi_overlap - erf_overlap / k_physical)
print(f"erf = {k_physical} * (tsi - {tsi_ref:.4f})")

fitted = k_physical * (tsi_overlap - tsi_ref)
resid = erf_overlap - fitted
print("R2 with fixed slope:", 1 - np.sum(resid**2) / np.sum((erf_overlap - erf_overlap.mean())**2))

solar_new = k_physical * (tsi_obs.loc[2023:2025] - tsi_ref)
print(solar_new)


df.to_csv("volcanic_solar_updated.csv", index=False)
