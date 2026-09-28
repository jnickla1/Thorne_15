import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# ---------- load ----------
best = pd.read_csv("ERF_best_1750-2024.csv", index_col=0)
best.index = best.index.astype(int)

df = pd.read_csv("volcanic_solar.csv", encoding="utf-8-sig")
year_cols = [c for c in df.columns if c not in ["Scenario", "Variable", "Region", "Unit"]]
volc_old = pd.Series({int(y): float(df[df["Variable"] == "Volcanic"].iloc[0][y]) for y in year_cols if y != "2501"})
solar_old = pd.read_csv("old_forcing/solar_erf_timebounds.csv", index_col=0)["erf"]
tsi_obs = pd.read_csv("tsi_cdr_yearly.csv", index_col=0)["tsi"]  # your corrected TSI pull

glossac = pd.read_csv("glossac_aod550.csv", index_col=0).squeeze() / 1000.0  # save your pasted data here
print("GloSSAC 2023-2025 (sanity check, not used for values):")
print(glossac.loc[2023:2025])

# ---------- 1. verify the offset, print evidence ----------
print("=== offset diagnostic (unshifted vs shifted-by-1 mean abs diff vs ERF_best) ===")
for series_name, series in [("SOLAR", solar_old), ("VOLCANIC", volc_old)]:
    for start, end in [(1750,1899),(1900,1999),(2000,2022)]:
        yrs = [y for y in range(start,end+1) if y in best.index and y in series.index]
        unshift = np.abs(series.loc[yrs] - best.loc[yrs, series_name.lower()]).mean()
        syrs = [y for y in yrs if y-1 in best.index]
        shift = np.abs(series.loc[syrs] - best.loc[[y-1 for y in syrs], series_name.lower()].values).mean()
        print(f"{series_name} {start}-{end}: unshifted={unshift:.4f}  shifted={shift:.4f}")

# ---------- 2. replace SOLAR: full 1750-2024 (confirmed clean throughout) ----------
for y in range(1750, 2025):
    df.loc[df["Variable"] == "Solar", str(y)] = best.loc[y, "solar"]

# ---------- 3. replace VOLCANIC: confirmed-clean era only ----------
VOLC_CLEAN_START = 2000  # pre-2000 is a different reconstruction vintage, not a label bug --
                          # decide separately, alongside the full recalibration
for y in range(VOLC_CLEAN_START, 2025):
    df.loc[df["Variable"] == "Volcanic", str(y)] = best.loc[y, "volcanic"]

# ---------- 4. solar 2025: TSI regression, verified ----------
overlap = range(2003, 2015)
tsi_overlap = tsi_obs.loc[overlap].values
erf_overlap = best.loc[overlap, "solar"].values
k_solar, b_solar = np.polyfit(tsi_overlap, erf_overlap, 1)
r2_solar = 1 - np.sum((erf_overlap - (k_solar*tsi_overlap + b_solar))**2) / np.sum((erf_overlap - erf_overlap.mean())**2)
print(f"\nSOLAR: erf = {k_solar:.5f}*tsi + {b_solar:.5f}   R2={r2_solar:.4f}")
solar_2025 = k_solar * tsi_obs.loc[2025] + b_solar
print("2025 ->", solar_2025, " (2024:", best.loc[2024,"solar"], ")")
df.loc[df["Variable"] == "Solar", "2025"] = solar_2025

# ---------- 5. volcanic 2025: GloSSAC power law, verified, weak fit ----------
C = best.loc[2000:2002, "volcanic"].mean()
fit_years = glossac.loc[1979:2024][glossac.loc[1979:2024] < 0.02].index
aod_fit, erf_fit = glossac.loc[fit_years].values, best.loc[fit_years, "volcanic"].values

def model(aod, a, b): return C - a * np.power(aod, b)
popt, _ = curve_fit(model, aod_fit, erf_fit, p0=[5,1], bounds=([0,0],[1000,5]), maxfev=50000)
r2_volc = 1 - np.sum((erf_fit - model(aod_fit,*popt))**2) / np.sum((erf_fit-erf_fit.mean())**2)
print(f"\nVOLCANIC: C={C:.4f} a={popt[0]:.2f} b={popt[1]:.3f}   R2={r2_volc:.4f}  (n={len(fit_years)})")

print("Residuals, recent years (actual - predicted) -- check for the post-Hunga-Tonga pattern:")
for y in [2020,2021,2022,2023,2024]:
    print(" ", y, round(best.loc[y,"volcanic"] - model(glossac[y],*popt), 4))

volcanic_2025 = model(glossac[2025], *popt)
#print("2025 ->", volcanic_2025, " (2024:", best.loc[2024,"volcanic"], ")")
#df.loc[df["Variable"] == "Volcanic", "2025"] = volcanic_2025
print("2025 just has persistenance")
df.loc[df["Variable"] == "Volcanic", "2025"] = df.loc[df["Variable"] == "Volcanic", "2024"]

last_real_year, last_real_val = 2025, float(df.loc[df["Variable"] == "Volcanic", "2024"].iloc[0])   # persistence
taper_years = 9                                             # same span the original file used (2025->2033)
target_year = last_real_year + taper_years                  # 2034 -- reaches 0 here
rate = last_real_val / taper_years

for i, y in enumerate(range(last_real_year + 1, target_year + 1)):
    df.loc[df["Variable"] == "Volcanic", str(y)] = max(last_real_val - rate * (i + 1), 0)

df.to_csv("volcanic_solar_2025.csv", index=False)
print("\nWrote volcanic_solar_2025.csv")

# ---------- 6. long-format files read by sampling/09 (1750-2101, 352 timebounds) ----------
ws = df.set_index("Variable")

solar_long = ws.loc["Solar", [str(y) for y in range(1750, 2102)]].astype(float)
solar_long.index = solar_long.index.astype(int)
solar_long.to_frame("erf").rename_axis("year").to_csv("solar_erf_timebounds.csv")

existing_volc = pd.read_csv("old_forcing/volcanic_ERF_1750-2101_timebounds.csv", index_col=0)["erf"]
new_volc = ws.loc["Volcanic", [str(y) for y in range(VOLC_CLEAN_START, 2102)]].astype(float)
existing_volc.loc[VOLC_CLEAN_START:2101] = new_volc.values
existing_volc.rename_axis("timebounds").to_frame("erf").to_csv("volcanic_ERF_1750-2101_timebounds.csv")

# check exactly what 09 does
s = pd.read_csv("solar_erf_timebounds.csv", index_col="year")["erf"].loc[1750:2101].values
v = pd.read_csv("volcanic_ERF_1750-2101_timebounds.csv", index_col="timebounds")["erf"].loc[1750:2101].values
print("solar length:", len(s), " volcanic length:", len(v), " (both must be 352)")
