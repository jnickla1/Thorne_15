import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import theilslopes

IN_PATH = Path.cwd() / "historical_emissions_1750-2024.csv"
OUT_PATH = Path.cwd() / "historical_emissions_1750-2025.csv"
META_COLS = ["model", "scenario", "region", "variable", "unit"]

df = pd.read_csv(IN_PATH)

# 1. year columns -> FaIR's half-year timepoint convention (1750 -> 1750.5, etc.)
year_cols = [c for c in df.columns if c not in META_COLS]
df = df.rename(columns={c: float(c) + 0.5 for c in year_cols})
year_cols = sorted(float(c) + 0.5 for c in year_cols)  # wait: already renamed, so:
year_cols = [c for c in df.columns if c not in META_COLS]

# 2. index by species so extension + override land on the right rows
df = df.set_index("variable")

def extend_one_year(df, year_cols, target_year=2025, window=5):
    years = np.array(sorted(year_cols))
    recent_years = years[-window:]
    out = {}
    for var, row in df.iterrows():
        y = row[recent_years].values.astype(float)
        slope, intercept, _, _ = theilslopes(y, recent_years)
        out[var] = intercept + slope * (target_year + 0.5)
    return pd.Series(out, name=target_year + 0.5)

new_col = extend_one_year(df, year_cols)
new_col.loc["CO2 FFI"] = 38.1     # GCB 2025 preliminary β€” swap in the exact ESSD figure
new_col.loc["CO2 AFOLU"] = 4.1    # GCB 2025 preliminary β€” swap in the exact ESSD figure
df[2025.5] = new_col

df.loc["HFC-4310mee", "unit"] = "kt HFC4310mee/yr"  # cosmetic; see note below

# 3. restore the FaIR-ready shape (scenario, region, variable, unit, years...) β€” no `model` column,
# matching historical_1750-2023.csv's layout. Double-check this against your actual
# data2/emissions/v1.4.1/historical_1750-2024.csv, which I haven't seen directly.
out_df = (df.reset_index()
            .drop(columns=["model"])
            [["scenario", "region", "variable", "unit"] + sorted(year_cols) + [2025.5]])

# scenario is currently "historical" in the raw file β€” check what your ready v1.4.1 file
# actually uses (likely "all", per the older historical_1750-2023.csv convention) and rename
# if needed: out_df["scenario"] = "all"

out_df.to_csv(OUT_PATH, index=False)
