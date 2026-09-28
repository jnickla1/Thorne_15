import pandas as pd

cmip7 = pd.read_csv("historical_emissions_1750-2025.csv")
META_COLS = ["model", "scenario", "region", "variable", "unit"]
year_cols = [c for c in cmip7.columns if c not in META_COLS]

cmip7 = cmip7[cmip7["variable"] != "CO2"].copy()  # drop the redundant combined-CO2 row

iamc_name = {v: f"Emissions|{v}" for v in cmip7["variable"].unique()}
iamc_name["CO2 FFI"] = "Emissions|CO2|Energy and Industrial Processes"
iamc_name["CO2 AFOLU"] = "Emissions|CO2|AFOLU"

cmip7["Variable"] = cmip7["variable"].map(iamc_name)
cmip7["Unit"] = cmip7["unit"].replace("kt HFC4310/yr", "kt HFC4310mee/yr")
cmip7["Model"] = "Historical"
cmip7["Scenario"] = "GCP+CEDS+PRIMAP+GFED"  # same label 10_combine-emissions.py already uses
cmip7["Region"] = "World"

year_rename = {c: str(int(float(c))) for c in year_cols}
cmip7 = cmip7.rename(columns=year_rename)
plain_year_cols = sorted(year_rename.values(), key=int)

out = cmip7[["Model", "Scenario", "Region", "Variable", "Unit"] + plain_year_cols]
out.to_csv("all_1750-2025.csv", index=False)
