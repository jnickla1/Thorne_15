import numpy as np
import pandas as pd

# --------------------------------------------------------
# Load CSV results (evaluation script)
# --------------------------------------------------------

# Make sure the method filter is correct:
method = "EBMKF_ta4"


# --------------------------------------------------------
# Load NPZ results (distributed iterative script)
# --------------------------------------------------------
npz_vals = []
csv_vals = []
for ens in range(60):
    fname = f"../Results3/nmes2run_exp3_r{ens}.npz"
    try:
        data = np.load(fname)
        # rmse75 is typically stored as array with EMBKF_ta4 as first entry
        rmse75 = float(data["rmse75"][0])  
        npz_vals.append(rmse75)
    except FileNotFoundError:
        npz_vals.append(np.nan)
        print(f"Missing: {fname}")

    csv = pd.read_csv("../Results/current_fut_statistics_fut_NorESM_RCP45_Volc"+str(ens)+".csv")
    csv_vals.append(csv[csv["method_name"] == method]["75RMS"].reset_index(drop=True)[0])
npz_vals = np.array(npz_vals)
csv_vals = np.array(csv_vals)
# --------------------------------------------------------
# Build comparison dataframe
# --------------------------------------------------------
df_compare = pd.DataFrame({
    "ens": range(60),
    "npz_rmse75": npz_vals,
    "csv_rmse75": csv_vals,     # in case CSV has more rows
})

df_compare["abs_diff"] = np.abs(df_compare["npz_rmse75"] - df_compare["csv_rmse75"])

print(df_compare)

# --------------------------------------------------------
# Optional: print only mismatches
# --------------------------------------------------------
print("\nEntries with differences > 1e-8:")
print(df_compare[df_compare["abs_diff"] > 1e-8])
ensmean=np.sqrt(np.mean(npz_vals**2))
print(ensmean)
