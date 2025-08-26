import numpy as np, pandas as pd, glob
from scipy.stats import gaussian_kde
import sys

scenarios = ["fut_ESM1-2-LR_SSP126_constVolc","fut_ESM1-2-LR_SSP245_constVolc","fut_ESM1-2-LR_SSP370_constVolc",
                   "fut_NorESM_RCP45_Volc","fut_NorESM_RCP45_VolcConst"]                     # list of 5 names
approaches = ["inverse_variance","sharpened_blended"]
nruns = 60                           # total runs per scenario
lhund = -100                          # your existing var
nnmethods = int(sys.argv[1]) 
comparison_type = 'ff'

# Collect per-scenario stacks
firstcross15_sum = np.full((5,2), np.nan)
firstcross2_sum = np.full((5,2), np.nan)
firstcross5_sum = np.full((5,2), np.nan)
ncrosses_sum     = np.full((5,2), np.nan)
rmse_sum         = np.full((5,2), np.nan)
kl_array         = np.full((5,2), np.nan)
rmse75_sum         = np.full((5,2), np.nan)
kl75_array         = np.full((5,2), np.nan)

for exp_index, exp in enumerate(scenarios):
    # load all runs for this scenario
    files = sorted(glob.glob(f"Results3/nmes{nnmethods}run_exp{exp_index}_r*.npz"))
    ncrosses_stack = []; rmse_stack = []; fc_diffs = []; kl_stack=[]
    rmse75_stack=[]; kl75_stack=[]
    print(exp+" "+str(len(files)))
    for f in files:
        z = np.load(f)
        ncrosses_stack.append(z["ncrosses"])            # (2,)
        rmse_stack.append(z["rmse"])                    # (2,)
        rmse75_stack.append(z["rmse75"]) 
        fc_diffs.append(z["firstcross15_diff"])         # (years,2)
        #breakpoint()
        kl_stack.append(z["kl"])
        kl75_stack.append(z["kl75"])
    ncrosses_stack = np.stack(ncrosses_stack, axis=0)   # (nruns,2)
    rmse_stack     = np.stack(rmse_stack, axis=0)       # (nruns,2)
    rmse75_stack     = np.stack(rmse75_stack, axis=0)
    fc_diffs       = np.stack(fc_diffs, axis=0)         # (nruns,2)
    kl_stack       = np.stack(kl_stack,axis=0)
    kl75_stack       = np.stack(kl75_stack,axis=0)
    # KDE over runs for first-cross metric
    runs_len = 50 if exp.startswith("fut_ESM1") else 60
    for m in range(2):
        kde = gaussian_kde(fc_diffs[:runs_len, m].ravel())
        firstcross15_sum[exp_index, m] = kde.integrate_box_1d(-1, 1)
        firstcross2_sum[exp_index, m] = kde.integrate_box_1d(-2, 2)
        firstcross5_sum[exp_index, m] = kde.integrate_box_1d(-5, 5)
    # Averages across runs (your original logic)
    ncrosses_sum[exp_index] = np.nanmean(ncrosses_stack, axis=0)
    rmse_sum[exp_index]     = np.sqrt(np.nanmean(rmse_stack**2, axis=0))
    kl_array[exp_index]     = np.nanmean(kl_stack, axis=0) / (-10 - lhund)
    rmse75_sum[exp_index]     = np.sqrt(np.nanmean(rmse75_stack**2, axis=0))
    kl75_array[exp_index]     = np.nanmean(kl75_stack, axis=0) / (-10 +75)
    
# Tidy → wide CSVs
arrays = {"firstcross1": firstcross15_sum,"firstcross2":firstcross2_sum,"firstcross5":firstcross5_sum, "ncrosses": ncrosses_sum,
          "rmse": rmse_sum, "kl": kl_array,"rmse75": rmse75_sum, "kl75": kl75_array}
records = []
for metric, arr in arrays.items():
    for i, exp in enumerate(scenarios):
        for j, approach in enumerate(approaches):
            records.append({"experiment": exp, "approach": approach,
                            "metric": metric, "value": arr[i, j]})
df = pd.DataFrame(records)
wide = df.pivot_table(index=["experiment","approach"], columns="metric",
                      values="value").reset_index()
wide.columns.name = None
wide.to_csv(f"final_combined_{nnmethods}methods_wide_{comparison_type}.csv", index=False)

