"""Read-only check: effect of the 1850-1899 -> 1850-1900 (inclusive) preindustrial baseline on Figure 24.

Changing the baseline subtracts a constant c from a member's whole temperature series, so the 1.5 °C crossing
in the new frame is the (1.5 + c) crossing in the stored (old-frame) series. For a method whose output shifts by
exactly c when its input does (running means, LOWESS, GAM, Kalman, ...), comparing method and benchmark at
threshold 1.5 + c is exact. FaIR/GWI (fixed bias constants) and EBMKF are not exactly shift-equivariant, so for them
this is a first-order approximation.

Shifts:
  MPI-ESM1-2-LR member m: c_m = mean(hist[0:51]) - mean(hist[0:50]) of that member's raw annual series.
  NorESM (aligned to HadCRUT5 obs over 1980-2000): c_obs = same quantity for HadCRUT5.csv, identical for every member.
  Ensemble-mean benchmark (used by eFdyrs15): shifts by the mean of the member shifts.

Fdyrs15 / eFdyrs15 are recomputed exactly as evaluate_future_scenarios.py does (p >= 0.5 first crossing at 1/12-yr
resolution, 25-yr window), from the stored current central estimates and SEs in Results/resultsfut_*.pickle.
The reported change is F(threshold + c) - F(threshold) with the same code, added to the stored CSV value.
Sampled methods (FaIR, GWI, CGWL) are stored as mean/SE, so their p-curve is Gaussian-approximated in both terms.
Writes only into this folder.
"""
import sys, os, pickle
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from netCDF4 import Dataset

ROOT = Path('/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats'); OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT)); import config
SEL = ['CGWL10y_sfUKCP', 'EBMKF_ta4', 'FaIR_comb_unB', 'GAM_AR1', 'GWI_tot_CGWL', 'Kal_flexLin', 'lowess1dt36wnc']
SELA = ['FaIR_nonat_unB', 'GWI_anthro_CGWL']
EXPS = [('fut_ESM1-2-LR_SSP126_constVolc', 50), ('fut_ESM1-2-LR_SSP245_constVolc', 50), ('fut_ESM1-2-LR_SSP370_constVolc', 50),
        ('fut_NorESM_RCP45_Volc', 60), ('fut_NorESM_RCP45_VolcConst', 60)]
INUM, WIN = 12, 25
THR = np.arange(1.1, 3, .1)[4]   # the exact float the pipeline uses for "1.5"

def annual(x): return np.array([np.mean(x[i:i + 12]) for i in range(0, len(x), 12)])

obs = pd.read_csv(ROOT/'Common_Data/HadCRUT5.csv')['Anomaly'].to_numpy()
c_obs = np.mean(obs[0:51]) - np.mean(obs[0:50])

def target_cross(years, std, thr):
    fine = np.arange(years[0], years[-1] + 1/INUM, 1/INUM); s0 = np.interp(fine, years, std)
    s = s0[~np.isnan(s0)]; fc = fine[~np.isnan(s0)][1:]
    hit = fc[np.logical_and(s[:-1] < thr, s[1:] >= thr)]
    return -1/INUM/2 + hit[0] if len(hit) else np.nan

def method_first(central, se, thr, closest, lastyr):
    evalyrs = np.arange(int(closest) - WIN, min(int(closest) + WIN, lastyr))
    p = stats.norm.cdf((central[evalyrs - 1850] - thr) / se[evalyrs - 1850])
    fnn = np.argmax(~np.isnan(p)); p[:fnn] = 0
    lnn = len(p) - np.argmax(~np.isnan(p)[::-1]) - 1; p[lnn + 1:] = 1
    fine = np.arange(evalyrs[0], evalyrs[-1] + 1/INUM, 1/INUM); pi = np.interp(fine, evalyrs, p)
    cr = np.logical_and(pi[:-1] < 0.5, pi[1:] >= 0.5); mid = fine[:-1]/2 + fine[1:]/2
    d = mid[cr] - closest
    return d[0] if len(d) else 15.0

rows = []
for exp, n in EXPS:
    mpi = 'ESM1-2-LR' in exp
    if mpi:
        h = np.asarray(Dataset(f"{config.CLIMATE_DATA_PATH}/ESM1-2-LR/combined/historical_aave_tas.nc")["tas"][:].__array__(), dtype=np.float64)
        cs = np.array([np.mean(annual(h[m])[0:51]) - np.mean(annual(h[m])[0:50]) for m in range(n)])
    else:
        cs = np.full(n, c_obs)
    c_ens = cs.mean()
    ens = np.nanmean(pd.read_csv(ROOT/f'Results/ensemble_mean_{exp}.csv', header=None).to_numpy(), axis=0)
    lastyr = 2101 if mpi else 2100
    for m in range(n):
        r = pickle.load(open(ROOT/f'Results/results{exp}{m}.pickle', 'rb'))
        std = r['cent20y']['LT_trend'][2]; years = 1850 + np.arange(len(std))
        csv = pd.read_csv(ROOT/f'Results/current_fut_statistics_{exp}{m}.csv').set_index('method_name')
        for meth, col, base in [(x, 'Fdyrs15', std) for x in SEL] + [(x, 'eFdyrs15', ens) for x in SELA]:
            cen, se = np.asarray(r[meth]['LT_trend'][0], float), np.maximum(np.asarray(r[meth]['LT_trend'][1], float), 1e-12)
            shift = cs[m] if col == 'Fdyrs15' else c_ens
            t0, t1 = target_cross(years, base, THR), target_cross(years, base, THR + shift)
            f0 = method_first(cen, se, THR + 0, t0, lastyr)          # old baseline (reproduction)
            f1 = method_first(cen, se, THR + cs[m], t1, lastyr)       # new baseline: method sees its own member shift
            stored = float(csv.loc[meth, col])
            rows.append(dict(experiment=exp, member=m, method=meth, column=col, c_member=cs[m], c_target=shift,
                             stored=stored, repro_old=f0, recomputed_new=f1, delta=f1 - f0, new_value=stored + (f1 - f0)))
    print(exp, f'c mean {cs.mean():+.5f} min {cs.min():+.5f} max {cs.max():+.5f}', flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT/'fig24_baseline1900_per_member.csv', index=False)
print('c_obs (NorESM) =', round(c_obs, 6))
