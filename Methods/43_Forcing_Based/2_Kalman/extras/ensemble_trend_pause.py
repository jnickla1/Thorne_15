#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
from scipy.stats import t as student_t

experiment_type = "fut_ESM1-2-LR_SSP370_constVolc"


cdataprefix = os.path.expanduser('~/')+"climate_data" + "/"
def average_every_n(arr, n):
    return np.array([np.mean(arr[i:i+n]) for i in range(0, len(arr), n)])


def load_mpi_tas_annual( ssp: str,nmem: int):
    """
    Returns:
      years (int array)
      tas_annual_anom (float array): concatenated historical+future, baseline 1850-1899 = 0
    """
    fut_path  = f"{cdataprefix}ESM1-2-LR/combined/{ssp.lower()}_aave_tas.nc"
    hist_path = f"{cdataprefix}ESM1-2-LR/combined/historical_aave_tas.nc"

    # --- future ---
    with Dataset(fut_path, "r") as ds:
        tas_f_m = np.asarray(ds.variables["tas"][:, :], dtype=np.float64)
        t_f = np.asarray(ds.variables["time"][:], dtype=np.float64) / 365.0 + 1850.0

    # --- historical ---
    with Dataset(hist_path, "r") as ds:
        tas_h_m = np.asarray(ds.variables["tas"][:, :], dtype=np.float64)
        t_h = np.asarray(ds.variables["time"][:], dtype=np.float64) / 365.0 + 1850.0

    years_f = np.floor(average_every_n(t_f, 12)).astype(int)
    years_h = np.floor(average_every_n(t_h, 12)).astype(int)
    years_all = np.concatenate([years_h, years_f])
    tas_ens = np.zeros((nmem,len(years_all)))

    # monthly -> annual
    for i in range(nmem):
        tas_f_y = average_every_n(tas_f_m[i,:], 12)
        tas_h_y = average_every_n(tas_h_m[i,:], 12)
        # baseline 1850–1899 using *historical* first 50 annual values
        baseline = tas_h_y[:50].mean()
        tas_ens[i,:] = np.concatenate([tas_h_y, tas_f_y]) - baseline
        

    return years_all, tas_ens

for this_ssp in [245,370]:
    times, ens370 = load_mpi_tas_annual( f"ssp{this_ssp}" ,50)
    ens370avg = np.mean(ens370,axis=0)

    experiment_type = f"fut_ESM1-2-LR_SSP{this_ssp}_constVolc"
    runmeans = pd.read_csv("../../../../Results/ensemble_mean_"+experiment_type+".csv",header=None)
    runmeans_ens = np.mean(runmeans,axis=0)

    fig, axis = plt.subplots(figsize=(6,6))
    axis.plot(times,ens370avg,label="Ensemble Average")
    axis.plot(times,runmeans_ens,label="20-yr Centered Mean, Ens. Avg", zorder=3,lw=2)
    axis.set_xlim([1970,2050])
    axis.set_ylim([0.25,2])
    axis.set_xlabel("Year")
    axis.set_ylabel("Temperature Δ, (°C)\nrelative to 1850-1900")
    axis.set_title(f"MPI-ESM1-2-LR (SSP{this_ssp}) Pause in Temperature Rise")

    dec_steps = runmeans_ens[(1960-1850):(2080-1850):10]
    print(np.diff(dec_steps))

    plt.plot([2015,2015],[0.75,1.6], color='k', linestyle=':')
    plt.text(2015, 0.6, 'Switch from\nHistorical to\nSSP simulation', rotation=0, 
             verticalalignment='center', horizontalalignment='center', fontsize=12, color='k')

    # 1) Fit window
    searly = 1970
    wind=35
    fit_mask = (times >= searly ) & (times < searly +wind)
    t_fit = np.asarray(times[fit_mask], float)
    y_fit = np.asarray(runmeans_ens[fit_mask], float)

    # 2) Center time for numerical stability
    t0 = searly +wind/2  # midpoint of 1970–2005
    x_fit = t_fit - t0

    # 3) Design matrix for quadratic: y = b0 + b1*x + b2*x^2
    X = np.column_stack([np.ones_like(x_fit), x_fit, x_fit**2])

    # 4) OLS fit
    beta, *_ = np.linalg.lstsq(X, y_fit, rcond=None)

    # Residuals + variance
    y_fit_hat = X @ beta
    resid = y_fit - y_fit_hat
    n = len(y_fit)
    p = X.shape[1]          # 3 params
    dof = n - p
    s2 = (resid @ resid) / dof
    s = np.sqrt(s2)

    # Covariance of beta
    XtX_inv = np.linalg.inv(X.T @ X)
    cov_beta = s2 * XtX_inv

    # 5) Projection grid (1970–2050)
    t_proj = np.arange(1970, 2051, 1)
    x_proj = t_proj - t0
    Xp = np.column_stack([np.ones_like(x_proj), x_proj, x_proj**2])

    y_proj = Xp @ beta

    # Standard error of the MEAN prediction at each x
    # se_mean^2 = diag(Xp * cov_beta * Xp^T)
    se_mean = np.sqrt(np.einsum("ij,jk,ik->i", Xp, cov_beta, Xp))

    # (Optional) Standard error for PREDICTION interval (includes residual variance)
    se_pred = np.sqrt(se_mean**2 + s2)

    # t critical value
    alpha = 0.05
    tcrit = student_t.ppf(1 - alpha/2, dof)

    # 95% CI for mean curve
    ci_lo = y_proj - tcrit * se_mean
    ci_hi = y_proj + tcrit * se_mean

    # 95% prediction interval (optional; wider)
    pi_lo = y_proj - tcrit * se_pred
    pi_hi = y_proj + tcrit * se_pred

    # fitted curve + CI shading
    plt.plot(t_proj, y_proj, color="k", lw=2, label="Quad. fit: 20-yr cent Ens Avg (1970-2004)")
    plt.fill_between(t_proj, ci_lo, ci_hi, alpha=0.25,color='k')

    #DO IT AGAIN AT LATER TIME

    # 1) Fit window
    slate = 2025
    fit_mask = (times >= slate) & (times < slate+wind)
    t_fit = np.asarray(times[fit_mask], float)
    y_fit = np.asarray(runmeans_ens[fit_mask], float)

    # 2) Center time for numerical stability
    t0 = slate +wind/2 # midpoint of 1970–2005
    x_fit = t_fit - t0

    # 3) Design matrix for quadratic: y = b0 + b1*x + b2*x^2
    X = np.column_stack([np.ones_like(x_fit), x_fit, x_fit**2])

    # 4) OLS fit
    beta, *_ = np.linalg.lstsq(X, y_fit, rcond=None)

    # Residuals + variance
    y_fit_hat = X @ beta
    resid = y_fit - y_fit_hat
    n = len(y_fit)
    p = X.shape[1]          # 3 params
    dof = n - p
    s2 = (resid @ resid) / dof
    s = np.sqrt(s2)

    # Covariance of beta
    XtX_inv = np.linalg.inv(X.T @ X)
    cov_beta = s2 * XtX_inv

    # 5) Projection grid (1970–2050)
    t_proj = np.arange(1970, 2051, 1)
    x_proj = t_proj - t0
    Xp = np.column_stack([np.ones_like(x_proj), x_proj, x_proj**2])

    y_proj = Xp @ beta

    # Standard error of the MEAN prediction at each x
    # se_mean^2 = diag(Xp * cov_beta * Xp^T)
    se_mean = np.sqrt(np.einsum("ij,jk,ik->i", Xp, cov_beta, Xp))

    # (Optional) Standard error for PREDICTION interval (includes residual variance)
    se_pred = np.sqrt(se_mean**2 + s2)

    # t critical value
    alpha = 0.05
    tcrit = student_t.ppf(1 - alpha/2, dof)

    # 95% CI for mean curve
    ci_lo = y_proj - tcrit * se_mean
    ci_hi = y_proj + tcrit * se_mean

    # 95% prediction interval (optional; wider)
    pi_lo = y_proj - tcrit * se_pred
    pi_hi = y_proj + tcrit * se_pred

    # fitted curve + CI shading
    plt.plot(t_proj, y_proj, color="r", lw=2, label="Quad. fit: 20-yr cent Ens Avg (2025-2059)")
    plt.fill_between(t_proj, ci_lo, ci_hi, alpha=0.25,color='r')
    plt.legend()
    fig.savefig(f"MPIESM_pause_{this_ssp}.png", dpi=600, bbox_inches='tight')
plt.show()
