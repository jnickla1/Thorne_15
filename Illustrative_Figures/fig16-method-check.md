# Fig. 16 method check (Sect. 4.4), Oct 2026

Checked each method's schematic against its papers and its operational code (the version feeding Fig. 15 / IGCC 2024). The current figure build is v8.

## Gear colour rule (v8)
- A gear's own colour shows where the forcing, or the response to it, comes from: black for observed or assessed series, purple for CMIP6 model output.
- A red ring means its amplitude is fitted to observations.
- An all-red gear (EBM-KF "?") is fitted noise with no data source.
- How the rule applies to each method:
  - Linear: black gear with red ring (observed CO2, fitted slope).
  - EBM-KF: black gears with no ring (literature parameters, not fitted).
  - FaIR and GWI: black gears with red rings (scaling factors or regression coefficients).
  - ROF: purple gears with red rings (β scaling of ESM fingerprints).
  - KCC: purple gears with red rings for the ESM components. Black gears with red rings for the natural forcings, which come from observed forcing through an EBM and are then updated by observations.
- v8 fixes: the ROF sun gear changed from black to purple, and the KCC volcanic gear changed from purple to black.
- Open question for John: the dashed red ring on the EBM-KF cloud feedback is not a gear, so the rule does not cover it.

## Linear (Jarvis & Forster 2024)
- GMST is regressed on CO2 concentration above pre-1700 (linear in ppm, not log) using weighted least squares. Weights come from HadCRUT5 uncertainty.
- Residuals are modelled as AR(1) (Cochrane–Orcutt). Uncertainty is propagated with a 10^4-sample Monte Carlo.
- There is no explicit feedback or flux dynamics: non-CO2 forcing and natural variability are absorbed into the slope or the residual.
- Figure change: removed the Planck feedback and the upward longwave arrow.
- Text suggestion for 4.4.1: replace "simple linear regression" with "weighted least squares with AR(1) errors".

## EBM-KF (Nicklas et al. 2025)
- The schematic matches the paper:
  - aerosol optical depth transmittance
  - atmospheric/cloud albedo, including anthropogenic cloud forcing
  - surface-albedo feedback
  - T^4 longwave
  - water-vapour feedback
  - CMIP5 exchange coefficient
- Internal variability enters as R_const, the covariance of annual observations minus the 30-yr mean.
- Attribution in this paper assimilates single-forcing ESM runs (see the Fig. 15 caption). That is now noted in the figure.
- Resolved: the red "?" gears are Kalman process noise, which adds or removes heat from the climate state.
- Resolved: the water-vapour and albedo feedbacks come from IPCC AR6 Ch. 7 values, which are largely CMIP6-based, so purple stays.

## FaIR (fair-calibrate v1.4.1; attribution with FaIR 2.2.2)
- Contrails and aviation NOx are removed from the species list in v1.4.x. The contrail gear was moved to the land surface to stand for land-use albedo and light-absorbing particles on snow, both of which are scaled.
- `stochastic_run=True` in the prior ensemble: red-noise forcing (sigma_eta, gamma) plus temperature noise (sigma_xi), calibrated to CMIP6 abrupt-4xCO2 runs.
- Resolved: `04_dump-calibration.py` (FaIR 2.2.2) reads `calibrated_constrained_parameters222.csv`, which has `stochastic_run=True`, per the FaIR update chat of 12 Oct 2026.
- Constraint chain:
  1. ~1.5M prior members
  2. GMST RMSE < 0.17 K against GMST 1850–2022
  3. Reweighting to 8 assessed metrics: ECS, TCR, OHC, 2003–2022 warming, ERFari, ERFaci, ERFaer, CO2 concentration
  4. 841 posterior members

## GWI (Walsh code, Forster et al. 2025)
- `model_choice='FaIR_V2'`: a three-box impulse response with 28 CMIP6 emulations. Median timescales are about 1, 8 and 320 yr.
- The previous figure's "4.1 / 209 yr" are the Haustein et al. (2017) two-box values, used only by the old AR5_IR option.
- 3-way OLS on GHG, other human forcing and natural forcing, with a constant term.
- Internal variability: CMIP6 piControl segments, drift-filtered, regressed in parallel. Their coefficients are added to the observational fit.
- Monte Carlo ensemble of order 10^7.

## ROF (ESMValTool forster23 branch; IGCC 2024 supplement)
- Regressors: historical+SSP2-4.5 (ALL), hist-nat and hist-GHG, each with its ssp245 extension. OTH is computed as ALL − NAT − GHG.
- hist-aer is processed but not regressed. Gillett et al. (2021) used hist-aer in their primary analysis.
- Resolved: the paper's ROF data come from ClimateIndicator/anthropogenic-warming-assessment (IGCC), so hist-aer is dropped from the figure.
- Three models: CanESM5, IPSL-CM6A-LR and MIROC6.
- Data: 5-yr blended-masked GMST, compared with non-infilled HadCRUT5.
- Fit: total least squares with a Ledoit–Wolf regularised covariance.
- Internal variability comes from intra-ensemble anomalies, not piControl.
- Text suggestion for 4.4.4.3: add GHG to the list of responses used.

## KCC (Ribes et al. 2021; IGCC 2024 supplement)
- 22 CMIP6 models with a GAM-smoothed prior on ALL, NAT and GHG:
  - ALL from historical + SSP runs
  - NAT from a rescaled EBM response
  - GHG from hist-GHG, or rebuilt from 1%CO2 runs where hist-GHG is missing
- Internal variability: a mixture of two AR(1) processes, fitted to observations minus the CMIP6 multi-model mean.
- One-step Gaussian (kriging) update: observations constrain ALL directly, and NAT and GHG shift through the prior covariance.

## Draft caption (for v8)
Figure 16. Schematic comparison of the attributable-warming methods in Sect. 4.4. Each row is one method. The columns show, from left to right, how the forced response is modelled, how internal variability is represented, and how observations constrain the result.

Left column (dynamics): Gears are forcing time series. A gear's colour shows where the forcing, or the modelled response to it, comes from (black, observed or assessed; purple, CMIP6 models), and a red ring marks a forcing whose amplitude is fitted to observations. Circular arrows are feedbacks, with their sign. Straight arrows are energy fluxes (blue, shortwave; orange, longwave), wavy symbols are heat exchange between layers, and I-bars show layer heat capacities. Braces group forcings into simulation experiments (ROF, KCC) or regression components (GWI). Thin connecting lines mark assumed relationships between mechanisms, such as the chemistry links in FaIR. The red "?" gears in EBM-KF are Kalman-filter process noise, which adds or removes heat from the modelled climate state. Gridded blocks are CMIP6 coupled Earth system models, with stacked outlines marking multi-model ensembles. GWI uses FaIR's three-box thermal response with parameters emulating 28 CMIP6 models (median timescales shown). In ROF, the response to other human forcings is taken as ALL − NAT − GHG. In KCC, the natural response comes from an energy balance model, and the GHG response from hist-GHG runs or, where these are unavailable, 1 % CO2 runs.

Middle column (internal variability): dots are residuals of observations, dashed lines are synthetic noise from a statistical model, and solid lines are coupled-model runs (thin, individual members; thick, ensemble mean).

Right column (observational constraint): illustrative panels of each estimation step, with observations as grey dots: least-squares fitting, Kalman filtering (inset: forecast, observation and update), Monte Carlo sampling, and Bayesian updating of a prior (inset: prior and posterior distributions).

Colours throughout show the type of evidence: observed or measured (black), fitted to observations (red), derived from CMIP6 models (purple), or an assumed statistical form (grey).
