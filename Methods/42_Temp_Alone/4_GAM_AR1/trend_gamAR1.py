import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf
from statsmodels.tsa.filters.filtertools import convolution_filter
from pygam import LinearGAM, s

# Read in the data
dat = pd.read_csv("HadCRUT5.csv")

# Plot settings for two subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

# Fit GAM assuming residuals are uncorrelated
times = dat['Time'].values[:, np.newaxis]  # Reshape to 2D
G = LinearGAM(s(0)).fit(times, dat['Anomaly'])
pt_fit = G.predict(times)
pt_se = G.confidence_intervals(times, width=0.95)

# Plot the results for the uncorrelated GAM
axs[0].scatter(dat['Time'], dat['Anomaly'], color='gray', label='Temperature Anomaly', s=20)
axs[0].plot(dat['Time'], pt_fit, lw=2, label='GAM Fit')
axs[0].fill_between(dat['Time'], pt_se[:, 0], pt_se[:, 1], color='red', alpha=0.3, label='95% CI')
axs[0].set_xlabel("Year")
axs[0].set_ylabel("Temp Anomaly")
axs[0].set_title("GAM with AR1=0")
axs[0].legend()

# Add 20-year running mean
running_mean = convolution_filter(dat['Anomaly'].values, np.ones(20)/20, nsides=2)
axs[0].plot(dat['Time'], running_mean, color='blue', label='20-year Running Mean')

# Estimate AR(1) coefficient from residuals
residuals_G = dat['Anomaly'] - pt_fit
r1 = acf(residuals_G, nlags=1, fft=False)[1]

# Create a GAM object to be used in the gamAR1 implementation
def gamAR1(G,y0, rho, method="ML"):
    """ Function to handle AR1 autocorrelation correction for GAM. """

    n = len(y0)
    X0=G.feature_matrix(y0)

    # Inverse correlation matrix coefficients
    a = 1 / (1 - rho ** 2)
    b = 1 + 2 * rho ** 2 * a
    c = -rho * a

    # Coeffs for cholesky factor R
    f = -np.sqrt(a - 1)
    e = c / f
    f = f / e
    e = e / e

    # Apply transformations to the data
    y = np.concatenate([y0[:n - 1] * e + y0[1:] * f, [y0[n - 1]]])
    X = np.vstack([X0[:n - 1] * e + X0[1:] * f, X0[n - 1]]) 

    # Fit the model
    #G._model_data.y = y
    #G._model_data.X = X
    G2= LinearGAM(s(0)).fit(y, X)

    # Adjust the log-likelihood
    if method in ["ML", "REML"]:
        G2.statistics_['pseudo_r2']['explained_deviance'] -= n * np.log(e)

    return G2

# Apply the AR(1) correction and re-fit the model
G_ar1 = gamAR1(G, times,r1)
pt_fit_ar1 = G_ar1.predict(dat['Time'])
pt_se_ar1 = G_ar1.confidence_intervals(dat['Time'], width=0.95)

# Plot the results for GAM with AR(1)
axs[1].scatter(dat['Time'], dat['Anomaly'], color='gray', label='Temperature Anomaly', s=20)
axs[1].plot(dat['Time'], pt_fit_ar1, lw=2, label='GAM Fit with AR1')
axs[1].fill_between(dat['Time'], pt_se_ar1[:, 0], pt_se_ar1[:, 1], color='red', alpha=0.3, label='95% CI')
axs[1].set_xlabel("Year")
axs[1].set_ylabel("Temp Anomaly")
axs[1].set_title("GAM with estimated AR1")
axs[1].legend()

# Add 20-year running mean to AR(1) plot
axs[1].plot(dat['Time'], running_mean, color='blue', label='20-year Running Mean')

plt.tight_layout()
plt.show()
