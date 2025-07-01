import numpy as np
import matplotlib.pyplot as plt

# Generate example data
np.random.seed(42)
x = np.linspace(1, 20, 20)
y = 2 * x + np.random.normal(scale=4, size=len(x))  # Linear trend with noise

# Create figure and subplots
fig, axes = plt.subplots(4, 1, figsize=(4, 10), sharex=True, sharey=True)
#fig.suptitle("Different Statistical Methods for Fitting a Line")

methods = ["Kalman Filter: Evolve","Monte Carlo: Samples","Linear Algebra: Least Squares Fit","Bayesian: Refine Prior"]

# Plot scatterplot in each panel
for i, ax in enumerate(axes):
    ax.scatter(x, y, color='darkgrey', label="Data")
    ax.set_title(methods[i])
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
# Label axes
#axes[0].set_ylabel("Y values")
#for ax in axes:
#    ax.set_xlabel("X values")


import sys
sys.path.append(os.path.expanduser("~/")+"data/jnickla1/Thorne_15_codefigurestats/Methods/42_Temp_Alone/5_Kalman/")
import Kal_flexLin_method as KF

z0=[0,KF.q00,KF.R000]
kf_output,_,_,_,_ = KF.KF_compute(x,y,y,z0)
axes[0].plot(x,kf_output)



# Monte Carlo Simulation: Generate 50 random lines
num_samples = 20
slopes = np.random.normal(2, 1, num_samples)
intercepts = np.random.normal(0, 5, num_samples)
errors = []

for m, b in zip(slopes, intercepts):
    y_pred = m * x + b
    error = np.sum((y - y_pred) ** 2)
    errors.append(error)

# Normalize errors for shading (lower error = darker line)
norm_errors = (errors - np.min(errors)) / (np.max(errors) - np.min(errors))
colors = [(1, e, e) for e in norm_errors]  # White to red gradient
for j in range(num_samples):
    axes[1].plot(x, slopes[j] * x + intercepts[j], color=colors[j], alpha=0.5,zorder=0)


from scipy import stats
#Perform Ordinary Least Squares (OLS) Regression
def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor

n=len(x)
regres = stats.linregress(x, y) #can also give stdrr of slope, intercept
slope = regres.slope
intercept =regres.intercept
sstdrr = regres.stderr
y_pred = slope * x + intercept
residuals = y - y_pred
# Calculate standard error of the estimate (s_e)
s_e = np.sqrt(np.sum(residuals**2) / (n - 2))

# Critical t-value for 1se confidence level
t_crit = stats.t.ppf(.5 + 0.3413, df=n-2) #one standard error, two sided

cis = confidence_interval(x, y_pred, np.mean(x), s_e, t_crit, n )

axes[2].plot(x, y_pred, color='green', label="OLS Fit")
axes[2].plot(x, y_pred+2*cis, color='green',lw=0.5, label="95% Confidence Interval")
axes[2].plot(x, y_pred-2*cis, color='green',lw=0.5, label="95% Confidence Interval")

axes[0].set_ylim(0,40)
plt.tight_layout()

##plt.figure()
##import pymc3 as pm
##with pm.Model() as model:
##    # Priors for the parameters
##    slope = pm.Normal('slope', mu=1, sd=10)
##    intercept = pm.Normal('intercept', mu=-5, sd=10)
##    sigma = pm.HalfNormal('sigma', sd=1)
## 
##    # Expected value of the outcome
##    mu = intercept + slope * x
## 
##    # Likelihood (sampling distribution) of the observations
##    Y_obs = pm.Normal('Y_obs', mu=mu, sd=sigma, observed=y)
## 
##    # Run the MCMC sampling
##    trace = pm.sample(2000, tune=1000)
## 
### Plot the posterior distributions
##pm.plot_posterior(trace, var_names=['slope', 'intercept', 'sigma'])

#prior
y_prior = .5 * x + 15
densf=2
axes[3].fill_between(x,y_prior-5*cis, y_prior+5*cis, color='orange',alpha=0.9/densf, zorder=0,linewidth=0.0)
axes[3].fill_between(x,y_prior-10*cis, y_prior+10*cis, color='orange',alpha=0.6/densf, zorder=0,linewidth=0.0)
axes[3].fill_between(x,y_prior-15*cis, y_prior+15*cis, color='orange',alpha=0.3/densf, zorder=0,linewidth=0.0)

y_pred2 = y_pred *4/5 + y_prior/5
axes[3].fill_between(x,y_pred2-2*cis, y_pred2+2*cis, color='sienna',alpha=0.9, zorder=0,linewidth=0.0)
axes[3].fill_between(x,y_pred2-4*cis, y_pred2+4*cis, color='sienna',alpha=0.6, zorder=0,linewidth=0.0)
axes[3].fill_between(x,y_pred2-6*cis, y_pred2+6*cis, color='sienna',alpha=0.3, zorder=0,linewidth=0.0)

plt.show()

plt.show()
