import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import minimize



# Estimate log-likelihoods for each model's predictions of y_true


# Define our true target distribution (to be approximated)
true_mean = 0.75
true_std = 0.2


# Define two predictive Gaussian models
model_means = [0.5, 1.15]
model_stds = [0.8, 0.7]

# Draw posterior predictive samples from each model
n_samples = 500

for p in range(30):

    y_true = np.random.normal(true_mean, true_std, size=n_samples)

    
    def log_likelihood(weight):
        weight = np.clip(weight, 0, 1)
        mix_pdf = weight * norm.pdf(y_true, model_means[0], model_stds[0]) + \
                  (1 - weight) * norm.pdf(y_true, model_means[1], model_stds[1])
        log_lik = np.sum(np.log(np.clip(mix_pdf, 1e-12, None)))
        return -log_lik  # negative for minimization

    # Find optimal weight to maximize the likelihood
    res = minimize(log_likelihood, x0=0.5, bounds=[(0.0, 1.0)])
    w_opt = res.x[0]


    print(w_opt)

# Plot the results
   # post_preds = [
   #     np.random.normal(loc=model_means[0], scale=model_stds[0], size=n_samples),
   #     np.random.normal(loc=model_means[1], scale=model_stds[1], size=n_samples)
   # ]
    # Draw new stacked predictive samples
    #stacked_samples = w_opt * post_preds[0] + (1 - w_opt) * post_preds[1] #this is the incorrect way to do it
       # y_true2 = np.random.normal(true_mean, true_std, size=n_samples)
if(False):
    plt.hist(y_true, bins=30, density=True, alpha=0.4, label='True')
    plt.hist(post_preds[0], bins=30, density=True, alpha=0.4, label='Model 1')
    plt.hist(post_preds[1], bins=30, density=True, alpha=0.4, label='Model 2')
    #plt.hist(stacked_samples, bins=30, density=True, alpha=0.6, label='Stacked')
    plt.legend()
    plt.title(f"Bayesian stacking, optimal weight = {w_opt:.2f}")

    plt.figure()
    xsp = np.linspace(0,1)
    plt.plot(xsp,[log_likelihood(x) for x in xsp])
    plt.show()
