import bayesblend as bb
import numpy as np

# Step 1: Simulate draws for two Gaussian predictive models at 1 time point
n_draws = 1000
mu1, sigma1 = 0.5, 0.8
mu2, sigma2 = 1.15,0.7
y_true = 0.75
sigma_true = 0.2

# Generate draws
post_pred1 = np.random.normal(mu1, sigma1, size=(1,n_draws))
post_pred2 = np.random.normal(mu2, sigma2, size=(1,n_draws))
y_samp = np.random.normal(y_true, sigma_true, size=(1,n_draws))
# Compute log-likelihoods of the true target under each draw
log_lik1 = -0.5 * np.log(2 * np.pi * sigma1**2) \
           - 0.5 * ((y_samp - mu1) ** 2) / sigma1**2
log_lik2 = -0.5 * np.log(2 * np.pi * sigma2**2) \
           - 0.5 * ((y_samp - mu2) ** 2) / sigma2**2

# Step 2: Build Draws objects
draws1 = bb.Draws(log_lik=log_lik1, post_pred=post_pred1)
draws2 = bb.Draws(log_lik=log_lik2, post_pred=post_pred2)

# Step 3: Fit stacking model
stack_model = bb.MleStacking(model_draws={'m1': draws1, 'm2': draws2})
stack_model.fit()

print("Stacking weights:", stack_model.weights)

# Step 4: Blend the predictives
blended = stack_model.predict()
#print("Blended LPD:", blended.lpd)
print("Blended post_pred mean, std:", 
      np.mean(blended.post_pred_2d), 
      np.std(blended.post_pred_2d))
