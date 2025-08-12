import numpy as np
from numpy.polynomial.hermite import hermgauss
from custom_stack_utils import predict_non_nan, infill_log_likelihoods, fit_stacking_weights_logp
import pdb;

# Step 1: Simulate draws for two Gaussian predictive models at 1 time point
n_draws = 5000
mu1, sigma1 = 0.5, 0.8
mu2, sigma2 = 1.15,0.7
y_true = 0.75
sigma_true = 0.2

# Generate draws
pred1 = np.random.normal(mu1, sigma1, size=(1,n_draws))
#post_pred2 = np.append(np.random.normal(mu2, sigma2, size=(1,n_draws)),np.nan).reshape(1,-1)
post_pred2 = np.random.normal(mu2, sigma2, size=(2,n_draws))

post_pred1 = np.concatenate([pred1, np.full((1,n_draws), np.nan)], axis=0)


x, w = hermgauss(100)
z = np.sqrt(2.0) * x
alpha = w / np.sqrt(np.pi)  # normalized weights

# Precompute y nodes for the target distribution
y_nodes = (y_true + sigma_true * z).reshape(1, -1)


# Compute log-likelihoods of the true target under each draw
loglik1 = -0.5 * np.log(2 * np.pi * sigma1**2) - 0.5 * ((y_nodes - mu1) ** 2) / sigma1**2
log_lik1 = np.concatenate([loglik1, np.full((1,100), np.nan)], axis=0)

loglik2 = -0.5 * np.log(2 * np.pi * sigma2**2)  - 0.5 * ((y_nodes - mu2) ** 2) / sigma2**2
log_lik2 = np.concatenate([loglik2, loglik2], axis=0)

comb_likes = np.stack((log_lik1, log_lik2),axis=0) #(2, 2, 500)
new_comb_likes = infill_log_likelihoods(comb_likes,penalty=0)


# Step 2: Build Draws objects
#draws1 = bb.Draws(log_lik=new_comb_likes[0].reshape(1, -1), post_pred=post_pred1.reshape(1, -1))
#draws2 = bb.Draws(log_lik=new_comb_likes[1].reshape(1, -1), post_pred=post_pred2.reshape(1, -1))


# Step 3: Fit stacking model
#stack_model = bb.MleStacking(model_draws={'m1': draws1, 'm2': draws2})
#stack_model.fit()
weights = fit_stacking_weights_logp(new_comb_likes,alpha)
print("Stacking weights:", weights)


# Step 4: Blend the predictives
pred_array = np.stack((post_pred1, post_pred2),axis=0) #(2, 2, 500)
blended = predict_non_nan(weights, pred_array)
#print("Blended LPD:", blended.lpd)
print("Blended 0 post_pred mean, std:", 
      np.mean(blended[0,:]), 
      np.std(blended[0,:]))

print("Blended 1 post_pred mean, std:", 
      np.mean(blended[1,:]), 
      np.std(blended[1,:]))
