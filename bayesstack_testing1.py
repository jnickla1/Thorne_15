import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.stats import norm
from scipy.optimize import minimize

# True target (narrow)
true_mean, true_std = 0.75, 0.2

# Two Gaussian models
model_means = np.array([0.5, 1.15])
model_stds  = np.array([0.8, 0.7])

# Gauss–Hermite nodes/weights for Z~N(0,1) expectation
# E[f(Z)] = (1/√π) * Σ w_i f(√2 * x_i)
#n_gh = 40

for n_gh in [30,40,50,75,100]:
    x, w = hermgauss(n_gh)
    z = np.sqrt(2.0) * x
    alpha = w / np.sqrt(np.pi)  # normalized weights

    # Precompute y nodes for the target distribution
    y_nodes = true_mean + true_std * z

    # Precompute model pdfs at nodes (shape: [n_models, n_nodes])
    p_at_nodes = np.vstack([
        norm.pdf(y_nodes, loc=model_means[0], scale=model_stds[0]),
        norm.pdf(y_nodes, loc=model_means[1], scale=model_stds[1]),
    ])

    def obj(w_scalar):
        w_scalar = np.clip(w_scalar, 0.0, 1.0)
        mix_pdf = w_scalar * p_at_nodes[0] + (1.0 - w_scalar) * p_at_nodes[1]
        # stable log with tiny floor
        return -np.sum(alpha * np.log(np.maximum(mix_pdf, 1e-300)))

    res = minimize(obj, x0=[0.5], bounds=[(0.0, 1.0)], method="L-BFGS-B")
    w_opt = float(res.x)
    print("Optimal weight on model 1:", w_opt)
