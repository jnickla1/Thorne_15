import numpy as np
import warnings

def _normalize_weights(weights: np.ndarray):
    """Normalize weights due to rounding error, witch strict value check"""
    #if not np.isclose(np.sum(weights), 1, atol=1e-7):
    #    raise ValueError(f"Weights do not sum to 1: {weights}.")
    return np.array([max(0, w) for w in weights / sum(weights)])




def predict_non_nan(weights, pred_array, seed=None):
    """
    Blend draws from multiple models using stacking weights while avoiding NaNs.
    Args:
        weights: Dictionary of {model_name: weight_array}.
        seed: Random seed.
    Returns:
        A blended Draws object.
    """

    #model_draws=stack_model.model_draws
    #weights=stack_model.weights
    from collections import defaultdict

    rng = np.random.default_rng(seed)

    #model_keys = list(model_draws.keys())
    M = len(weights)
    S = pred_array.shape[2] #nsamps
    N = pred_array.shape[1] #n_datapoints or years
    #SHAPE = next(iter(model_draws.values())).shape
    
    #weight_array = np.concatenate(list(weights.values())).T[0]  # shape ( M)
    weight_array = weights
    # Create index list: for each sample, which model it should draw from
    blend = np.full((N,S),np.nan)
    
    for i in range(N):
        #check which models have non-nan predictions
        modelsworking = (~np.isnan(pred_array[:,i,:])).all(axis=1)
        neweights = _normalize_weights(weight_array * modelsworking)

        ## generate indices
        inds  = rng.choice(list(range(M)), S, p=neweights)
        ## pick from pred_array forming a weighted mixture.
        blend[i,:] = [pred_array[inds[s],i,s] for s in range(S)]
        
    return blend


def infill_log_likelihoods(model_logliks,penalty=0.2,samps=False):
    """
    Fill in NaNs in a (method, year,samps) matrix of log-likelihoods using scaled sampling.

    model_logliks: np.ndarray of shape (n_models, n_years,samps)

    Returns:
        filled_logliks: np.ndarray with same shape, NaNs replaced
    """
    nsamps = model_logliks.shape[2]
    filled = model_logliks.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if samps:
            mean_loglik_by_modelandyear = np.nanmean(model_logliks, axis=2)
        else:
            mean_loglik_by_modelandyear = model_logliks[:,:,int(nsamps/2)]
    
    for m in range(filled.shape[0]): #for each method
        for y in range(filled.shape[1]):
            if np.isnan(filled[m, y]).all(): #assume all nans or no nans
                # pick a random valid year from same method
                valid_years = np.where(~np.isnan(mean_loglik_by_modelandyear[m]))[0] #pick a valid year
                if valid_years.size == 0:
                    raise ValueError(f"Cannot infill: Method {m} has no valid years.")
                chosen_y = np.random.choice(valid_years)
                #chosen_samp = np.random.randint(0,nsamps,size=nsamps)
                sample_val = filled[m, chosen_y,:]
                #breakpoint()
                # rescale by average across models
                vmodels_infillyear = np.where(~np.isnan(mean_loglik_by_modelandyear[:,y]))[0]
                scale = np.mean(mean_loglik_by_modelandyear[vmodels_infillyear,y]) - \
                        np.mean(mean_loglik_by_modelandyear[vmodels_infillyear,chosen_y])
                filled[m, y] = sample_val + scale - penalty  # optional penalty of 0.2
    return filled

# --- Example usage ---
# In your own script
# from custom_stack_utils import predict_non_nan  # if saved separately


from scipy.optimize import minimize
from scipy.special import logsumexp

def fit_stacking_weights_logp(logp_at_nodes, alpha, w0=None, eps=1e-12):
    """
    Maximize sum_{year,node} alpha[node] * log( sum_k w_k * exp(logp[k,year,node]) )
    Args:
      logp_at_nodes: array (K, T, M) of log densities at quadrature nodes
      alpha: array (M,) quadrature weights (e.g., Gauss–Hermite, normalized)
      w0: optional init (K,), defaults to uniform
    Returns:
      weights (K,) on the simplex
    """
    K, T, M = logp_at_nodes.shape
    A = np.tile(alpha, T)                         # (T*M,)
    X = logp_at_nodes.reshape(K, -1)              # (K, T*M)

    def obj(w):
        w = np.clip(w, eps, 1.0); lw = np.log(w)
        lse = logsumexp(X + lw[:, None], axis=0)  # (T*M,)
        return -np.sum(A * lse)

    def grad(w):
        w = np.clip(w, eps, 1.0); lw = np.log(w)
        logits = X + lw[:, None]
        lse = logsumexp(logits, axis=0)
        probs = np.exp(logits - lse)              # (K, T*M)
        g = -(probs * A).sum(axis=1) / w          # (K,)
        return g

    cons = dict(type='eq', fun=lambda w: np.sum(w) - 1.0)
    bnds = [(0.0, 1.0)] * K
    if w0 is None: w0 = np.full(K, 1.0 / K)
    res = minimize(obj, w0, jac=grad, method="SLSQP", bounds=bnds, constraints=cons)
    return res.x


