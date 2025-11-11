import numpy as np


def lowess(x, y, xwidth=10, degree=1, kernel='tricube', retP =False, ridge=0):
    """
    Enhanced LOWESS smoother with uncertainty and support for degree-1/2 polynomials and different kernels.
    James Brennan, https://james-brennan.github.io/posts/lowess_conf/
    and chatGPT to extend to degree 2
    Parameters:
        x : array-like
            Independent variable data.
        y : array-like
            Dependent variable data.
        xwidth : int, optional
            Nuber of datapoints to use in each local regression (bandwidth).
        degree : int, optional
            Degree of the polynomial fit (1 for linear, 2 for quadratic).
        kernel : str, optional
            Type of kernel weighting ('tricube' or 'gaussian').
        retP :
            return return the equivalent number of parameters p

        ridge: add for numerical stability, say 1e-12
    
    Returns:
        y_sm : array
            Smoothed values.
        y_stderr : array
            Standard error of the smoothed values.
    """
    # Get some parameters
    N = len(x)  # Number of observations
    # Don't assume the data is sorted
    order = np.argsort(x)
    # Storage for results
    y_sm = np.zeros_like(y)
    y_stderr = np.zeros_like(y)
    H_diag = np.zeros(N)  # Diagonal of the projection (hat) matrix
    
    # Define weighting functions
    if kernel == 'tricube':
        weight_func = lambda d: np.clip((1 - np.abs(d)**3)**3, 0, 1)
    elif kernel == 'gaussian':
        weight_func = lambda d: np.exp(-0.5 * (d)**2)
    else:
        raise ValueError(f"Unsupported kernel '{kernel}'. Use 'tricube' or 'gaussian'.")
    
    # Run the regression for each observation i
    for i in range(N):
        # weights for the local window around x_i
        dist = np.abs((x[order][i] - x[order])) / xwidth
        w = weight_func(dist)
        sqrtw = np.sqrt(w)

        # Build *unweighted* design X, then weight via sqrt(w)
        if degree == 1:
            X = np.column_stack([np.ones(N), x[order]])
        elif degree == 2:
            X = np.column_stack([np.ones(N), x[order], x[order]**2])
        else:
            raise ValueError("degree must be 1 or 2")

        A = sqrtw[:, None] * X          # = W^{1/2} X
        b = sqrtw * y[order]            # = W^{1/2} y

        # Solve (Xᵀ W X) β = Xᵀ W y  via normal equations on A
        XtWX = A.T @ A + ridge * np.eye(A.shape[1]) #ridge for stability
        XtWy = A.T @ b
        beta = np.linalg.solve(XtWX, XtWy)

        # Prediction at x_i (NO weights in xrow)
        if degree == 1:
            xrow = np.array([1.0, x[order][i]])
        else:
            xi = x[order][i]
            xrow = np.array([1.0, xi, xi*xi])

        yhat_i = xrow @ beta
        place = order[i]
        y_sm[place] = yhat_i

        # Residuals in data space, weighted by w
        yhat_all = (X @ beta)
        res = y[order] - yhat_all
        p = X.shape[1]
        sigma2 = (w * res**2).sum() / (N - p)

        # Standard error of the smoothed estimate at x_i
        XtWX_inv = np.linalg.inv(XtWX)
        y_stderr[place] = np.sqrt(sigma2 * (xrow @ XtWX_inv @ xrow))

        # Hat-matrix diagonal (optional)
        H_diag[place] = w[i] * (xrow @ XtWX_inv @ xrow)

        
    if not(retP):
        return y_sm, y_stderr
    else:
        return y_sm, y_stderr, H_diag


def run_method(years, temperature, uncert, model_run, experiment_type):
    avg_len_l=3

    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    for i in range(avg_len_l, len(years)):
        chunktmps=temperature[0:i+1]
        chunkyrs=years[0:i+1]
        lfit,lsde = lowess(chunkyrs, chunktmps, xwidth=10, degree=1, kernel='tricube')
        means[i] = lfit[-1]
        ses[i] = lsde[-1]
    return means, ses, lfit, lsde

