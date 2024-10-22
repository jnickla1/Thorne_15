import numpy as np


def lowess(x, y, xwidth=10, degree=1, kernel='tricube', retP =False):
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
        dist = np.abs((x[order][i] - x[order])) / xwidth
        w = weight_func(dist)
        
        # Form the system matrix based on the polynomial degree
        if degree == 1:
            A = np.stack([w, x[order] * w]).T
        elif degree == 2:
            A = np.stack([w, x[order] * w, (x[order]**2) * w]).T
        else:
            raise ValueError(f"Unsupported degree {degree}. Use 1 or 2.")
        
        b = w * y[order]
        
        # Solve the linear system
        ATA = A.T.dot(A)
        ATb = A.T.dot(b)
        sol = np.linalg.solve(ATA, ATb)
        
        # Predict for the observation i
        yest = A[i].dot(sol)  # Prediction at x_i
        place = order[i]
        y_sm[place] = yest
        
        # Estimate variance for standard error calculation
        sigma2 = np.sum((A.dot(sol) - y[order])**2) / (N - degree - 1)
        
        # Calculate standard error
        y_stderr[place] = np.sqrt(sigma2 * A[i].dot(np.linalg.inv(ATA)).dot(A[i]))
        # Contribution to the projection matrix diagonal (hat matrix H)
        H_diag[place] = A[i].dot(np.linalg.inv(ATA)).dot(A[i])
        
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

