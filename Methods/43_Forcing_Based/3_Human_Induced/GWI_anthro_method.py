import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
from scipy.integrate import quad


#import pdb;

from scipy.interpolate import CubicSpline, PPoly
from scipy.stats import rv_continuous

def gen_orig_number(new_member_number,sz_ens):
    #fix lexicographic reshuffling
    nums = np.arange(1, sz_ens+1)
    reshuffled = sorted([f"{n}|" for n in nums])
    recovered_order = [int(s.rstrip("|")) for s in reshuffled]
    if new_member_number==-1:
        return recovered_order
    else:
        return recovered_order[new_member_number]

#from time import time
#start = time()

##def invert_function_newton(f, f_prime, y, x0, tol=1e-4, max_iter=4):
##    """Inverts a function f using Newton's method.
##
##    Args:
##        f: The function to invert.
##        y: The value for which we want to find the inverse.
##        x0: Initial guess for the inverse.
##        tol: Tolerance for the error.
##        max_iter: Maximum number of iterations.
##
##    Returns:
##        The approximate inverse value.
##    """
##
##    #def f_prime(x):
##    #    # Approximate the derivative using a small step size
##    #    h = 1e-5
##    #    return (f(x + h) - f(x)) / h
##
##    x = x0
##    for _ in range(max_iter):
##        x_new = x - (f(x) - y) / f_prime(x)
##        if abs(x_new - x) < tol:
##            return x_new
##        x = x_new
##
##    return x_new
##    #raise ValueError("Newton's method did not converge")


def create_extended_spline(x, y):
    # Step 1: Fit the cubic spline to the data
    cs = CubicSpline(x, y, bc_type='natural')
    # Step 2: Extract the coefficients from the cubic spline
    coeffs = cs.c  # Shape is (4, n-1) for cubic spline segments
    # Step 3: Define the left and right extrapolation slopes and intercepts
    left_slope = cs(x[0],nu=1)  # Slope at the leftmost interval
    left_intercept = y[0] - left_slope * (1)  # y = m*x + b form so b = y-mx  
    right_slope = cs(x[-1],nu=1)  # Slope at the rightmost interval
    right_intercept = y[-1] #- right_slope * (x[-1]+1/2)
    # Step 4: Combine the cubic segments with two linear ones for extrapolation   
    left_poly = np.array([[0, 0, left_slope, left_intercept]])
    right_poly = np.array([[0, 0, right_slope, right_intercept]])
    # Step 5: Concatenate the coefficients for the full PPoly object
    extended_coeffs = np.hstack((left_poly.T, coeffs, right_poly.T))
    extended_x = np.concatenate(([x[0] - 1], x, [x[-1] + 1]))  # Extend x with extrapolation regions
    # Step 6: Create the PPoly object
    extended_ppoly = PPoly(extended_coeffs, extended_x)
    return extended_ppoly


percentiles = np.array([5, 17, 50, 83, 95]) / 100  # Convert to 0-1 range
# Define the custom distribution class that fits a cubic spline to get the cdf
class CustomSplineDistribution(rv_continuous):
    def __init__(self, ordvalues, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.median = ordvalues[2]
        self.est_stdev = (ordvalues[3] - ordvalues[1])/2
        self.values = ordvalues  # We need values to define the support for the CDF/PPF
        zgrid= stats.norm.ppf( percentiles)
        self.adjspline = create_extended_spline(zgrid, ((ordvalues-self.median)/self.est_stdev)) #now a natural cubic spline
        self.mean = self._genmean()
        self.std = self._genstd()
        #breakpoint()



    # Quantile function (PPF): Inverse of the CDF
    def ppf(self, q):
        y = (self.adjspline(stats.norm.ppf( q )))*self.est_stdev 
        return (y+ self.median)

    #def ppfprime(self, q):
    #    h = 1e-5
    #    y = (self.adjspline(stats.norm.ppf( q ),nu=1)*(stats.norm.ppf( q +h ) - stats.norm.ppf( q -h ))/2/h)*self.est_stdev 
    #    return y

    # CDF: Inverse of the quantile function (PPF)
    def cdf(self, x):
        # We now directly invert the spline to get the cumulative distribution
        if (np.isnan(x).all()):
            return np.nan
        else:
            xarr = np.atleast_1d(x) #allow this to take in a whole array
            x_inner = np.full(np.shape(xarr),np.nan)
            for i in range(len(xarr)):
                try:
                    x_inner[i] = self.adjspline.solve( (xarr[i]-self.median)/self.est_stdev , extrapolate=True)[0] #should have only one root
                except:
                    breakpoint()
            if np.iterable(x):
                return stats.norm.cdf(x_inner)
            else:
                return stats.norm.cdf(x_inner[0])
       # return np.clip(self.adjspline.solve(ielf.ppf, self.ppfprime,x, q_guess), 0, 1)  # use newton to invert the ppf function

    
##    def pdf0(self, x):
##        if (np.isnan(x)):
##            return np.nan
##        else:
##            try:
##                x_inner = self.adjspline.solve( (x-self.median)/self.est_stdev , extrapolate=True)[0] #should have only one root
##            except:
##                breakpoint()
##            return stats.norm.pdf(x_inner) did not properaly take the derivative of the inner part

       
# PDF: Derivative of the CDF (approximated numerically)
    def pdf(self, x,alp0=1e-5):
        alp= self.est_stdev * alp0
        if (np.isnan(x).all()):
            return np.nan
        else:
            return (self.cdf(x+alp)- self.cdf(x-alp))/2/alp
        #eps = 1e-3  # Small epsilon for numerical derivative

    def pdfn0(self, x):
        p=self.pdf(x)
        if(np.iterable(p)):
            replace_locs = (p< 1e-10)
            p[replace_locs] = stats.norm.pdf(x[replace_locs],loc = self.mean,scale = self.std)
            return p

        else:
            if(p< 1e-10):
                 return stats.norm.pdf(x,loc = self.mean,scale = self.std) + 1e-300
            else:
                return p

    def cdfn0(self, x):
        p=self.cdf(x)
        if(np.iterable(p)):
            replace_locs = ((p< 1e-10) or ((1-p)<1e-10))
            p[replace_locs] = stats.norm.cdf(x[replace_locs],loc = self.mean,scale = self.std)
            return p
        
        else:
            if( (p< 1e-10) or ((1-p)<1e-10)):
                return stats.norm.cdf(x,loc = self.mean,scale = self.std) + 1e-300
            else:
                return p


        

    def _genmean(self):
    
        # Define the integrand for expected value calculation
        def integrand(x):
            return x * self.pdf(x)
        
        # Integrate within a ±3 standard deviation range for efficiency
        lower_bound = self.median - 5 * self.est_stdev
        upper_bound = self.median + 5 * self.est_stdev
        
        # Perform the integration with adaptive quadrature
        result, _ = quad(integrand, lower_bound, upper_bound, limit=50)
        return result
        
    def _genstd(self):
        # Define the integrand for expected value calculation
        def variance_integrand(x):
            return ((x -self.mean) ** 2) * self.pdf(x)
        # Integrate within a ±3 standard deviation range for efficiency
        lower_bound = self.median - 5 * self.est_stdev
        upper_bound = self.median + 5 * self.est_stdev
        
                # Perform the integration to calculate variance
        variance, _ = quad(variance_integrand, lower_bound, upper_bound,limit=50)

        # Standard error is the square root of variance
        se = np.sqrt(variance)
        
        return se
    
cur_path = os.path.dirname(os.path.realpath(__file__))

#NOW IN GWI_anthro_orig
##gwi_levels_retro0 = pd.read_csv(cur_path+"/GWI_data/GWI_full_info.csv", header=[0, 1])
##gwi_levels_retro =gwi_levels_retro0.iloc[1:,] #blank row
##gwi_r =gwi_levels_retro['Ant'].to_numpy()
##gwi_levels_curr0 = pd.read_csv(cur_path+"/GWI_data/GWI_hist_only.csv", header=[0, 1])
##gwi_levels_curr =gwi_levels_curr0.iloc[1:,]
##gwi_c =gwi_levels_curr['Ant'].to_numpy()
##lyearr = np.shape(gwi_r)[0]
##lyearc = np.shape(gwi_c)[0]
##lyear = lyearr + 1850 
ord_ind = [0,1,4,2,3]
min_fact = 4 #the std error of some of these gets way too small in places, so when it goes below this times smaller than the temps_1std we correct it.

def run_method(years, temperature, uncert, model_run, experiment_type):
    
   # print(f"Run_method start {time() - start:.2f}s")

    #(temps_CIl, temps_CIu) = uncert
    #temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    
    #nsamps = len(forec.columns)-2
    #samp_cur = np.full((np.shape(years)[0],nsamps) ,np.nan)
#ONLY DEFINED NOW FOR CURRENT WARMING LEVEL, so k must be 0 or will return just NANs
    syear = 1950
    if experiment_type == 'historical':
        gwi_levels_curr0 = pd.read_csv(cur_path+"/Thorne2025_GWI_Results/ANNUAL_ESM1-2-LR/"+
                "GWI_results_ANNUAL_HISTORICAL-ONLY_SCENARIO--observed-SSP245_ENSEMBLE-MEMBER-"+
                                       "-all_VARIABLES--GHG-Nat-OHF___REGRESSED-YEARS--1850-1950_to_1850-2024.csv", header=[0, 1])
        curbias = 0.0082        
    else:
        #future case, grabbing the ANNUAL resutls
        exp_attr = experiment_type.split("_") #fut_ESM1-2-LR_SSP126_constVolc #

        if (exp_attr[1]=='ESM1-2-LR'):
            gwi_levels_curr0 = pd.read_csv(cur_path+"/Thorne2025_GWI_Results/ANNUAL_ESM1-2-LR/"+
                "GWI_results_ANNUAL_HISTORICAL-ONLY_SCENARIO--SMILE_ESM-"+exp_attr[2]+"_ENSEMBLE-MEMBER--"+
                                           str(model_run)+"_VARIABLES--GHG-Nat-OHF___REGRESSED-YEARS--1850-1950_to_1850-2100.csv", header=[0, 1])
            biasdict = {"SSP126": -0.052410192 - 0.0122/2, "SSP245": -0.057050231 -0.0258/2,"SSP370": -0.071617279 -0.0478/2 }
            curbias = biasdict[exp_attr[2]]

        elif (exp_attr[1]=='NorESM'):
            model_run_noresm = gen_orig_number(model_run,60)
            gwi_levels_curr0 = pd.read_csv(cur_path+"/Thorne2025_GWI_Results/ANNUAL_NorESM/"+
                "GWI_results_ANNUAL_HISTORICAL-ONLY_SCENARIO--NorESM_rcp45-"+exp_attr[3]+"_ENSEMBLE-MEMBER--"+
                                           str(model_run_noresm)+"_VARIABLES--GHG-Nat-OHF___REGRESSED-YEARS--1850-1950_to_1850-2099.csv", header=[0, 1])
            biasdict = {"Volc": -0.208317627207827/2 -0.18528/2 , "VolcConst": -0.250277422257974/2 -0.23926/2}
            curbias = biasdict[exp_attr[3]]
    gwi_levels_curr =gwi_levels_curr0.iloc[1:,]
    gwi_c =gwi_levels_curr['Ant'].to_numpy() - curbias
    lyearc = np.shape(gwi_c)[0]
    lyear = lyearc + syear
    
    rdists = []
    cdists = []
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4
    
    for i in range(0, lyearc):
        cdists.append(CustomSplineDistribution(ordvalues=gwi_c[i,ord_ind], a=-4, b=7))#these start in 1950

    def empirical_mean(year_idx,k):
        #print(f"emp_mean start {time() - start:.2f}s")

        means = np.full(np.shape(year_idx),np.nan)
        for yr in year_idx:
            if (yr >=lyear or yr <1850):
                means[yr - 1850] =  np.nan
            elif (k==0):
                if (yr >=syear):
                    means[yr - 1850] = cdists[yr-syear].mean
                else:
                    means[yr - 1850] = np.nan

        #print(f"emp_mean finished {time() - start:.2f}s")

        return means

    def empirical_se(year_idx,k):
        #print(f"emp_se start {time() - start:.2f}s")
        ses = np.full(np.shape(year_idx),np.nan)
        for yr in year_idx:
            if (yr >=lyear or yr <1850):
                ses[yr - 1850] =  np.nan
            elif (k==0):
                if (yr >=syear):
                    ses[yr - 1850] = cdists[yr-syear].std
                else:
                    ses[yr - 1850] = np.nan

            if (ses[yr - 1850] < temps_1std[yr - 1850]/min_fact):
                ses[yr - 1850] = temps_1std[yr - 1850]/min_fact
        #print(f"emp_se finished {time() - start:.2f}s")
        return ses

    def empirical_pvalue(year_idx, point,k, two_sided=True):
        year_idx = np.atleast_1d(year_idx)
        point = np.atleast_1d(point) #should be the same dimension
        # Initialize an array to store empirical p-values
        empirical_p = np.full(np.shape(point), np.nan)
        for i, yr in enumerate(year_idx):
            #i = yr-1850
            if (yr >=lyear or yr <1850):
                empirical_p[i] =  np.nan
            elif (k==0):
                if (yr >=syear):
                    try:
                        cdfpt = cdists[yr-syear].cdfn0(point[i])
                    except:
                        breakpoint()
                    if two_sided:
                        empirical_p[i] = 2*min(1 - cdfpt,cdfpt)
                    else:
                        empirical_p[i] =1-cdfpt

                else:
                    empirical_p[i] =  np.nan

        return empirical_p

    def empirical_log_likelihood(year_idx, point,k):
        year_idx = np.atleast_1d(year_idx)
        point = np.atleast_1d(point)
        # Initialize an array to store empirical p-values
        empirical_l = np.full(np.shape(point), np.nan)
        for i, yr in enumerate(year_idx):
            #i = yr-1850
            if (yr >=lyear or yr <1850):
                empirical_l[i] =  np.nan
            elif (k==0):
                if (yr >=syear):
                    empirical_l[i] = (cdists[yr-syear].pdfn0(point[i]))
        if((empirical_l==0).any()):
            breakpoint()
        return np.log(empirical_l )

    def ppf_resample(year_idx, n_samples,k):
        year_idx = np.atleast_1d(year_idx)
        resamps = np.full((len(year_idx),n_samples), np.nan)
        for i, yr in enumerate(year_idx):
            if (yr >=lyear or yr <1850):
                resamps[i] =  np.nan
            elif (k==0):
                if (yr >=syear):
                    uniform_samples = np.random.uniform(0, 1, size=n_samples)
                    resamps[i] = (cdists[yr-syear].ppf(uniform_samples))
        return resamps

    return {
        'mean': empirical_mean,
        'se': empirical_se,
        'pvalue': empirical_pvalue,
        'log_likelihood': empirical_log_likelihood,
        'resample': ppf_resample,

    }

    #return means, ses, empser.copy(), empser.copy()

