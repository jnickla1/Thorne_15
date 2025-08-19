import pandas as pd
import numpy as np
from scipy.stats import norm, gaussian_kde
from scipy.optimize import minimize
from numpy.polynomial.hermite import hermgauss
from custom_stack_utils import predict_non_nan, infill_log_likelihoods, fit_stacking_weights_logp
import pdb;
#import sys

if __name__ == '__main__':
    comparison_type = 'hh' #sys.argv[1] 'hf' or 'ff'

#recreate tables: read in files and combine
# 3 options: historical (all) -> historical NOTE: not purely retrospective by weights and scalings
#            historical -> future
#            future (loo) and historical -> future
#this part specified as sys argument
    nsamples=1000
    nnodes = 100
        
    if comparison_type[0]=='h':
        outputfilename = 'historical'
        df_results = pd.read_csv(f'Results2/{outputfilename}_names_var_scale.csv', index_col=0)
        all_methods = df_results['method_name']
        err_vars = df_results['err_var'].to_numpy()
        err_vars100 = df_results['err_var100'].to_numpy()
        best_scales = df_results['best_alter_scale'].to_numpy()
        central_est = np.load(f'Results2/{outputfilename}_central_est.npy')
        node_lls = np.load(f'Results2/{outputfilename}_hermguass_lls.npy')
        newsamples = np.load(f'Results2/{outputfilename}_newsamples.npy')
        nyears = central_est.shape[1]
        years = np.arange(1850,1850+nyears)
        
    if comparison_type[1]=='f':
        scenarios=["fut_ESM1-2-LR_SSP126_constVolc","fut_ESM1-2-LR_SSP245_constVolc","fut_ESM1-2-LR_SSP370_constVolc",
                   "fut_NorESM_RCP45_Volc","fut_NorESM_RCP45_VolcConst"]
        nmethods=55
        nyears=250
        const_shape = (5,60,nmethods)
        tseries_shape = (5,60,nmethods,250)
        sample_shape = (5,60,nmethods,250,nsamples)
        err_varsf    = np.fill(const_shape,np.nan)
        err_vars100f = np.fill(const_shape,np.nan)
        best_scalesf = np.fill(const_shape,np.nan)
        central_estfut = np.fill(tseries_shape,np.nan)
        node_llsf = np.fill((5,60,nmethods,nnodes),np.nan)
        newsamplesfut =  np.fill(sample_shape,np.nan)
        
        for exp_index, experiment_type in enumerate(scenarios):
            for model_run in range(60):
                if(experiment_type[0:8]=="fut_ESM1" and model_run>=50):
                    continue
                    #fill everything with nans
                else:
                    outputfilename = f'{experiment_type}{model_run}'
                    df_results = pd.read_csv(f'Results2/{outputfilename}_names_var_scale.csv', index_col=0)
                    if(exp_index==0 and model_run==0):
                        all_methods = df_results['method_name'] #only need this once
            
                    err_varsf[exp_index,model_run,:] = df_results['err_var'].to_numpy()
                    err_vars100f[exp_index,model_run,:] = df_results['err_var100'].to_numpy()
                    best_scalesf[exp_index,model_run,:] = df_results['best_alter_scale'].to_numpy()
                    central_estfut[exp_index,model_run,:,:] = np.load(f'Results2/{outputfilename}_central_est.npy')
                    node_llsf[exp_index,model_run,:,:,:] = np.load(f'Results2/{outputfilename}_hermgauss_lls.npy')
                    newsamplesfut[exp_index,model_run,:,:,:] = np.load(f'Results2/{outputfilename}_newsamples.npy')
        
                #if(experiment_type[0:8]=="fut_NorE"):
                    #append nans to everything along last year dimension to make every array the same size

        if comparison_type[0]=='f':
            err_vars=err_varsf
            err_vars100=err_vars100f
            node_lls=node_llsf
            #in the ff case our basis to construct the weights is the future runs too

#subset the read-in tables to set list of methods
#            only 8 vs ~16 vs ~30 methods
#            do this in a big outer loop
    nmethods_list = [10] #8 16 later 30

    for nmethods in nmethods_list:
        #big outer loop
        mask_a=[]
        avail_methods_list = []
        if nmethods ==10:
            # --- a) Subset to just the 8 selected methods ---
            avail_methods_list = [ "CGWL10y_for_halfU","EBMKF_ta2","EBMKF_ta4","GAM_AR1",
                 "lowess1dg20wnc","Kal_flexLin","FaIR_comb_unB","FaIR_nonat_unB","GWI_anthro_SR15","CGWL10y_sfUKCP"]
            mask_a = all_methods.isin(avail_methods_list).values
            avail_methods = all_methods[mask_a]
            #df_cropped = df[mask_a].reset_index(drop=True)
            

            
        elif nmethods ==30:
            rmse_df = pd.read_csv("current_methods_statistics_250818True.csv")
            bad_rmse_methods = rmse_df.loc[rmse_df['RMSE'] > 0.06, 'method_name']
            manual_exclude = ['CGWL10y_forec', 'some_other_method']  # replace with actual names
            exclude_methods = set(bad_rmse_methods).union(manual_exclude)
            mask_a = ~all_methods.isin(exclude_methods).values
            avail_methods = all_methods[mask_a]


        # --- c) Further filtering (you fill in) ---
        # mask_c = ...


        #TODO: PROBLEM WITH hf CASE referencing a different list of methods
        
        if comparison_type[0]=='h':
            vars100 = err_vars100[mask_a]
            #scales =  best_scales[mask_a] #dont need from this point on to use the scales
            lls = node_lls[mask_a] #np.shape(lls) =(10, 175, 100)
            
        elif comparison_type[0]=='f':
            vars100 = np.nanmean(err_vars100[:,:,mask_a],axis=(0,1))
            #combine all to one variance estimate
            #scales =  best_scales[:,:,mask_a] #DO WE NEED to combine the scales?
            lls0 = node_lls[:,:,mask_a,:,:] #np.shape(lls) =(10, 175, 100)
            lls = lls0.transpose(2, 0, 1, 3, 4).reshape(nmethods, nyears * 5 * 60, nnodes)
            #TODO - for ff drop the particular one we are leaving out here

        if comparison_type[0]=='h':
            centralsh = central_est[mask_a]
            samplesh = newsamples[mask_a]

        if comparison_type[1]=='f':
            centralsf = central_estfut[:,:,mask_a,:]
            samplesf = newsamplesfut[:,:,mask_a,:,:]
            #note that these two have different dimesions and that's ok
            #both might be needed to calibrate the final se adjustment, then verify accuracy
    
#generate weights/distributions for each: IVW and stacking
        ivweights = 1/vars100
        x, w = hermgauss(100)
        z = np.sqrt(2.0) * x
        alpha = w / np.sqrt(np.pi)
        infilled_lls = infill_log_likelihoods(lls[:,10:-10,:],penalty=0)
        stack_weights = fit_stacking_weights_logp(infilled_lls,alpha)
        
#new center and dist generation function called once or many times for LOO
        if comparison_type[0]=='h':
            ivcenters = np.nansum(centralsh * ivweights[:,None],axis=0) / np.nansum( ~np.isnan(centralsh) * ivweights[:,None],axis=0)
            blended = predict_non_nan(stack_weights, samplesh)

        if comparison_type[0]=='f':
            ivcenters = np.nansum(centralsf * ivweights[None,None,:,None],axis=2) / np.nansum( ~np.isnan(centralsf) * ivweights[None,None,:,None],axis=2)
            blended = np.full((5,60,250,nsamples), np.nan)
            for exp_index, experiment_type in enumerate(scenarios):
                for model_run in range(60):
                    if(experiment_type[0:8]=="fut_ESM1" and model_run>=50):
                        continue
                    #TODO: also continue if we hit the LOO one that 
                    blended[exp_index,model_run,:,:] = predict_non_nan(stack_weights, samplesf[exp_index,model_run,:,:,:])


#training standard to fix variances
        if comparison_type[0]=='h':
            standard = np.load(f'Results2/{outputfilename}_standard.npy')
            standard_se = np.load(f'Results2/{outputfilename}_standard_se.npy')
            #def rescale_log_likelihood(scale_alt):
            #    log_lik=stats.norm.logpdf(standard[-100:],loc=ivcenters[-100:],scale=scale_alt)
            #    return -np.nansum(log_lik)
            def rescale_KL(scalearr):
                scale=scalearr[0]
                yearKL = np.log(scale/standard_se[-100:-10]) + (standard_se[-100:-10]**2 + (standard[-100:-10]-ivcenters[-100:-10])**2 )/2/(scale)**2 - 0.5
                return np.nansum(yearKL)
            best_scale_iv = minimize(rescale_KL, x0=0.04, bounds=[(0.001, 1)]).x[0] #gets a constant uncertainty
            best_iv_KL = rescale_KL(best_scale_iv)

            def sharpen_samples(blended_in: np.ndarray, nsharp: float, seed: int ) -> np.ndarray:
                """Resample-with-replacement means (with fractional last draw) per year."""
                nyears, nsamples = blended_in.shape
                k = int(np.floor(nsharp))
                frac = float(nsharp - k)
                draws = k + (1 if frac > 0 else 0)
                rng = np.random.default_rng(seed)
                out = np.empty_like(blended_in)
                for y in range(nyears):
                    idx = rng.integers(0, nsamples, size=(nsamples, draws))
                    picks = blended_in[y, idx]                         # (nsamples, draws)
                    if frac > 0:
                        num = picks[:, :k].sum(axis=1) + frac * picks[:, -1]
                        den = k + frac
                    else:
                        num = picks.sum(axis=1)
                        den = k if k > 0 else 1.0
                    out[y] = num / den
                return out

            def kl_total_true_vs_empirical_gh(sharp, stand, stand_se, eps=1e-300):
                """
                Sum_y KL( P_y || Q_y ),  P_y = N(standard[y], standard_se[y]^2),
                Q_y = empirical via scipy.stats.gaussian_kde on sharp[y, :], using Gauss–Hermite quadrature.
                """
                total = 0.0
                for y in range(len(stand)-100,len(stand)-10):
                    if np.isnan(stand[y]):
                        continue
                    mu, sd = float(stand[y]), float(max(stand_se[y], 1e-12))
                    y_nodes = mu + sd * z                         # transform nodes to P_y support
                    kde = gaussian_kde(sharp[y])                  # KDE for Q_y
                    log_p = norm.logpdf(y_nodes, loc=mu, scale=sd)
                    log_q = np.log(np.maximum(kde.evaluate(y_nodes), eps))
                    total += float(np.sum(alpha * (log_p - log_q)))
                return total

            
            def rescale_KL_stack(scale: float):
                sharpened = sharpen_samples(blended, scale[0], seed = 42)
                return kl_total_true_vs_empirical_gh(sharpened, standard, standard_se)

            best_nsharp = minimize(rescale_KL_stack, x0=2, bounds=[(1, 10)], tol=0.005).x[0]
            sharpened_blended = sharpen_samples(blended, best_nsharp, seed = 2)
            best_stack_KL = kl_total_true_vs_empirical_gh(sharpened_blended, standard, standard_se)
                
    #fut standard can read in from "Results/ensemble_mean_fut_ESM1-2-LR_SSP370_constVolc.csv"
#test: compare results to standard (loop thorugh all experiment runs)
            
    inum=12 #internal interpolation within years
    fineyrs_all = np.arange(years[0],years[-1]+1/inum,1/inum)
    std_intp0 = np.interp(fineyrs_all,years,standard)
    std_intp = std_intp0[~np.isnan(std_intp0)]
    fineyrs_c0 = fineyrs_all[~np.isnan(std_intp0)]
    fineyrs_c = fineyrs_c0[1:]
    thrshs= np.arange(0.5,np.nanmax(standard) ,.1) #thresholds #CHANGE FOR NONHISTORICAL
    print(f"evaluating {len(thrshs)} of 0.1°C thresholds, starting at 0.5°C")
    closest_years = [-1/inum/2+fineyrs_c[np.logical_and(std_intp[0:-1]<i, std_intp[1:]>=i)][0] for i in thrshs]
                   #will have a variable number of steps, at least including 0.5
    closest_yrs_rnd = np.round(closest_years)
    print(np.round(closest_years,decimals=2))
    win_sz=25

    aedyrs=np.zeros((2,len(thrshs)))
    maedyrs=np.zeros((2,len(thrshs)))
    faedyrs=np.zeros((2,len(thrshs)))
    ncrosses = np.zeros(2)
    for m in range(2):
        ncross=0
        for j in range(len(closest_years)):
            evalmin=int(closest_years[j])-win_sz
            evalmax=min(int(closest_years[j])+win_sz,2024) #CHANGE FOR NONHISTORICAL
            evalyrs = np.arange(evalmin,evalmax)
            fineevalyrs = np.arange(evalyrs[0],evalyrs[-1]+1/inum,1/inum) 
            this_method_p_steps = np.full(np.shape(evalyrs),np.nan)
            if m==1:
                #this_method_p_steps = result['pvalue'](evalyrs, np.full(np.shape(evalyrs),thrshs[j]),k, two_sided=False)
                this_method_p_steps=np.zeros(np.shape(evalyrs))
                for i,e in enumerate(evalyrs):
                    kde = gaussian_kde(sharpened_blended[e-1850])                  
                    this_method_p_steps[i] = kde.integrate_box_1d(thrshs[j],15) #integrate probability above threshold
            else:
                this_method_p_steps = norm.cdf(((ivcenters[evalyrs-1850]-thrshs[j])/ best_scale_iv))

            first_non_nan = np.argmax(~np.isnan(this_method_p_steps))
            this_method_p_steps[:first_non_nan] = 0
            # Replace NaNs at the end with 1s if they exist
            last_non_nan = len(this_method_p_steps) - np.argmax(~np.isnan(this_method_p_steps)[::-1]) - 1
            this_method_p_steps[last_non_nan + 1:] = 1
            psteps_intp = np.interp(fineevalyrs,evalyrs,this_method_p_steps)

            now_cross = np.logical_and(psteps_intp[0:-1]<0.5, psteps_intp[1:]>=0.5)
            ncross= ncross + np.sum(now_cross)
            if j==0 and np.sum(now_cross)>1 : #special tag at 0.5 #CHANGE FOR NONHISTORICAL
                ncross = ncross + 0.1*np.sum(now_cross)
            ncrosses[m]=ncross
            fineevalyrsh=fineevalyrs[0:-1]/2 + fineevalyrs[1:]/2
            diffcross = (fineevalyrsh[now_cross] - closest_years[j])
            evalmean = np.nanmean(diffcross) #if crossing multiple times take the mean
            if np.isnan(evalmean):
                evalmean=win_sz #just put 15 yrs difference
            aedyrs[m,j] = evalmean
            if len(diffcross)>0:
                mevalmean = diffcross[np.nanargmax(abs(diffcross))] #if crossing multiple times take the worst one
                feval = diffcross[0]
            else:
                mevalmean=15 #just put 15 yrs difference
                feval = 15
            maedyrs[m,j] = mevalmean
            faedyrs[m,j]=feval
    print("\nMethods: ")
    print(avail_methods)
    print("\nWeights: ")
    print(f"ivarw: {np.round(ivweights/np.sum(ivweights),decimals=2)}")
    print(f"stack: {np.round(stack_weights,decimals=2)}")
    print("All threshold RMSE in years:"+ str(np.sqrt(np.square(aedyrs).mean(axis=1))))
    print(f"0.5°C threshold error years:{faedyrs[:,0]}")
    print(f"1.0°C threshold error years:{faedyrs[:,5]}")
    #array of central_est, approx_se,
    sharp_blend_central = np.mean(sharpened_blended,axis=1)
    sharp_blend_se = np.std(sharpened_blended,axis=1)
    #list of RMSE
    lhund=-100
    ivarw_rmse = np.sqrt(np.nanmean((sharp_blend_central[lhund:-1]-standard[lhund:-1])**2))
    stack_rmse = np.sqrt(np.nanmean((ivcenters[lhund:-1]-standard[lhund:-1])**2))
    print(f"RMSEs: ivarw:{ivarw_rmse}, stack:{stack_rmse}") #yes there is some improvement for both but it's slight #store these RMSEs
    #compute probability of within 1 year:
    for m in range(2):
        kdeyear = gaussian_kde(maedyrs[m,:])
        ptot = kdeyear.integrate_box_1d(-1,1)
        print(f'{m} total p wihin 1 yr:{ptot}') #store these probabilities
    # KL divergence
    #best_iv_KL, best_stack_KL,
#save output tables: 
