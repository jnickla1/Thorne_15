import pandas as pd
import numpy as np
from scipy.stats import norm, gaussian_kde
from scipy.optimize import minimize
from numpy.polynomial.hermite import hermgauss
from custom_stack_utils import predict_non_nan, infill_log_likelihoods, fit_stacking_weights_logp
import pdb;
import sys

if __name__ == '__main__':
    comparison_type = 'ff'

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
        all_methodsh=all_methods.copy()
        err_vars = df_results['err_var'].to_numpy()
        err_vars100 = df_results['err_var100'].to_numpy()
        best_scales = df_results['best_alter_scale'].to_numpy()
        central_est = np.load(f'Results2/{outputfilename}_central_est.npy')
        node_lls = np.load(f'Results2/{outputfilename}_hermguass_lls.npy')
        newsamples = np.load(f'Results2/{outputfilename}_newsamples.npy')
        nyears = central_est.shape[1]
        
    if comparison_type[1]=='f':
        scenarios=["fut_ESM1-2-LR_SSP126_constVolc","fut_ESM1-2-LR_SSP245_constVolc","fut_ESM1-2-LR_SSP370_constVolc",
                   "fut_NorESM_RCP45_Volc","fut_NorESM_RCP45_VolcConst"]
        nmethods=55
        nyears=250
        const_shape = (5,60,nmethods)
        tseries_shape = (5,60,nmethods,250)
        sample_shape = (5,60,nmethods,250,nsamples)
        err_varsf    = np.full(const_shape,np.nan)
        err_vars100f = np.full(const_shape,np.nan)
        best_scalesf = np.full(const_shape,np.nan)
        central_estfut = np.full(tseries_shape,np.nan)
        node_llsf = np.full((5,60,nmethods,nyears,nnodes),np.nan)

        standardsf = np.full((5,60,nyears),np.nan)
        standards_sef = np.full((5,60,nyears),np.nan)
        
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
                    central_estfut[exp_index,model_run,:,:] = np.load(f'Results2/{outputfilename}_central_est.npy',mmap_mode="r")[:,0:250]
                    node_llsf[exp_index,model_run,:,:,:] = np.load(f'Results2/{outputfilename}_hermguass_lls.npy',mmap_mode="r")[:,0:250,:]
                    standardsf[exp_index,model_run,:] = np.load(f'Results2/{outputfilename}_standard.npy')[0:250]
                    standards_sef[exp_index,model_run,:] = np.load(f'Results2/{outputfilename}_standard_se.npy')[0:250]
        

        if comparison_type[0]=='f':
            err_vars=err_varsf
            err_vars100=err_vars100f
            node_lls=node_llsf
            #in the ff case our basis to construct the weights is the future runs too

#subset the read-in tables to set list of methods
#            only 8 vs ~16 vs ~30 methods
#            do this in a big outer loop
    nmethods_list = [int(sys.argv[1])] #8 16 later 30
    years = np.arange(1850,1850+nyears)
    for nnmethods in nmethods_list:
        #big outer loop
        mask_a=[]
        avail_methods_list = []
        if nnmethods ==10 or nnmethods == 8:
            # --- a) Subset to just the 8 selected methods ---
            avail_methods_list = [ "CGWL10y_for_halfU","EBMKF_ta2","EBMKF_ta4","GAM_AR1",
                 "lowess1dg20wnc","Kal_flexLin","FaIR_comb_unB","GWI_tot_CGWL","GWI_tot_SR15","CGWL10y_sfUKCP"]
            mask_a = all_methods.isin(avail_methods_list).values
            avail_methods = all_methods[mask_a]
            #df_cropped = df[mask_a].reset_index(drop=True)
            

            
        elif nnmethods ==30:
            rmse_df = pd.read_csv("current_methods_statistics_250818True.csv")
            bad_rmse_methods = rmse_df.loc[rmse_df['RMS'] > 0.06, 'method_name']
            #manual_exclude = ['CGWL10y_forec', 'some_other_method']  # replace with actual names
            manual_exclude = ['']
            exclude_methods = set(bad_rmse_methods).union(manual_exclude)
            mask_a = ~all_methods.isin(exclude_methods).values
            avail_methods = all_methods[mask_a]


        elif nnmethods ==16:
            rmse_df = pd.read_csv("current_methods_statistics_250818True.csv")
            bad_rmse_methods = rmse_df.loc[rmse_df['RMS'] > 0.06, 'method_name']
            #manual_exclude = ['CGWL10y_forec', 'some_other_method']  # replace with actual names
            manual_exclude = ['']
            df126 = pd.read_csv("averaged_runs126.csv")
            bad_126 = df126.loc[df126['100RMS'] > 0.075, 'method_name']
            df245=pd.read_csv("averaged_runs245.csv")
            bad_245 = df245.loc[df245['100RMS'] > 0.075, 'method_name']
            df370=pd.read_csv("averaged_runs370.csv")
            bad_370 = df370.loc[df370['100RMS'] > 0.075, 'method_name']
            dfvolc = pd.read_csv("averaged_runsVolc.csv")
            bad_volc = dfvolc.loc[ np.logical_and(dfvolc['100RMS'] > 0.085, dfvolc['e100RMS'] > 0.085), 'method_name']
            exclude_methods = set(bad_rmse_methods).union(manual_exclude).union(bad_126).union(bad_245).union(bad_370).union(bad_volc)
            mask_a = ~all_methods.isin(exclude_methods).values
            avail_methods = all_methods[mask_a]
            #print(avail_methods)
        print(str(len(avail_methods))+" methods")
        nmethods=len(avail_methods)
##        if comparison_type=='hh':
##            vars100 = err_vars100[mask_a]
##            #scales =  best_scales[mask_a] #dont need from this point on to use the scales
##            lls = node_lls[mask_a] #np.shape(lls) =(10, 175, 100)
##            centralsh = central_est[mask_a]
##            samplesh = newsamples[mask_a]
##
##        elif comparison_type=='hf':
##            mask_b = all_methodsh.isin(avail_methods).values
##            vars100 = err_vars100[mask_b]
##            lls = node_lls[mask_b]
##            centralsh = central_est[mask_b]
##            samplesh = newsamples[mask_b]
##            centralsf = central_estfut[:,:,mask_a,:]

        nruns = 60 #SETTING UP THE EVALUATION LOOP
        endyear=2099
        summethods_shape = (5,60,2)
        #ncrosses_array=np.full(summethods_shape,np.nan)
        #firstcross15_diff=np.full(summethods_shape,np.nan)
        #rmse_array= np.full(summethods_shape,np.nan)
        #kl_array= np.full(summethods_shape,np.nan)
        sum2methods_shape = (5,2)
        #firstcross15_sum=np.full(sum2methods_shape,np.nan)
        lhund=-100
        for exp_index in [int(sys.argv[2])]:
            experiment_type = scenarios[exp_index]
            for model_run in [int(sys.argv[3])]:
                if(experiment_type[0:8]=="fut_ESM1" and model_run>=50):
                    continue
                err_vars100a=err_vars100[:,:,mask_a].copy()
            #    err_vars100a[exp_index,model_run,:]=np.nan #masking out this run's err_vars
                vars100 = np.nanmean(err_vars100a,axis=(0,1))
                #combine all to one variance estimate
                #scales =  best_scales[:,:,mask_a] #DONT NEED to combine the scales
                print("setting up lls")
                lls0 = node_lls[:,:,mask_a,lhund:-10,:].copy() 
             #   lls0[exp_index,model_run,:,:]=np.nan
                lls_wnan = lls0.transpose(2, 0, 1, 3, 4).reshape(nmethods, (-10-lhund)* 5 * 60, nnodes)
                finite = np.isfinite(lls_wnan)                           # (K, T, M)
                valid_T = np.all(finite, axis=(0, 2))               # keep timepoints with no NaNs across models & nodes
                lls = lls_wnan[:, valid_T, :]                      # (K, T_keep, M)
                #TODO - for ff drop the particular one we are leaving out here
                central_estfuta=central_estfut[:,:,mask_a,lhund:-10].copy()
                central_estfuta[exp_index,model_run,:]=np.nan
                centralsf0 = central_estfuta.transpose(2, 0, 1, 3).reshape(nmethods, (-10-lhund) * 5 * 60)
                centralsf = centralsf0[:, valid_T] 
                print("created lls") 
        #generate weights/distributions for each: IVW and stacking
                ivweights = 1/vars100
                x, w = hermgauss(100)
                z = np.sqrt(2.0) * x
                alpha = w / np.sqrt(np.pi)
                infilled_lls = infill_log_likelihoods(lls,penalty=0)
                stack_weights = fit_stacking_weights_logp(infilled_lls,alpha)

                print("weights fitted")
                

                if comparison_type[0]=='f':
                    #ivcenters = np.nansum(central_estfut[:,:,mask_a,lhund:-10] * ivweights[None,None,:,None],axis=2) / np.nansum( ~np.isnan(central_estfut[:,:,mask_a,lhund:-10]) * ivweights[None,None,:,None],axis=2)
                    valid = np.any(~np.isnan(central_estfut[:, :, mask_a, lhund:-10]), axis=(2))
                    num = np.nansum(central_estfut[:, :, mask_a, lhund:-10] * ivweights[None, None, :, None], axis=2)
                    den = np.nansum((~np.isnan(central_estfut[:, :, mask_a, lhund:-10])) * ivweights[None, None, :, None], axis=2)
                    ivcenters = np.full(num.shape, np.nan)
                    np.divide(num, den, out=ivcenters, where=valid)
                    #these constructed for all
                    blended_trial = np.full((10,-10-lhund,nsamples), np.nan)
                    standards_trial = np.full((10,-10-lhund), np.nan)
                    standards_se_trial = np.full((10,-10-lhund), np.nan)
                    for exp_indext, experiment_typet in enumerate(scenarios):
                        for model_runt in range(2):
                            model_runt2 = model_runt*2
                            if (exp_index == exp_indext and model_run==model_runt2):
                                model_runt2 = model_runt+1
                            else:
                                model_runt2 = model_runt
                            outputfilenamet = f'{experiment_typet}{model_runt2}' 
                            tempsamplesf = np.load(f'Results2/{outputfilenamet}_newsamples.npy',mmap_mode="r")[:,lhund:-10,:] #newsamplesfut[exp_index,model_run,:,:,:] =
                            blended_trial[(exp_indext*2 + model_runt),:,:] = predict_non_nan(stack_weights, tempsamplesf[mask_a,:,:])
                            standards_trial[(exp_indext*2 + model_runt),:] = standardsf[exp_indext,model_runt2,lhund:-10]
                            standards_se_trial[(exp_indext*2 + model_runt),:] = standards_sef[exp_indext,model_runt2,lhund:-10]

        #training standard to fix variances
                if True:
                    #def rescale_log_likelihood(scale_alt):
                    #    log_lik=stats.norm.logpdf(standard[-100:],loc=ivcenters[-100:],scale=scale_alt)
                    #    return -np.nansum(log_lik)
                    def rescale_KL_anystd(scale,standardi,standard_sei,eeoffset=0,ivcentersi=ivcenters):
                        yearKL = np.log(scale/standard_sei) + (standard_sei**2 + (standardi-ivcentersi)**2 )/2/(scale)**2 - 0.5
                        return np.nansum(yearKL,axis=None)
                    
                    def rescale_KL(scale): #now this is a passthrough function
                        return rescale_KL_anystd(scale,standardsf[:,:,lhund:-10],standards_sef[:,:,lhund:-10])
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

                    def kl_total_true_vs_empirical_gh(sharp, stand, stand_se, eps=1e-300,eeoffset=0):
                        """
                        Sum_y KL( P_y || Q_y ),  P_y = N(standard[y], standard_se[y]^2),
                        Q_y = empirical via scipy.stats.gaussian_kde on sharp[y, :], using Gauss–Hermite quadrature.
                        """
                        total = 0.0
                        for y in range(len(stand)+lhund,len(stand)-10-eeoffset):
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
                        sharpened = sharpen_samples(blended_trial.reshape(-1, nsamples), scale[0], seed = 42)
                        return kl_total_true_vs_empirical_gh(sharpened, standards_trial.reshape(-1), standards_se_trial.reshape(-1))

                    best_nsharp = minimize(rescale_KL_stack, x0=2, bounds=[(1, 10)], tol=0.005).x[0]
                    #sharpened_blended = sharpen_samples(blended_trial.reshape(-1, nsamples), best_nsharp, seed = 2)
                    #best_stack_KL = kl_total_true_vs_empirical_gh(sharpened_blended, standards_trial.reshape(-1), standards_se_trial.reshape(-1))
                    #not needed    
                print("finished fitting weights and variances")
            #fut standard can read in from "Results/ensemble_mean_fut_ESM1-2-LR_SSP370_constVolc.csv"
        #test: compare results to standard (loop thorugh all experiment runs)

                # Normalize ivweights
                iv_row = ivweights / np.sum(ivweights)

                # Build dataframe
                df = pd.DataFrame(
                [iv_row, stack_weights],columns=avail_methods)
                df.insert(0, "optimized var/nsamps", [best_scale_iv, best_nsharp])

                # Save to CSV
                df.to_csv(f"ff_{nnmethods}_weights_summary.csv", index=False)
        
