###ADAPTED FROM
#https://gitlab.com/saidqasmi/kcc_notebook/-/blob/master/KCC_notebook.ipynb

# Loading KCC package
library(KCC)
# and other useful package(s)
library(abind)
library(reticulate)
# Set random number generator to ensure reproducibility
set.seed(1)

# Sample size to derive normal distributions
Nres = 1000
sample_str = c("be", paste0("nres",1:Nres))
setwd("/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats/Methods/43_Forcing_Based/3_Human_Induced/KCC/kcc_notebook/")
year = 1850:2100
ny = length(year)

# Load forcings from the Priestley Center
load("FF_CMIP6.rda")
# Plot if needed
#plot(year,FF$FF$nat[FF$FF$year %in% year],type="l")

# Load EBM parameters fitted on available CMIP6 models
load("ebm_params.rda")

#ebm_params

# Compute the EBM response e, to be used in equation 7 in Qasmi and Ribes (2021)
e = ebm_response(FF,ebm_params,year,Nres)

# Plot the best-estimate of e
#plot(year,e[,"be"],type="l", ylab="GMST response to natural forcings")

# Load historical+ssp585 GMST ensemble means from CMIP6 models
X_glo_scen_tmp = loadRData("GMST_CMIP6_histssp585_ann_idx.Rdata")

# Extract the scenario and the chosen simulated years
X_glo_scen_full = X_glo_scen_tmp[as.character(year),"histssp585",,]

# Consider the same models in historical+ssp585 and piControl simulations
models_scen_full = dimnames(X_glo_scen_full)$model
# Check if any NA in times series
models = models_scen_full[apply(is.na(X_glo_scen_full[,,"value"]), 2, sum) == 0] 
# The following CMIP6 models will be used to derive the prior

Nmod = length(models)
X_glo_scen = X_glo_scen_full[,models,]
dec_glo = T #

if (dec_glo == T) {
  # Extract the forced response to ALL and NAT forcings and fit to a Gaussian distribution
  X_fit_glo = array(NA, dim = c(ny, 2, Nmod),
                    dimnames = list(year = as.character(year),
                                    forcing = c("all","nat"),
                                    model = models))
  # degrees of freedom for the spline function
  df = 6                            
  message("Decomposition for each model (take some time!)...")
  X_fit_glo[,c("all","nat"),] = x_fit(X_glo_scen[,,"value"], e, df, ant=F)
} else {
  X_fit_glo = loadRData("X_fit_GMST_CMIP6_histssp585_ann.Rdata")
}

year_obsNew = 1850:2024
year_obsO = 1850:2019
Xo_glo_full = read.csv("HadCRUT.5.0.2.0.analysis.ensemble_series.global.annual.csv",row.names = 1,)
Xo_gloNew = Xo_glo_full[-c(1,2)]
colnames(Xo_gloNew)=1:200

#REPEAT from this section down
styr=51 #101
for (i in styr: length(year_obsNew)) {
end_idx = 1850+ i
print(end_idx)
year_obsN = 1850:(end_idx-1)
if (end_idx>1990){
    ref_obs = 1961:1990 }
else if (end_idx>1960){
    ref_obs = 1931:1960 }
else if (end_idx>1930){
  ref_obs = 1901:1930 }
else if (end_idx>1900){
  ref_obs = 1871:1900 }
else {ref_obs = 1850:1879 }


ny_obs = length(year_obsN)
Xo_glo_fullCW = loadRData("GMST_CW_ann.Rdata")

Xo_glo <- data.matrix(Xo_gloNew[1:i,], rownames.force = NA)

#RETAIN OLD DATASET TO COMPARE
Xo_glo_CW = Xo_glo_fullCW[as.character(year_obsO),]

# Check the ensemble size: 170 years, 100 members
#dim(Xo_glo)

# Estimate the response to all external forcings by the multimodel ensemble mean
raw_mmm = apply(X_fit_glo[,"all",], 1, mean, na.rm=T)
# Calculate obs residuals, our estimate of internal variability
raw_mmm_c = raw_mmm - mean(raw_mmm[year %in% ref_obs]) # raw_mmm_c must be anomalies wrt 1961-1990
Xo_glo_med = apply(Xo_glo, 1, median)
Xo_glo_c = Xo_glo_med - mean(Xo_glo_med[year_obsN %in% ref_obs])
res_glo = Xo_glo_c - raw_mmm_c[as.character(year_obsN)]
# Fit the parameters of the MAR models on residuals
message("Fitting MAR parameters (may take some time!)...")
theta_obs_glo = estim_mar2_link(res_glo)
# Compute the associated covariance matrix
Sigma2_obs_iv_glo = Sigma_mar2(theta_obs_glo,res_glo)

if(FALSE){ #not plotting this each time
acf_res = acf(res_glo, plot=F)
par(cex=1.1,font=2, font.lab=2,font.axis=2,lwd=2,mgp=c(2.8,.7,0),mar=c(4,4,1,4),las=2,tcl=-.4,cex.lab=1.2)
plot(acf_res,ylab="ACF of GMST residuals wrt CMIP6 multimodel mean",ci=0,lwd=2)
colors_acf = c("red", "black")
lines(0:20, Sigma2_obs_iv_glo[1,1:21]/Sigma2_obs_iv_glo[1,1], col=colors_acf[1], lwd=2)
legend("topright", legend = c("MAR fit", "OBS"), col = colors_acf, lwd=2, lty=1)}

Sigma2_obs_glo = Sigma2_obs_iv_glo + var(t(Xo_glo))

year = as.integer(dimnames(X_fit_glo)$year)

#print(year_S_str)
names(dimnames(Xo_glo))<-names(dimnames(Xo_glo_CW))
print(names(dimnames(Xo_glo)))


X_krig_glo_list = prior2posterior(X_fit_glo, Xo_glo, Sigma2_obs_glo, Nres, centering_CX=T, ref_CX=1850:1900)

# Convert output lists to arrays
# X_unconstrained
X_uncons_glo = mvgauss_to_Xarray(X_krig_glo_list$uncons$mean,X_krig_glo_list$uncons$var,Nres)
# X_constrained
X_cons_glo = mvgauss_to_Xarray(X_krig_glo_list$cons$mean,X_krig_glo_list$cons$var,Nres)
X_cons_glo_all = (X_cons_glo[,,"all"])
write.csv(X_cons_glo_all, paste0("kcc_historical_all_",end_idx,".csv"), row.names=FALSE)
X_cons_glo_nat = (X_cons_glo[,,"nat"])
write.csv(X_cons_glo_nat, paste0("kcc_historical_nat_",end_idx,".csv"), row.names=FALSE)
}

