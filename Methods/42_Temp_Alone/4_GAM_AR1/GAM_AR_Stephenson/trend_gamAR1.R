###########################################################################
# Trend estimate using GAM+AR1 noise courtesy of Simon Wood's
# routine gamAR1 at end of script.  

#####################################################################
gamAR1 <- function(G,rho,method="ML") {
  ## G is a pre-fit object from `gam' (model set up, but
  ## not yet estimated). This routine manipulates this object
  ## to account for AR1 errors, and then estimates model.
  y0 <- G$y;X0 <- G$X;n <- G$n
  ## coeffs of inverse correlation matrix, Ci
  ## diag(Ci) = c(a,b,b,b,....,b,b,a), while the
  ## sub &  super diagonals are vectors containing `c'
  a <- 1/(1-rho^2);b <- 1 + 2*rho^2*a
  c <- - rho*a;
  ## Coeffs of the choleski factor, R, of Ci, s.t. R'R=C
  ## leading diagonal of R is e, super diagonal is f.
  f <- -sqrt(a-1);e <- c/f
  ## test...
  f <- f/e ; e <- e/e
  ## form Ry
  G$y <- c(y0[1:(n-1)]*e + y0[2:n]*f,y0[n])
  ## ... and RX
  G$X <- rbind(X0[1:(n-1),]*e + X0[2:n,]*f,X0[n,])
  
  b <- gam(G=G,method=method)
  
  if (method=="ML"||method=="REML") {
    # need to add 0.5 log |C| = log|R| to likelihood
    b$gcv.ubre <- b$gcv.ubre - n*log(e)
  }
  b
}

##########################################################


library(mgcv)

# Read in some yearly temperature data with variables Time and Anomaly
#dat=data.frame(read.table("~/Downloads/Thorne_figures/HadCRUT5.csv",header=TRUE,sep=","))
setwd("~/Downloads/Thorne_15_codefigurestats/Methods/42_Temp_Alone/4_GAM_AR1/GAM_AR_Stephenson/")
dat0=data.frame(read.table("~/Downloads/Thorne_15_codefigurestats/Common_Data/HadCRUT5.csv",header=TRUE,sep=","))
#colnames(dat0) <- c("Time","Anomaly","eCO2",	"opt_depth",	"R_Tvar",	"OHCA",	"ROC_tvar",	"anthro_clouds"	,"TSI")
par(mfrow=c(1,2))
styr=100
n_partitions = dim(dat0)[1]-styr+1
fits_matrix <- matrix(NA, nrow=n_partitions, ncol=nrow(dat0))
se_matrix <- matrix(NA, nrow=n_partitions, ncol=nrow(dat0))
fits_matrix0 <- matrix(NA, nrow=n_partitions, ncol=nrow(dat0))
se_matrix0 <- matrix(NA, nrow=n_partitions, ncol=nrow(dat0))

# Loop through each partition
for (i in styr:dim(dat0)[1]) {
  dat<- dat0[1:i, ]
# GAM fit assuming residuals are uncorrelated
G=gam(Anomaly~s(Time),data=dat)
pt=predict(G,dat,se=TRUE)

if(i==dim(dat0)[1]){
plot(dat$Time,dat$Anomaly,type="p",pch=20,col="gray",xlab="year",ylab="temp anomaly")
lines(dat$Time,pt$fit,lwd=2)
matplot(dat$Time,cbind(pt$fit-1.96*pt$se.fit,pt$fit+1.96*pt$se.fit), lwd=2,lty=c(1,1), col=c(2,2),type="l",add=TRUE)
title("GAM with AR1=0")
# Add 20-year running mean
lines(dat$Time,filter(dat$Anomaly,rep(1,20)/20),col="blue")
#write.csv(pt$fit, "gam_AR0_last_fit.csv", row.names=FALSE)
}
fits_matrix0[i-styr+1,1:i] <- pt$fit
se_matrix0[i-styr+1,1:i] <- pt$se.fit

# GAM fit that assumes residuals are AR(1) with autocorrelation estimated from lag-1 moment
r1=acf(residuals(G),plot=FALSE)$acf[2]

G=gam(Anomaly~s(Time),data=dat,fit=FALSE)
b <- gamAR1(G,rho=r1,method="ML")
pt=predict(b,dat,se=TRUE)
if(i==dim(dat0)[1]){
plot(dat$Time,dat$Anomaly,type="p",pch=20,col="gray",xlab="year",ylab="temp anomaly")
lines(dat$Time,pt$fit,lwd=2)
matplot(dat$Time,cbind(pt$fit-1.96*pt$se.fit,pt$fit+1.96*pt$se.fit), lwd=2,lty=c(1,1), col=c(2,2),type="l",add=TRUE)
title("GAM with estimated AR1")
lines(dat$Time,filter(dat$Anomaly,rep(1,20)/20),col="blue")}
fits_matrix[i-styr+1,1:i] <- pt$fit
se_matrix[i-styr+1,1:i] <- pt$se.fit
#write.csv(model.matrix(b), "gam_modelmatrix.csv", row.names=FALSE)
}

write.csv(fits_matrix, "gamAR1_fits_historical.csv", row.names=FALSE)
write.csv(se_matrix, "gamAR1_se_fits_historical.csv", row.names=FALSE)

write.csv(fits_matrix0, "gamAR0_fits_historical.csv", row.names=FALSE)
write.csv(se_matrix0, "gamAR0_se_fits_historical.csv", row.names=FALSE)


# Note: intervals are slightly wider and trend curve is slightly less bendy than for uncorrelated residuals case

# Add 20-year running mean









