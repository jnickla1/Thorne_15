## EM algorithm using Kalman Filtering, example 6.7 of Shumway and Stoffer
## Page 342-344

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import os
#regenerate the land surface and ocean surface estimates:
#T = 0.292*LandRec +0.708*OceanRec
# LandRec = (T - 0.708*OceanRec)/ 0.292

dir_path = os.path.dirname(os.path.realpath(__file__))
data = np.genfromtxt(open(dir_path+'/../../../Common_Data/toyKFmodelData8c.csv', "rb"),dtype=float, delimiter=',')
years0=data[:,0]
years0[0]=1850
temperature = data[:, 1]
n_iter=len(temperature) #174
#data_ocn = pd.read_csv(dir_path+'/../../../Common_Data/HadSST.4.0.1.0_annual_GLOBE.csv', sep=',')
#OceanRec = data_ocn['anomaly'].values



##moving_aves=np.empty(len(temperature))
##moving_aves[:] = np.nan
##N=4;cN=int(np.ceil(N/2));fN=int(np.floor(N/2));
##for i in range((fN),(len(temperature)-cN+1)):
##    lasta=i+cN;firsta=i-fN
##    moving_aves[i] = np.mean(temperature[firsta:lasta]);




##1. Initialize the procedure by selecting starting values for the parameters
u0= np.array([[-0.35],[0]])
sig0 = 0.2
phi = np.array([[1,1],[0,1]]) ## we are keeping this as a random walk, so NOT updating this
#Q =  np.array([[0.005,0.001],[0.001,0.0001]])
#R = np.array([[0.36,0.11],[0.11,0.16]])              #[0.25,0.006],[0.006,0.19]])
q00=1
##xt is 1d, yt is 2d
R000=0.36
R001=0.11
R011=0.16

A = np.array([[1,0]]) #emissions matrix


# allocate space for arrays





lYthetas = [-1000]
dellYtheta = 100

means = np.full(np.shape(years0),np.nan)
avg_len_l=10
avg_len_u=11 #actually 10 yrs below, that year, 10 yrs after
for i in range(avg_len_l, len(years0) - avg_len_u+1):
    chunk=temperature[i-avg_len_l:i+avg_len_u]
    means[i] = np.mean(chunk)

diff_lin = (means[:-2]+means[2:])/2 - means[1:-1]
diff_mean = (temperature - means)/np.sqrt(avg_len_l+avg_len_u)
Qo=np.cov(diff_mean[(avg_len_u+1):-(avg_len_u+1)],diff_lin[avg_len_u:-(avg_len_u)]) #keeping this


z1 = np.array([-0.02065589 , 5.94990778 , 0.23777725]) #precomputed

def run_method(yrs, temps, uncert, model_run, experiment_type):
    data_orig = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data_orig.loc[:,"Anomaly"].to_numpy()
    preind_base = np.mean(temps_obs[0:50])
    
    xhat, P, xhathat, Phathat, thisllike = KF_compute(yrs, temps+preind_base,np.zeros(np.shape(temps)), z1)
    return xhat-preind_base, P, xhathat-preind_base, Phathat


#while dellYtheta > 1:
def KF_compute(years, temperature, match_series, z0):
    
    n_iter = len(temperature)

    
    sz = (n_iter,2,1) # size of array
    sz2d=(n_iter,2,2) # state covariance is also 1d
    Q = Qo* z0[1] #constant scaling
    R = np.array([[ z0[2]]])
    
    xhat=np.zeros(sz)      # a posteri estimate of x
    P=np.zeros(sz2d)         # a posteri error estimate


    eta = np.zeros((n_iter,1,1))      # a posteri estimate of x

    xhatminus=np.zeros(sz) # a priori estimate of x
    Pminus=np.zeros(sz2d)    # a priori error estimate
    K=np.zeros((n_iter,2,1))         # gain or blending factor
    S=np.zeros((n_iter,1,1))

    xhathat=np.zeros(sz)   # smoothed a priori estimate of x
    Phat=np.zeros(sz2d)      # smoothed posteri error estimate 
    Khat=np.zeros(sz2d)      # smoothed gain or blending factor

    Phathat=np.zeros(sz2d)      # smoothed posteri error estimate 

#E-step: 6.1 (Kalman filter), 6.2 (Kalman Smoother), 6.3 (third covar est) for t=1,..n
    xhat[0] = u0
    P[0,:,:] = sig0
    Pminus[0,:,:] = P[0,:,:]
    thislYtheta=0
    thisllike=0
    
#forward loop:
    for k in range(1,n_iter):
        # time update
        
        xhatminus[k] = np.matmul(phi, xhat[k-1])
        
        Pminus[k] = np.matmul(np.matmul(phi,P[k-1]),np.transpose(phi))+Q   #np.matmul(np.matmul(F[k] ,P[k-1]), np.transpose(F[k])) +Q 

        S[k]= np.matmul(np.matmul(A ,Pminus[k]), np.transpose(A)) + R
        K[k] = np.matmul(Pminus[k],np.matmul(np.transpose(A),np.linalg.inv(S[k])))
        eta[k]=np.array([[temperature[k]]])-np.matmul(A,(xhatminus[k]))
        xhat[k] = xhatminus[k]+np.matmul(K[k],eta[k])
        P[k] = np.matmul((np.eye(2)- np.matmul(K[k],A) ),Pminus[k])
        
        thislYtheta = thislYtheta + 0.5 * np.log(np.linalg.det(S[k])) + np.matmul( np.matmul( np.transpose(eta[k]), np.linalg.inv(S[k])), eta[k])
        if(~np.isnan(match_series[k])):
            thisllike= thisllike+ stats.norm.logpdf(xhat[k][0],loc=match_series[k],scale=P[k][0,0])
            
    xhathat[n_iter-1]=xhat[n_iter-1]
    Phat[n_iter-1]=P[n_iter-1] 
    xhathat[0]=xhat[0]
    Phat[0]=P[0]
#backwards loop:
    for ik in range(2,n_iter+1):
        k=n_iter-ik
        try:
            Khat[k] = np.matmul(np.matmul(P[k],np.transpose(phi)),np.linalg.inv(Pminus[k+1])) #compute inverse for higher dimensions
        except:
            Khat[k] =Khat[k+1]
        xhathat[k] = xhat[k]+np.matmul(Khat[k],(xhathat[k+1]- np.matmul(phi, xhat[k])))
        Phat[k] = P[k] + np.matmul(np.matmul(Khat[k],(Phat[k+1]- Pminus[k])),np.transpose(Khat[k]))
        


    Phathat[n_iter-1]=  np.matmul(np.matmul( ( np.eye(2) - np.matmul( K[n_iter-1], A)), phi),P[n_iter-2])
    Phathat[0]=Phat[0]
#backwards loop to make Phathat:
    for ik in range(2,n_iter):
        k=n_iter-ik
        Phathat[k] = np.matmul( P[k],np.transpose(Khat[k-1])) + \
                np.matmul(np.matmul( Khat[k], (Phat[k] - np.matmul(phi,P[k]))),np.transpose(Khat[k-1]) ) #Eq. 6.56

#compute incomplete-data likelihood

    
    dellYtheta = thislYtheta[0,0] - lYthetas[-1]
    
    lYthetas.append(thislYtheta[0,0])

    #Used smoothed values to calculate S11, S10, S00 given 6.67-6.69 pg 
    S11 =np.array([[0,0],[0,0]])
    S10 =np.array([[0,0],[0,0]])
    S00 =np.array([[0,0],[0,0]])
    for k in range(1,n_iter):
        S11 = S11 + np.matmul(xhathat[k],np.transpose(xhathat[k])) + Phat[k]
        S10 = S10 + np.matmul(xhathat[k],np.transpose(xhathat[k-1]))  + Phathat[k]
        S00 = S00 + np.matmul(xhathat[k-1],np.transpose(xhathat[k-1]))  + Phat[k-1]
    
#M-step: update the estimates of u0, sig0, phi, Q, R using 6.70-6.73
    #ignore update of phi?

    phi_new = np.matmul(S10,np.linalg.inv(S00))
    Qnew = 1/n_iter * ( S11 - np.matmul(phi_new,np.transpose(S10)))
    Rsum=np.zeros((2,2))
    for k in range(1,n_iter):
        etaN = np.array([[temperature[k]]])-A*(xhathat[k])
        Rsum = Rsum + np.matmul(etaN ,np.transpose(etaN)) + np.matmul(np.matmul(A,Phat[k]),np.transpose(A))

    Rnew = Rsum / n_iter
    return xhat[:,0,0], P[:,0,0], xhathat[:,0,0], Phathat[:,0,0], thisllike


if  __name__ == "__main__":
    z0=[0,q00,R000]
    _,_,_,_,retlike = KF_compute(years0, temperature, means, z0)
    print(retlike)

    def fwd_optimize(x0):
        _,_,_,_,optlik=KF_compute(years0, temperature, means, x0)
        return -optlik

    from scipy.optimize import minimize
    res = minimize(fwd_optimize, z0, method='nelder-mead',
                   options={'xatol': 1e-3, 'disp': True})
    print(res.x)
    xhat,P,xhathat,_,optlik=KF_compute(years0, temperature, means, res.x)
    plt.figure()
    plt.plot(years0,xhat,color='black')
    plt.plot(years0,xhathat,color='orange',zorder=4)
    plt.fill_between(years0, xhat-P, xhat+P, alpha=0.6)
    plt.plot(years0,means,'go')
    plt.show()
