## EM algorithm using Kalman Filtering, example 6.7 of Shumway and Stoffer
## Page 342-344

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import config

#regenerate the land surface and ocean surface estimates:
#T = 0.292*LandRec +0.708*OceanRec
# LandRec = (T - 0.708*OceanRec)/ 0.292

data2 = pd.read_csv(config.CODEBASE_PATH+"/Common_Data/HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv")
years= data2["Time"].to_numpy()
#years[0]=1850
nyrs = len(years)
temperature = data2["Anomaly (deg C)"].to_numpy()
n_iter=len(temperature) -1 #175
temperature=temperature[0:n_iter]
data_ocn = pd.read_csv(config.CODEBASE_PATH+"/Common_Data/HadSST.4.0.1.0_annual_GLOBE.csv", sep=',')
OceanRec = data_ocn['anomaly'].values
OceanRec_sc =  OceanRec *0.9/0.7
LandRec = (temperature- 0.708*OceanRec_sc[0:n_iter])/ 0.292

moving_aves=np.empty(len(temperature))
moving_aves[:] = np.nan
N=4;cN=int(np.ceil(N/2));fN=int(np.floor(N/2));
for i in range((fN),(len(temperature)-cN+1)):
    lasta=i+cN;firsta=i-fN
    moving_aves[i] = np.mean(temperature[firsta:lasta]);
diff_lin = (moving_aves[:-2]+moving_aves[2:])/2 - moving_aves[1:-1]
diff_mean = (temperature - moving_aves)/np.sqrt(N)
Q=np.cov(diff_mean[(fN+1):-(fN+1)],diff_lin[fN:-(fN)])



##1. Initialize the procedure by selecting starting values for the parameters
u0= np.array([[-0.35],[0]])
sig0 = 0.01
phi = np.array([[1,1],[0,1]]) ## we are keeping this as a random walk, so NOT updating this
#Q =  np.array([[0.005,0.001],[0.001,0.0001]])
R = np.array([[0.36,0.11],[0.11,0.16]])              #[0.25,0.006],[0.006,0.19]])
##xt is 1d, yt is 2d
A = np.array([[1,0],[1,0]]) #emissions matrix


# allocate space for arrays

sz = (n_iter,2,1) # size of array
sz2d=(n_iter,2,2) # state covariance is also 1d



lYthetas = [-1000]
dellYtheta = 100

#iterature until we achieve convergence:

while dellYtheta > 1:
    xhat=np.zeros(sz)      # a posteri estimate of x
    P=np.zeros(sz2d)         # a posteri error estimate


    eta = np.zeros((n_iter,2,1))      # a posteri estimate of x

    xhatminus=np.zeros(sz) # a priori estimate of x
    Pminus=np.zeros(sz2d)    # a priori error estimate
    K=np.zeros(sz2d)         # gain or blending factor
    S=np.zeros((n_iter,2,2))

    xhathat=np.zeros(sz)   # smoothed a priori estimate of x
    Phat=np.zeros(sz2d)      # smoothed posteri error estimate 
    Khat=np.zeros(sz2d)      # smoothed gain or blending factor

    Phathat=np.zeros(sz2d)      # smoothed posteri error estimate 

#E-step: 6.1 (Kalman filter), 6.2 (Kalman Smoother), 6.3 (third covar est) for t=1,..n
    xhat[0] = u0
    P[0,:,:] = sig0
    Pminus[0,:,:] = P[0,:,:]
    thislYtheta=0
#forward loop:
    for k in range(1,n_iter):
        # time update
        
        xhatminus[k] = np.matmul(phi, xhat[k-1])
        
        Pminus[k] = np.matmul(np.matmul(phi,P[k-1]),np.transpose(phi))+Q   #np.matmul(np.matmul(F[k] ,P[k-1]), np.transpose(F[k])) +Q 

        # measurement update if(Rvary):

        S[k]= np.matmul(np.matmul(A ,Pminus[k]), np.transpose(A)) + R

        K[k] = np.matmul(Pminus[k],np.matmul(np.transpose(A),np.linalg.inv(S[k])))
        eta[k]=np.array([[LandRec[k]],[OceanRec_sc[k]]])-np.matmul(A,(xhatminus[k]))

        xhat[k] = xhatminus[k]+np.matmul(K[k],eta[k])
        P[k] = np.matmul((np.eye(2)- np.matmul(K[k],A) ),Pminus[k])
        
        thislYtheta = thislYtheta + 0.5 * np.log(np.linalg.det(S[k])) + np.matmul( np.matmul( np.transpose(eta[k]), np.linalg.inv(S[k])), eta[k])

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
        etaN = np.array([[LandRec[k]],[OceanRec_sc[k]]])-np.matmul(A,(xhathat[k]))
        #np.array([[EBM.temps[k]+EBM.DO_mean_temp],[EBM.IV_meas_fwd[k]]])-np.matmul(A,xhathat[k])
        Rsum = Rsum + np.matmul(etaN ,np.transpose(etaN)) + np.matmul(np.matmul(A,Phat[k]),np.transpose(A))

    Rnew = Rsum / n_iter
    #R = Rsum[::-1,::-1] / n_iter
    #-173.0596179622334 -2.691956143107717
#[[ 1.01421857  0.41972761]
# [-0.00126575  1.10357265]]
#[[-1.37020693e-01  1.41227742e-01]
# [-1.39551416e-01  5.66752453e-05]]
#[[0.36625514 0.08973213]
# [0.08973213 0.18262669]]
    
    #phi=phi_new
    if __name__ == "__main__":
        print(thislYtheta[0,0])
        print(thislYtheta[0,0], dellYtheta)
        print(phi_new)
        print(Q_new)
        print(R_new)

if __name__ == "__main__":
    plt.figure()
    plt.plot(years,LandRec,color='green')
    plt.plot(years,OceanRec_sc[0:n_iter],color='navy')
    plt.plot(years,temperature,color='grey')
    #plt.show()
    plt.plot(years,xhathat[:,0,0],color='black')
    plt.plot(years,xhat[:,0,0],color='magenta')
    plt.show()
