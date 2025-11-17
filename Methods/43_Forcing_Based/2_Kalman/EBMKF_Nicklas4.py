import numpy as np
import sys

import matplotlib.pyplot as plt
import matplotlib        as mpl
import scipy.stats as stats
from sklearn.metrics import r2_score
#plt.rcParams['text.usetex'] = True
#plt.rcParams["font.family"] = "Times New Roman"
#mpl.rcParams['mathtext.fontset'] = "cm"
import matplotlib.transforms as mtransforms
import pdb

#from os.path import expanduser
import config
#import os
#class configclass:
#    def __init__(self):
#        self.CODEBASE_PATH = os.path.expanduser('~/')+"Downloads/Thorne_15_codefigurestats"
#        self.CLIMATE_DATA_PATH = os.path.expanduser('~/')+"climate_data" 
#config = configclass()
common_dataprefix = config.CODEBASE_PATH+"/Common_Data/"

colorekf=(26./255, 44./255, 105./255)
colorstate=(62./255, 207./255, 117./255)
coloruncert=(52./255, 235./255, 235./255)
colorgrey=(0.5,0.5,0.5)

n_iters = 174 

zJ_from_W = 5.1006447*3.154 #*0.71 #CHANGE NOV25
zJtomm=0.121
pcolor='deeppink'
colorrts='goldenrod'
nbins=24
sdate=1850
data = np.genfromtxt(open(common_dataprefix+"toyKFmodelData8c.csv", "rb"),dtype=float, delimiter=',') #run as main
#data = np.genfromtxt(open("../Common_Data/toyKFmodelData8c.csv", "rb"),dtype=float, delimiter=',') #run from other file
dates=data[:,0]
dates[0]=sdate


#temps=data[:, 1]+287
lCo2=np.log10(data[:,2]) 
opt_depth=data[:,3]*0.001 #*0.053) #/1000?0.000100025
anthro_clouds=(data[:,7]+1) #ONLY FOR HISTORICAL RUN

sys.path.append(config.CODEBASE_PATH+"/OHC_ensemble/")
from sample_int_variance import get_average_ohca
comb_ohca, comb_ohca_uncer = get_average_ohca()
comb_ohca_uncer2=comb_ohca_uncer[0:n_iters]
Roc_tvar=np.square(comb_ohca_uncer2)/zJ_from_W/zJ_from_W
ocean_heat_measured =  comb_ohca[0:n_iters]

tsi=data[:,8]
#critexp = np.ones(ln(transmd))

#data2 = np.genfromtxt(open("HadCRUT5.global.annual.csv", "rb"),dtype=float, delimiter=',')
temps=data[:, 1]
R_tvar=np.square(data[:,4]) #still in temperature units



JonesOffset=13.85
offset= -np.mean(temps[(1960-sdate):(1990-sdate)]) +JonesOffset + 273.15 #Jones2013 13.7 to 14 (Jones1999)
#print(offset)
temps=temps+offset
#pindavg= np.mean(temps[(1850-sdate):(1930-sdate)])

heatCp=17 #change to 23 (0.037140) makes it worse NOV25

sig=5.6704e-8

a_refl= 0.834 #constant clearsky albedo
g_refl=0.909
sw_in=340.2




#print("Teq1850",Teq1850)

shaldepth=86
deepdepth = 1141

fdbkS = 0.35 
fdbkA = 0.42 # - 0.12
fdbkA_sigma = 0.316 #dfaA=fdbkA/(sw_in*a_refl*g_refl)
fdbkW = 1.3

def precompute_coeffs(printThem):
    global dfaS, dfaA, powp1, inbndf, rad1850, B1, B0, B1B0, outbndf, Cs, Cd, T02, Teq1850
    T02=np.mean(temps[(2000-1850):(2005-1850)])  #~287.55 #in 2002 #(JMN rev 2024)
    Teq1850=np.mean(temps[0:50]) #CHANGE NOV25 to 286.75 np.mean(temps[1:6])
    #pindavg=Teq1850 #286.7 #preindustrial avg
    dfaS=fdbkS/(sw_in*a_refl*g_refl)
    dfaA=fdbkA/(sw_in*a_refl*g_refl)
    powp1 = fdbkW*4/3.22 #1.3
    B1B0 = 12.74/sig/np.power( T02, 4-powp1 ) #12.74 15.45 #19.45
    inbndf= a_refl*g_refl*9.068 #sw_in* 137.49
    rad1850 = (sw_in*0.9318*a_refl*(1+dfaA*(Teq1850-T02))+anthro_clouds[0])*g_refl*(1+dfaS*(Teq1850-T02))
    B1 = (rad1850 / 5.670e-8 / np.power( Teq1850, 4-powp1 ))+B1B0 *lCo2[0] #594
    B0 = B1B0/B1
    outbndf= sig*B1
    Cs= 136.5*shaldepth/1000 #<17 like maybe 14
    Cd= 136.5*deepdepth/1000
    if(printThem):
        print(T02)
        print(np.mean(temps[-31:-1])-Teq1850)
        print("inbndf = " + str(inbndf*tsi[0]))
        print("B1B0 = " + str(B1B0))
        print("rad1850 = " + str(rad1850))
        print("B1 = " + str(B1))
        print("B0 = " + str(B0))
        print("outbndf = " + str(outbndf))

        
precompute_coeffs(False)


global gad
gad = 0.67
#epd = 1   #1.3
oc1850= 3.5 + 273.15 #absolute temperature of deep ocean
oc_meas = ocean_heat_measured/zJ_from_W #/deepdepth/1.55 + oc1850 #convert into a deep ocean temp


# intial parameters

H=np.eye(2) #emissions matrix

#covariance of the process noise

N = 30
moving_aves=np.empty(len(temps))
moving_aves[:] = np.nan
ocean_aves=np.empty(len(temps))
ocean_aves[:] = np.nan
std_aves=np.empty(len(temps))
std_aves[:] = np.nan
cN=int(np.ceil(N/2));
fN=int(np.floor(N/2));
for i in range((fN),(len(temps)-cN+1)):
    lasta=i+cN;firsta=i-fN
    moving_aves[i] = np.mean(temps[firsta:lasta]);
    ocean_aves[i] = np.mean(oc_meas[firsta:lasta]);
    std_aves[i]= np.std(temps[firsta:lasta]);

resid_temp = (temps - moving_aves);
resid_ocean = (oc_meas - ocean_aves);
covM= np.cov(np.array([resid_temp[cN:-fN],resid_ocean[cN:-fN]]));
#print(covM);



def compute_slope(x2,ki, optdn=-1, lCo2n=-1,anthro_cloud=-101):
    k=int(ki)
    tsik=sw_in
    if (anthro_cloud<-100):
        anthro_cloud = anthro_clouds[k]
    if (optdn<0):
        optdn = opt_depth[k]
    if (lCo2n<=0):
        lCo2n= lCo2[k]
        tsik=tsi[k]
    [x,oc]= x2
    inboundd=tsik*inbndf/(optdn+9.7279)*((dfaS+dfaA)+2*(dfaS*dfaA)*(x-T02)+dfaS*anthro_cloud/sw_in/0.9318/a_refl)
    outgoingd= outbndf*(4-powp1)*np.power(x, 3-powp1)*(1-B0*lCo2n) #0.0655
#following jacobian matrix conventions
    return np.array([[(1 + (inboundd-outgoingd - gad*(1-Cs/Cd))/heatCp ) , gad/Cd/heatCp], \
                     [ (inboundd-outgoingd-gad*(1-Cs/Cd))*Cs/heatCp +gad*(1-Cs/Cd)+Cs  ,  ( 1 - (1-Cs/heatCp)*gad/Cd)]])

def compute_slope_hidden(x2,ki, optdn=-1, lCo2n=-1,anthro_cloud=-101):
    global gad
    k=int(ki)
    tsik=sw_in
    if (anthro_cloud<-100):
        anthro_cloud = anthro_clouds[k]
    if (optdn<0):
        optdn = opt_depth[k]
    if (lCo2n<=0):
        lCo2n= lCo2[k]
        tsik=tsi[k]
    [x,H,eng_bal,inbndf_frac]= x2.T[0] #unpack completely
    oc=(H - Cs*(x-Teq1850))/Cd+oc1850
    inboundd=tsik*inbndf*inbndf_frac/(optdn+9.7279)*((dfaS+dfaA)+2*(dfaS*dfaA)*(x-T02)+dfaS*anthro_cloud/sw_in/0.9318/a_refl)
    outgoingd= outbndf*(4-powp1)*np.power(x, 3-powp1)*(1-B0*lCo2n) #0.0655

    inbound_dinbndf_frac=tsik*inbndf*(1+dfaA*(x-T02)+anthro_cloud/sw_in/0.9318/a_refl)*(1+dfaS*(x-T02))/(optdn+9.7279)
    
#following jacobian matrix conventions

    return np.array([[(1 + (inboundd-outgoingd - gad*(1-Cs/Cd))/heatCp ) , gad/Cd/heatCp, 0, (inbound_dinbndf_frac)/heatCp], \
                     [ (inboundd-outgoingd-gad*(1-Cs/Cd))*Cs/heatCp +gad*(1-Cs/Cd)+Cs  ,  ( 1 - (1-Cs/heatCp)*gad/Cd),0 , (inbound_dinbndf_frac)*Cs/heatCp],\
                     [inboundd-outgoingd,0,0,inbound_dinbndf_frac],\
                     [0,0,0,1]]) #inputf_float just randomly walks
                     #[0,0,0,0,1]]) #gad_float just randomly walks -(x-Teq1850 - oc + oc1850) /heatCp, ,(x-Teq1850 - oc + oc1850)

def compute_update(x2,ki, optdn=-1, lCo2n=-1,anthro_cloud=-101): 
    k=int(ki)
    [x,H]= x2
    tsik=sw_in
    if (anthro_cloud<-100):
        anthro_cloud = anthro_clouds[k]
    if (optdn<0):
        optdn = opt_depth[k]
    if (lCo2n<=0):
        lCo2n= lCo2[k]
        tsik=tsi[k]
    #print(anthro_cloud);
    inbound=tsik*inbndf*(1+dfaA*(x-T02)+anthro_cloud/sw_in/0.9318/a_refl)*(1+dfaS*(x-T02))/(optdn+9.7279)
    outgoing= outbndf*np.power(x, 4-powp1 )*(1-B0*lCo2n)
    oc=(H - Cs*(x-Teq1850))/Cd+oc1850
    if (x)<280:
        raise Exception("Too Cold "+str(k+1850))
    #by convention addition of prior x, oc is not included
    Tchange=(inbound - outgoing - gad*(x-Teq1850 - oc + oc1850) )/heatCp
    return np.array([Tchange, Tchange*Cs + gad*(x-Teq1850 - oc + oc1850) , inbound-outgoing])

def compute_update_hidden(x2,ki, optdn=-1, lCo2n=-1,anthro_cloud=-101):
    global gad
    k=int(ki)
    #breakpoint()
    [x,H,eng_bal,inbndf_frac]= x2.T[0]    #we don't care about energy balance at last time point, so 3rd column is 0s
    #print(x)
    tsik=sw_in
    if (anthro_cloud<-100):
        anthro_cloud = anthro_clouds[k]
    if (optdn<0):
        optdn = opt_depth[k]
    if (lCo2n<=0):
        lCo2n= lCo2[k]
        tsik=tsi[k]
    #print(anthro_cloud);
    inbound=tsik*inbndf*inbndf_frac*(1+dfaA*(x-T02)+anthro_cloud/sw_in/0.9318/a_refl)*(1+dfaS*(x-T02))/(optdn+9.7279)
    outgoing= outbndf*np.power(x, 4-powp1 )*(1-B0*lCo2n)
    oc=(H - Cs*(x-Teq1850))/Cd+oc1850
    if (x)<280:
        raise Exception("Too Cold "+str(k+1850))
    #by convention addition of prior x, oc is not included

    Tchange=(inbound - outgoing - gad*(x-Teq1850 - oc + oc1850) )/heatCp
    return np.array([[x + Tchange, H + Tchange*Cs + gad*(x-Teq1850 - oc + oc1850) , inbound-outgoing,inbndf_frac]]).T


# allocate space for arrays
        
A = np.array([[1,0,0,0],
              [0,1,0,0],
              [0,0,1,0]])
Q=(covM/N);


sz = (n_iters,4,1) # size of array
sz2d=(n_iters,4,4)

TOA_swout_high = np.zeros(n_iters)
TOA_swout_low = np.zeros(n_iters)
TOA_swout = np.zeros(n_iters)
#TOA_swout[:] = np.nan
xhat=np.zeros(sz)      # a posteri estimate of x
P=np.zeros(sz2d)         # a posteri error estimate
F=np.zeros(sz2d)         # state transitions
xhatminus=np.zeros(sz) # a priori estimate of x
Pminus=np.zeros(sz2d)    # a priori error estimate
K=np.zeros((n_iters,4,3))         # gain or blending factor

xhathat=np.zeros(sz)   # smoothed a priori estimate of x
Phat=np.zeros(sz2d)      # smoothed posteri error estimate 
Khat=np.zeros((n_iters,4,4))      # smoothed gain or blending factor
Shat=np.zeros((n_iters,3,3))
S=np.zeros((n_iters,3,3))
xblind=np.zeros((n_iters,2,1))


lml=np.zeros(sz)
lsml=0
y=np.zeros((n_iters,3,1))
qqy=np.zeros((n_iters,3))

qqyh=[]
qqyk=[]

xnorm = np.linspace(-5.5, 5.5, 200)
qqyker=np.zeros((2,200))
qqykerall=[]
ynorm  = stats.norm.pdf(xnorm)
involcavg = np.mean(1/(opt_depth+9.7279))

import xarray as xr
import os
def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])

plt.figure()
plt.plot(dates, TOA_swout)
toa_file = config.CLIMATE_DATA_PATH+"/NASA/global_rt_mon_CERES-EBAF-4-2_RSS_gn_200003-202310.nc"
with xr.open_dataset(toa_file) as ds:
    TOA_measur = ds['rt'].values
    TOA_time = ds['time'].values

TOA_crop = average_every_n(TOA_measur[9: int(np.floor((len(TOA_measur)-9)/12))*12 + 9],12) #starts in month 3 of 2000
from scipy import stats
regX = np.arange(2001,2023)
regres = stats.linregress(regX, TOA_crop) #can also give stdrr of slope, intercept
y_pred = regres.slope * regX + regres.intercept
residuals_TOA = TOA_crop - y_pred
TOA_crop_var = np.var(residuals_TOA)

regres = stats.linregress(np.arange(1975,2023), temps[1975 - 1850: 2023 - 1850]) #repeat for temp with linear reg
y_pred = regres.slope * regX + regres.intercept
residuals_templ = temps[2001 - 1850: 2023 - 1850] - y_pred

regres = stats.linregress(np.arange(1975,2023), oc_meas[1975 - 1850: 2023 - 1850]) #repeat for ocean with linear reg
y_pred = regres.slope * regX + regres.intercept
residuals_oceanl = oc_meas[2001 - 1850: 2023 - 1850] - y_pred


cM2= np.cov(np.array([residuals_templ,residuals_oceanl, residuals_TOA]))/N *4;
dfrA_float_var = (fdbkA_sigma/(sw_in*a_refl*g_refl))**2


xblind[0,0:2,0]= np.array([Teq1850,oc_meas[0]])
for k in range(1,n_iters):
    temp_update = compute_update_hidden(
        np.array([[xblind[k-1,0,0],xblind[k-1,1,0],0,1]]).T,k) #first compute the expected values of TOA radiation
    xblind[k]=temp_update[0:2]
    TOA_swout[k] = temp_update[2,0]
    TOA_swout_high[k] = compute_update_hidden(
        np.array([[xblind[k-1,0,0],xblind[k-1,1,0],0,1 + fdbkA_sigma/(sw_in*a_refl*g_refl)]]).T,k)[2,0]
    TOA_swout_low[k] = compute_update_hidden(
        np.array([[xblind[k-1,0,0],xblind[k-1,1,0],0,1 - fdbkA_sigma/(sw_in*a_refl*g_refl)]]).T,k)[2,0]

#recombining
TOA_meas_artif = np.copy(TOA_swout)
TOA_swout2 = np.copy(TOA_swout)
TOA_meas_artif_var = (TOA_swout_low - TOA_swout_high)**2/16
TOA_meas_artif_var[125:]=TOA_meas_artif_var[125]
TOA_meas_artif[2001-1850:2023-1850] = TOA_crop
TOA_meas_artif_var[2001-1850:2023-1850] = TOA_crop_var/9
#plt.plot(, TOA_crop)
#plt.title("TOA radiative imbalance")



def ekf_run(z,n_iter,retPs=False):
    global gad
    xblind2=np.zeros((n_iter,3,1))
    xblind2[0,0:3,0]= np.array([Teq1850,oc_meas[0],0])
    for k in range(1,n_iter):
        temp_update = compute_update_hidden(
           np.array([[xblind2[k-1,0,0],xblind2[k-1,1,0],0,1]]).T,k) ##xblind not including correction to forcing
        xblind2[k]=temp_update[0:3]

    if(n_iter>(2030-1850)):
        gad = gad * ( .5 + .5 * z[(2020-1850),1,0]/ xblind2[(2020-1850),1,0]) #scale gad to reflect z

    # intial guesses
    xhat[0] = np.array([[Teq1850,oc_meas[0],0,1]]).T #start with inbndf_frac=1
    P[0,:,:] = 1.0
    P[0,1,1] = 20
    Pminus[0,:,:] = P[0,:,:]


    #F[0]=compute_slope(xhat[0,:,0],0) #necessary for last step of RTS smoother

    
    
    for k in range(1,n_iter):
        # time update

        
            
        F[k]=compute_slope_hidden(xhat[k-1],k)   #dx vary along columns

        xhatminus[k] = compute_update_hidden(xhat[k-1],k)
        
        #TOA_swout[k]=  xhatminus[k][2,0]
        sc_old=1
        Pminus[k] = np.matmul(np.matmul(F[k] ,P[k-1]), np.transpose(F[k])) +  np.array([[Q[0,0]   ,Q[0,1]   ,cM2[0,2] , 0],
                                                                                        [Q[1,0]  ,Q[1,1]    ,cM2[1,2] , 0],
                                                                                        [cM2[2,0]  ,cM2[2,1] ,cM2[2,2] , 0],
                                                                                        [0       ,0       ,0,    dfrA_float_var/N]]) #,dfrA_float_var/200 ###30
                                                                                        #[0       ,0       ,0       ,0        , 0.14**2/40]])  # stuff  dfrA_float_var


        # change 30 -> 10 to reduce additional measurement noise
        # change dfrA_float_var/60 due to overdispersion of CMIP

        S[k]= np.matmul(np.matmul(A ,Pminus[k]), np.transpose(A))+ N* np.array([[Q[0,0]  ,Q[0,1]  ,cM2[0,2]],
                                                                                 [Q[1,0]  ,Q[1,1]  ,cM2[1,2]],
                                                                                 [cM2[2,0],cM2[2,1],cM2[2,2]]]) \
                                + np.matrix([[R_tvar[k],0      ,  0],
                                             [0   ,   Roc_tvar[k],0],
                                             [0   ,     0 ,       9*TOA_meas_artif_var[k]]]) #3x3 stuff, time-varying


            
        K[k] = np.matmul(Pminus[k],np.matmul(np.transpose(A),np.linalg.inv(S[k])))
        y[k]=z[k]-np.matmul(A,(xhatminus[k]))
        xhat[k] = xhatminus[k]+np.matmul(K[k],y[k])
        P[k] = np.matmul((np.eye(4)- np.matmul(K[k],A) ),Pminus[k])
        stdevS=np.sqrt(np.abs(np.diag(S[k])))
        qqy[k]=(y[k].T/stdevS[0:3])
        gmstpdf=stats.norm.pdf(xnorm,loc = qqy[k,0], scale = np.sqrt(R_tvar[k])/stdevS[0])
        ohcapdf=stats.norm.pdf(xnorm,loc = qqy[k,1], scale = np.sqrt(Roc_tvar[k])/stdevS[1])
        qqykerall.append([gmstpdf,ohcapdf])
#        if k==120: STOP ON A PARTICULAR YEAR FOR DEBUGGING
#            breakpoint() #second iteration
        if k>1:
            qqyker[0] += gmstpdf/(n_iters-1)
            qqyker[1] += ohcapdf/(n_iters-1)
        


    xhathat[n_iter-1]=xhat[n_iter-1]
    Phat[n_iter-1]=P[n_iter-1]
    xhathat[0]=xhat[0]
    Phat[0]=P[0]

    lsml=0
    ybark=0
    

    if(True): #do RTS now
        for ik in range(2,n_iter+1):
        # RTS Smoother
            k=n_iter-ik
        # measurement update
            try:
                Khat[k] = np.matmul(np.matmul(P[k],np.transpose(F[k+1])),np.linalg.inv(Pminus[k+1])) #compute inverse for higher dimensions
            except:
                Khat[k] =Khat[k+1]
            #breakpoint()
            xhathat[k] = xhat[k]+np.matmul(Khat[k],(xhathat[k+1]- compute_update_hidden(xhat[k],k)))
            Phat[k] = P[k] + np.matmul(np.matmul(Khat[k],(Phat[k+1]- Pminus[k])),np.transpose(Khat[k]))
            yrts=z[k]-np.matmul(A,(xhathat[k])) #z[k]-H *xhathat[k]
            Shat[k]= np.matmul(np.matmul(A ,Phat[k]), np.transpose(A)) + 30* np.array([[Q[0,0]  ,Q[0,1]  ,cM2[0,2]],
                                                                                 [Q[1,0]  ,Q[1,1]  ,cM2[1,2]],
                                                                                 [cM2[2,0],cM2[2,1],cM2[2,2]]]) \
                                + np.matrix([[R_tvar[k],0      ,  0],
                                             [0   ,   Roc_tvar[k],0],
                                             [0   ,     0 ,       TOA_meas_artif_var[k]]]) #3x3 stuff, time-varying

    if (retPs==True):
        return xhat[0:n_iter], P[0:n_iter]
    elif (retPs==2):
        return xhat[0:n_iter], P[0:n_iter], S[0:n_iter], xhatminus[0:n_iter,0]
    elif (retPs==3):
        #breakpoint()
        if True:
            
            dates = np.arange(1850+30,1850+n_iter)
            datesAll = np.arange(1850,1850+n_iter)
            fig = plt.figure()
##            xh1gad=xhat[:,4].T[0] 
##            plt.plot(dates,xh1gad[30:],color='purple')
##            stdPgad=np.sqrt(np.abs(P))[:,4,4]
##            plt.fill_between(dates,xh1gad[30:] + 2*stdPgad[30:], xh1gad[30:] - 2*stdPgad[30:], color = 'blue')
            xh1fdbk=(xhat[:,3].T[0]-1) *  (sw_in*a_refl*g_refl)
            plt.plot(datesAll,xh1fdbk)
            stdPfdbk=np.sqrt(np.abs(P))[:,3,3] *(sw_in*a_refl*g_refl)
            plt.fill_between(dates,xh1fdbk[30:] + 2*stdPfdbk[30:], xh1fdbk[30:] - 2*stdPfdbk[30:], color = 'red')
            #plt.axhline(y=fdbkA,color='k')

            xblind3=np.zeros((n_iter,3,1))
            xblind3[0,0:3,0]= np.array([Teq1850,oc_meas[0],0])
            for k in range(1,n_iter):
                temp_update = compute_update_hidden(
                np.array([[xblind3[k-1,0,0],xblind3[k-1,1,0],0,xh1fdbk[k]/(sw_in*a_refl*g_refl)+1]]).T,k) #xblind while yes including correction to forcing
                xblind3[k]=temp_update[0:3]
            
            plt.figure()
            xh1TOA=xhat[:,2].T[0]
            stdPtoa=np.sqrt(np.abs(P))[30:,2,2]
            stdStoa=np.sqrt(np.abs(S))[30:,2,2]
            plt.plot(datesAll,xh1TOA,'-',label='posterior OHCA EBM-KF-uf state estimate $\^{H }_{t}$', color=colorekf)
            plt.fill_between(dates, (xh1TOA[30:]-2*stdStoa), (xh1TOA[30:]+2*stdStoa),label="95% CI ($\pm 2\sqrt{\hat{s}^H_t}$) of forecast OHCA $\^{H }_{t|t-1}$", color=coloruncert)
            plt.fill_between(dates, (xh1TOA[30:]-2*stdPtoa), (xh1TOA[30:]+2*stdPtoa),label="95% CI ($\pm 2\sqrt{\hat{p}^H_t}$) of OHCA state $\^{H }_{t}$", color=colorstate)
            TOA_meas_artif=z[:,2].T[0]
            plt.plot(datesAll, TOA_meas_artif,'o',markersize=2,color=colorgrey)
            toa_xblind0=xblind2[:,2].T[0]
            plt.plot(datesAll, toa_xblind0,'o',markersize=2,color='red')
            toa_xblind=xblind3[:,2].T[0]
            plt.plot(datesAll, toa_xblind,color='orange')
            plt.fill_between(datesAll, TOA_meas_artif-2*np.sqrt(TOA_meas_artif_var), TOA_meas_artif+2*np.sqrt(TOA_meas_artif_var), color="grey",zorder=5,alpha=0.2,lw=0)

            plt.figure()
            xh1ohca=xhat[:,1].T[0]
            stdPohca=np.sqrt(np.abs(P))[30:,1,1]
            stdSohca=np.sqrt(np.abs(S))[30:,1,1]
            plt.plot(datesAll,xh1ohca,'-',label='posterior OHCA EBM-KF-uf state estimate $\^{H }_{t}$', color=colorekf)
            plt.fill_between(dates, (xh1ohca[30:]-2*stdSohca), (xh1ohca[30:]+2*stdSohca),label="95% CI ($\pm 2\sqrt{\hat{s}^H_t}$) of forecast OHCA $\^{H }_{t|t-1}$", color=coloruncert)
            plt.fill_between(dates, (xh1ohca[30:]-2*stdPohca), (xh1ohca[30:]+2*stdPohca),label="95% CI ($\pm 2\sqrt{\hat{p}^H_t}$) of OHCA state $\^{H }_{t}$", color=colorstate)
            ohca_meas=z[:,1].T[0]
            plt.plot(datesAll, ohca_meas,'o',markersize=2,color=colorgrey)
            ohca_xblind=xblind2[:,1].T[0]
            plt.plot(datesAll, ohca_xblind,'o',markersize=2,color='red')
            ohca_xblind3=xblind3[:,1].T[0]
            plt.plot(datesAll, ohca_xblind3,'o',markersize=2,color='orange')

            plt.figure()
            xh1T=xhat[:,0].T[0]
            stdPT=np.sqrt(np.abs(P))[30:,0,0]
            stdST=np.sqrt(np.abs(S))[30:,0,0]
            plt.plot(datesAll,xh1T,'-',label='posterior OHCA EBM-KF-uf state estimate $\^{H }_{t}$', color=colorekf)
            plt.fill_between(dates, (xh1T[30:]-2*stdST), (xh1T[30:]+2*stdST),label="95% CI ($\pm 2\sqrt{\hat{s}^H_t}$) of forecast OHCA $\^{H }_{t|t-1}$", color=coloruncert)
            plt.fill_between(dates, (xh1T[30:]-2*stdPT), (xh1T[30:]+2*stdPT),label="95% CI ($\pm 2\sqrt{\hat{p}^H_t}$) of OHCA state $\^{H }_{t}$", color=colorstate)
            T_meas=z[:,0].T[0]
            plt.plot(datesAll, T_meas,'o',markersize=2,color=colorgrey)
            T_xblind=xblind2[:,0].T[0]
            plt.plot(datesAll, T_xblind,'o',markersize=2,color='red')
            T_xblind3=xblind3[:,0].T[0]
            plt.plot(datesAll, T_xblind3,'o',markersize=2,color='orange')
            
            #plt.plot(dates, TOA_swout,'o',markersize=2,color='red')
            #plt.fill_between(dates, TOA_meas_artif-2*np.sqrt(TOA_meas_artif_var[30:]), TOA_meas_artif+2*np.sqrt(TOA_meas_artif_var[30:]), color="grey",zorder=5,alpha=0.2,lw=0)
        return xhat[0:n_iter,0,0], P[0:n_iter,0,0],xhathat[0:n_iter,0,0], Phat[0:n_iter,0,0]
    elif (retPs==4):
        return xhat[0:n_iter,0,0], P[0:n_iter,0,0],xblind[0:n_iter,0,0], Phat[0:n_iter,0,0]
    else:
        return xhat[0:n_iter]


observ = np.array([[temps,oc_meas,TOA_meas_artif]]).T
#print(observ)
this_xhat=ekf_run(observ,n_iters)
xh1s=this_xhat[:,0].T[0]
xh0s=xhatminus[:,0].T[0]
stdS=np.sqrt(np.abs(S))[:,0,0]
stdP=np.sqrt(np.abs(P))[:,0,0]
Plastretain=np.abs(P)[-1,:,:]
xhh1s=xhathat[:,0].T[0]
xlastretain=this_xhat[-1,:]
stdPh=np.sqrt(np.abs(Phat))[:,0,0]
stdSh=np.sqrt(np.abs(Shat))[:,0,0]

xh1d=this_xhat[:,1].T[0]
xh0d=xhatminus[:,1].T[0]
stdSd=np.sqrt(np.abs(S))[:,1,1]
stdPd=np.sqrt(np.abs(P))[:,1,1]

#print(1-stats.norm.cdf(1.5+pindavg,xh1s[-1],stdP[-1]))
#print(1-stats.norm.cdf(1.5+pindavg,xh1s[-1],stdS[-1]))
#print(temps[-1]-pindavg)

Rtvara=np.mean(R_tvar[-11:-1])+30*Q[0,0]
Roctvara=np.mean(Roc_tvar[-11:-1])+30*Q[1,1]

def plot_boilerplate(ax=plt.gca()):
    ax.set_ylim(286.1,288.3)
   # plt.yticks(np.arange(286.2,288.4,0.2))
    ax.set_xticks(np.arange(1850,2025+1,25))
    ax.set_xlim(1850,2025)
    ax.set_xlabel('Year')
    ax.set_ylabel('Temperature (K)')
    ax.tick_params( direction = 'in',bottom=True, top=True, left=True, right=True )

if (__name__ == "__main__") and True:
    r2 = r2_score(observ[:,0], xblind[:,0])
    print('r2 score for GMST to blind is', r2)
    r2 = r2_score(observ[:,1], xblind[:,1])
    print('r2 score for OCHA to blind is', r2)

    r2 = r2_score(moving_aves[int(np.floor(N/2)):int(-(np.floor(N/2)+1))], xh1s[int(np.floor(N/2)):int(-(np.floor(N/2)+1))]) #
    print('r2 score for 30-mean GMST to EBM-KF is', r2)
    r2 = r2_score(ocean_aves[int(np.floor(N/2)):int(-(np.floor(N/2)+1))], xh1d[int(np.floor(N/2)):int(-(np.floor(N/2)+1))])
    print('r2 score for 30-mean OCHA to EBM-KF is', r2)

    #plt.rcParams['figure.figsize'] = (7, 6)

   # plt.figure(1)
    fig, (ax1,ax4)= plt.subplots(2, 1, figsize=(7,10), gridspec_kw={ "hspace": 0.3})
    plot_boilerplate(ax1)
    plot_boilerplate(ax4)
    ax1.plot(dates,temps,'o',label='noisy HadCRUT5 GMST measurements $Y_{t}$',markersize=2,color=colorgrey)
    ax1.fill_between(dates, temps-2*data[:,4], temps+2*data[:,4],label="associated 95% uncertainty of $Y_{t}$", color="lightgrey")
   # ax1.plot(dates,this_xhat[:,0],'-',label='posterior GMST EBM-KF-uf state estimate $\^{T }_{t}}$', color=colorekf)
    ax1.plot(dates,xblind[:,0],'--',label='blind model $\~{T}_{t+1} = \mathbf{F}( \~{T}_{t}, \~{H}_{t} ; [eCO_2]_t ,AOD_t,AC_{t},(\\frac{1}{4}G_{SC})_t )$',color='darkgoldenrod')
    #ax1.legend(loc="upper left", fontsize="10.5")
    handles, labels = ax1.get_legend_handles_labels()
    order = [1,0,2]#[1,2,0,3]
    ax1.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc="upper left",prop={'size': 10.5})
    #fig.suptitle('Kalman Filter is a "weighted average" between Model Projections and Noisy Data')
    ax1.set_title('Global Mean Surface Temperature State')
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_ylim(286.1,288.3)
    ax1.set_yticks(np.arange(286.2,288.3,0.2))
    ax1.set_xlim(1850,2025)
    ax1b = ax1.twinx()
    mn, mx = ax1.get_ylim()
    ax1b.set_ylim(mn-273.15, mx-273.15)
    ax1b.set_yticks(np.arange(13,15.2,0.2))
    ax1b.set_ylabel('Temperature (°C)')
    #plt.savefig("Figure_1rev.pdf", format="pdf") #dpi=300,format="png")


    #plt.rcParams['figure.figsize'] = (7, 6)
    #plt.figure(2)
    plot_boilerplate(ax4)
    ax4.plot(dates,ocean_heat_measured,'o',label='noisy Zanna (2019) measurements $\Psi_{t}$',markersize=2,color=colorgrey)
    ax4.fill_between(dates, ocean_heat_measured-2*comb_ohca_uncer2, ocean_heat_measured+2*comb_ohca_uncer2,label="associated 95% uncertainty of $ \Psi_{t}$", color="lightgrey")
    #ax4.plot(dates,this_xhat[:,1]*zJ_from_W,'-',label='posterior OHCA EBM-KF-uf state estimate $\^{H }_{t}}$', color=colorekf)
    ax4.plot(dates,xblind[:,1]*zJ_from_W,'--',label='blind model $ \~{H}_{t+1} = \mathbf{F}( \~{T}_{t}, \~{H}_{t} ; [eCO_2]_t ,AOD_t,AC_{t},(\\frac{1}{4}G_{SC})_t )$',color='darkgoldenrod')
    #ax4.legend(loc="upper left", fontsize="10.5")
    handles, labels = ax4.get_legend_handles_labels()
   # order = [1,2,0,3]
    ax4.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc="upper left",prop={'size': 10.5})
    ax4.set_title('Ocean Total Heat State')
    ax4.set_xlabel('Year')
    ax4.set_ylabel('Heat (ZJ)')
    ax4.set_ylim([-150,600])
    #ax4b = ax4.twinx()
    #mn, mx = ax4.get_ylim()
    #ax4b.set_yticks(np.arange(-1,8,1))
    #ax4b.set_ylim(mn*zJtomm/10, mx*zJtomm/10)
    
    #ax4b.set_ylabel('Thermosteric Sea Level Rise (cm)')

    axes=(ax1,ax4)
    axlabels=['a)','b)','c)','d)']
    for i in range(2):
        ax=axes[i]
        label=axlabels[i]
          # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-(40)/72, 15/72, fig.dpi_scale_trans) #-20/72, 7/72
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, fontsize='large', va='bottom')

    #fig.savefig("weighted_average_pres.pdf",format="pdf")
    #fig.savefig("weighted_average_pres.png", dpi=400,format="png")


    fig = plt.figure(figsize=(8.25,6))
    #(H - Cs*(x-Teq1850))/Cd+oc1850
    plt.plot(dates,(oc_meas - heatCp*(temps-Teq1850))/Cd+oc1850,'o',markersize=2,label='inferred from both Zanna and HADCRUT5 measurements',color=colorgrey)
    plt.plot(dates,(this_xhat[:,1]- heatCp*(this_xhat[:,0]-Teq1850))/Cd+oc1850,'-', label="posterior ocean cons. temp. EBM-KF state estimate $ \^{\\theta }_{t} $",color=colorekf)
                                                                                                                
    plt.plot(dates,(xblind[:,1]- heatCp*(xblind[:,0]-Teq1850))/Cd+oc1850,'--',label='blind model prediction $\~{\\theta }_{t}$',color='darkgoldenrod')
    plt.legend()
    ax=plt.gca()
    ax.set_ylabel('Temperature (K)')
    mn, mx = ax.get_ylim()
    ax1b = ax.twinx()
    ax1b.set_ylim(mn-273.15, mx-273.15)
    ax1b.set_ylabel('Temperature (°C)')
    ax.set_title("Deep Ocean Conservative Temperature $\\theta$")
    fig.savefig("DOPT.pdf",format="pdf")
    fig.savefig("DOPT.png", dpi=400,format="png")   

    fig = plt.figure(figsize=(8.25,6))
    ax_dict = fig.subplot_mosaic(
    """
    ab
    ac
    """,
    gridspec_kw={
        "width_ratios": [7, 1.25], "wspace": 0.4, "hspace": 0.3, "left":0.1, "right":0.96
    },)
    plt.sca(ax_dict["a"])
    plot_boilerplate(ax_dict["a"])
    plt.yticks(np.arange(286.2,288.3,0.2))
    #plt.plot(dates,xhatminus,'b-',label='a priori EBM-KF estimate', linewidth=3.0)
    ax_dict["a"].set_title('Estimated Climate State with Extended Kalman Filter (unfiltered)')
    ax_dict["a"].plot(dates,this_xhat[:,0],'-',label='posterior EBM-KF-uf state estimate $\^{T }_{t}$', color=colorekf)
  #  plt.plot(dates,xh0s, 'k.')
    ax_dict["a"].fill_between(dates, xh0s-2*stdS, xh0s+2*stdS,label="95% CI ($\pm 2\sqrt{\hat{s}^T_t}$) of forecast GMST $ \^{T }_{t|t-1} $", color=coloruncert)
    ##plt.fill_between(dates, xh1s-stdS, xh1s+stdS, color=(62./255, 140./255, 210./255))
    ax_dict["a"].fill_between(dates, xh1s-2*stdP, xh1s+2*stdP,label="95% CI ($\pm 2\sqrt{\hat{p}^T_t}$) of GMST state $\^{T }_{t}$", color=colorstate)
    ax_dict["a"].plot(dates,temps,'o',label='HadCRUT5 GMST measurements $Y_{t}$',markersize=2,color=colorgrey)
    ax_dict["a"].fill_between(dates, temps-2*data[:,4], temps+2*data[:,4],label="associated 95% uncertainty of $Y_{t}$", color="grey",alpha=0.2,zorder=5,lw=0)
    #plt.fill_between(dates, xhh1s-std, xhh1s+std, color=(171./255, 245./255, 206./255))
    #ax_dict["a"].legend(fontsize="10.5")
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [0,3,2,1,4]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],prop={'size': 10.5})
    ##plt.savefig("Figure_2.pdf", format = "pdf") #dpi=300,format="png")
    ax1b = ax_dict["a"].twinx()
    mn, mx = ax_dict["a"].get_ylim()
    ax1b.set_ylim(mn-273.15, mx-273.15)
    ax1b.set_yticks(np.arange(13,15.2,0.2))
    ax1b.set_ylabel('Temperature (°C)')
    #plt.rcParams['figure.figsize'] = (5, 4)
    ax1 = ax_dict["b"]
    ax2 = ax_dict["c"]
    #plt.subplots_adjust(wspace=0.4)
    stats.probplot(qqy[:,0], dist="norm", plot=ax2)
    ax1.set_title("Innovations")
    ax2.get_lines()[0].set_markerfacecolor(colorekf)
    ax2.get_lines()[0].set_markeredgewidth(0)
    ax2.get_lines()[1].set_color(pcolor)
    ax2.set_xlabel("Theoretical Quantiles",labelpad=0)
    ax2.set_title("")
    #ax1.hist(qqy[:,0], density=True,bins=nbins,color=colorekf)
    ax1.plot(xnorm,qqyker[0,:], color=colorekf)
    ax1.plot(xnorm,ynorm, color=pcolor)
    ax1.set_xlabel("Forecast Std Dev $\sqrt{\hat{s}^T_t}$",labelpad=-3)
    ax1.set_ylabel("Probability Density")
    ax1.set_xlim([-3,3])
  #  plt.rcParams['figure.figsize'] = (10, 8)
    axlabels=['a)','b)','c)']
    axs=(ax_dict["a"],ax1,ax2)
    for i in range(len(axs)):
        ax=axs[i]
        label=axlabels[i]
    # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-40/72, 1/72, fig.dpi_scale_trans) #-20/72, 7/72
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='large', va='bottom') #, fontfamily='serif')
    plt.savefig("mainGMST.pdf", format = "pdf")
    plt.savefig("mainGMST.png", dpi=400,format="png")




    fig = plt.figure(figsize=(8.25,6))
    ax_dict = fig.subplot_mosaic(
    """
    ab
    ac
    """,
    gridspec_kw={
        "width_ratios": [7, 1.25], "wspace": 0.4, "hspace": 0.3, "left":0.1, "right":0.96
    },)
    plt.sca(ax_dict["a"])
    plot_boilerplate(ax_dict["a"])
    
    #plt.plot(dates,xhatminus,'b-',label='a priori EBM-KF estimate', linewidth=3.0)
    ax_dict["a"].set_title('Estimated Ocean Total Heat State with EBM-KF-uf')
    ax_dict["a"].set_ylabel('Heat (ZJ)')
    ax_dict["a"].plot(dates,xh0d*zJ_from_W,'-',label='posterior OHCA EBM-KF-uf state estimate $\^{H }_{t}$', color=colorekf)
  #  plt.plot(dates,xh0s, 'k.')
    ax_dict["a"].fill_between(dates, (xh0d-2*stdSd)*zJ_from_W, (xh0d+2*stdSd)*zJ_from_W,label="95% CI ($\pm 2\sqrt{\hat{s}^H_t}$) of forecast OHCA $\^{H }_{t|t-1}$", color=coloruncert)
    ##plt.fill_between(dates, xh1s-stdS, xh1s+stdS, color=(62./255, 140./255, 210./255))
    ax_dict["a"].fill_between(dates, (xh1d-2*stdPd)*zJ_from_W, (xh1d+2*stdPd)*zJ_from_W,label="95% CI ($\pm 2\sqrt{\hat{p}^H_t}$) of OHCA state $\^{H }_{t}$", color=colorstate)
    ax_dict["a"].plot(dates,ocean_heat_measured,'o',label='Zanna (2019) measurements $\Psi_{t}$',markersize=2,color=colorgrey)
    ax_dict["a"].fill_between(dates, ocean_heat_measured-2*comb_ohca_uncer2, ocean_heat_measured+2*comb_ohca_uncer2,label="associated 95% uncertainty of $\Psi_{t}$", color="grey",zorder=5,alpha=0.2,lw=0)
    #plt.fill_between(dates, xhh1s-std, xhh1s+std, color=(171./255, 245./255, 206./255))
    #ax_dict["a"].legend(fontsize="10.5")
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [0,3,2,1,4]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],prop={'size': 10.5})
    ##plt.savefig("Figure_2.pdf", format = "pdf") #dpi=300,format="png")
    ax_dict["a"].set_ylim([-150,600])
    #ax4b = ax_dict["a"].twinx()
    #mn, mx = ax_dict["a"].get_ylim()
    #ax4b.set_ylim(mn*zJtomm/10, mx*zJtomm/10)
    #ax4b.set_yticks(np.arange(-1,8,1))
    #ax4b.set_ylabel('Thermosteric Sea Level Rise (cm)')

    #plt.rcParams['figure.figsize'] = (5, 4)
    ax1 = ax_dict["b"]
    ax2 = ax_dict["c"]
    #plt.subplots_adjust(wspace=0.4)
    stats.probplot(qqy[:,1], dist="norm", plot=ax2)
    ax1.set_title("Innovations")
    ax2.get_lines()[0].set_markerfacecolor(colorekf)
    ax2.get_lines()[0].set_markeredgewidth(0)
    ax2.get_lines()[1].set_color(pcolor)
    ax2.set_xlabel("Theoretical Quantiles",labelpad=0)
    ax2.set_title("")
    #ax1.hist(qqy[:,1], density=True,bins=nbins,color=colorekf)
    ax1.plot(xnorm,qqyker[1,:], color=colorekf)
    ax1.plot(xnorm,ynorm, color=pcolor)
    ax1.set_xlabel("Forecast Std Dev $\sqrt{\hat{s}^H_t}$",labelpad=-3)
    ax1.set_ylabel("Probability Density")
    ax1.set_xlim([-3,3])
  #  plt.rcParams['figure.figsize'] = (10, 8)
    axlabels=['a)','b)','c)']
    axs=(ax_dict["a"],ax1,ax2)
    for i in range(len(axs)):
        ax=axs[i]
        label=axlabels[i]
    # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-40/72, 1/72, fig.dpi_scale_trans) #-20/72, 7/72
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='large', va='bottom') #, fontfamily='serif')



    plt.rcParams['figure.figsize'] = (7, 6)
    plt.figure()
    
    plt.title('Slight Changes with RT Smoother')
    plt.plot(dates,xhathat[:,0],'-',label='RTS smoothed GMST estimate $ \^{\^{T }}_{t}$',color=colorrts)
    plt.plot(dates,xhat[:,0],'-',label='posterior GMST EBM-KF state estimate $\^{T }_{t}$',color=colorekf)


    plt.fill_between(dates[4:], xhh1s[4:]-2*stdPh[4:], xhh1s[4:]+2*stdPh[4:],label="95% CI $(\pm 2\sqrt{\hat{\hat {p}}^T_t})$ of RTS GMST state $\^{\^{T }}_{t}$", color='goldenrod', alpha=0.5)
    plt.fill_between(dates[4:], xh1s[4:]-2*stdP[4:], xh1s[4:]+2*stdP[4:],label="95% CI $( \pm 2 \sqrt{\hat{p}^T_t})$ of EBM-KF GMST state $\^{T }_{t}$", color=colorstate, alpha=0.5)

    plt.plot(dates,temps,'o',label='GMST HadCRUT5 measurements $Y_{t}$',markersize=2,color=colorgrey)
    

    #tvar_xhat=ekf_run(temps,n_iters,True)
    #xh1stv=tvar_xhat[:,0]
   # stdPtv=np.sqrt(np.abs(P))[:,0]
   # plt.plot(dates,xh1stv,'-',label='EBM-KF-TVR estimate $ \\breve {T }_{n}}$',color="darkgreen")
   # plt.fill_between(dates, xh1stv-stdPtv, xh1stv+stdPtv,label="68% CI of EBM-KF-TVR state $\sqrt{ \\breve {P}_{n}}$", color="skyblue", alpha=0.5)
    
    ax=plt.gca()
    ax.set_ylabel('Temperature (K)')
    ax.set_xlabel('Year')
    ax.set_yticks(np.arange(286.4,288.2,0.2))
    plot_boilerplate(ax)
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [1,0,4,3,2]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],prop={'size': 10.5})
    axb = ax.twinx()
    mn, mx = ax.get_ylim()
    axb.set_yticks(np.arange(13.2,15,0.2))
    axb.set_ylim(mn- 273.15, mx- 273.15)
    axb.set_ylabel('Temperature (°C)')
    
##    ax1.plot(xnorm,ynorm, color=pcolor)
##    ax1.set_xlabel("Model Predictive Standard Deviations")
##    ax1.set_ylabel("Probability Density")
  #  plt.savefig("RTS_smoother_changes.pdf", format = "pdf")
  #  plt.savefig("RTS_smoother_changes.png", dpi=400,format="png")

    



    def statsnp(seq):
        return [np.mean(seq), np.std(seq), np.min(seq), np.max(seq)]

    
    if(True):
        adjust = xh1s[5:]-xh0s[5:]
        print("Adjust"+str(statsnp(adjust)))
        blind_change=xblind[1:,0]-xblind[:-1,0]
        print("Blind_Change"+str(statsnp(blind_change)))
        cum_change=xh1s-xblind[:,0]
        print("Accumulated_Change"+str(statsnp(cum_change)))

        adjust = xh1d[5:]-xh0d[5:]
        print("Adjust D"+str(statsnp(adjust*zJ_from_W)))
        blind_change=xblind[1:,1]-xblind[:-1,1]
        print("Blind_Change D"+str(statsnp(blind_change*zJ_from_W)))
        cum_change=xh1d-xblind[:,1]
        print("Accumulated_Change D"+str(statsnp(cum_change*zJ_from_W)))
        print(np.argmax(cum_change[5:])+1850)
        

    fig = plt.figure()
    xh1fdbk=(this_xhat[:,3].T[0] - 1)*  (sw_in*a_refl*g_refl)
    plt.plot(dates[30:],xh1fdbk[30:])
    stdPfdbk=np.sqrt(np.abs(P))[:,3,3]*  (sw_in*a_refl*g_refl)
    plt.fill_between(dates[30:],xh1fdbk[30:] + 2*stdPfdbk[30:], xh1fdbk[30:] - 2*stdPfdbk[30:], color = 'red')

##    xh1gad=this_xhat[:,4].T[0] 
##    plt.plot(dates[30:],xh1gad[30:],color='purple')
##    stdPgad=np.sqrt(np.abs(P))[:,4,4]
##    plt.fill_between(dates[30:],xh1gad[30:] + 2*stdPgad[30:], xh1gad[30:] - 2*stdPgad[30:], color = 'blue')
    
##    plt.axhline(y=fdbkA,color='k')
##    fdbkA2 = 0.42 - 0.12
##    plt.axhline(y=fdbkA2,color='g')
##    lnenw=len(dates[30:])
##    plt.fill_between(dates[30:],[fdbkA+fdbkA_sigma]*lnenw,[fdbkA-fdbkA_sigma]*lnenw,color='grey',alpha=0.2)
   # plt.gca().set_ylim([0,0.75])

    #show the modified OHCA stuff
    #ax4.plot(dates,)
    plt.figure()
    xh1TOA=this_xhat[30:,2].T[0]
    stdPtoa=np.sqrt(np.abs(P))[30:,2,2]
    stdStoa=np.sqrt(np.abs(S))[30:,2,2]
    plt.plot(dates[30:],xh1TOA,'-',label='posterior OHCA EBM-KF-uf state estimate $\^{H }_{t}$', color=colorekf)
    plt.fill_between(dates[30:], (xh1TOA-2*stdStoa), (xh1TOA+2*stdStoa),label="95% CI ($\pm 2\sqrt{\hat{s}^H_t}$) of forecast OHCA $\^{H }_{t|t-1}$", color=coloruncert)
    plt.fill_between(dates[30:], (xh1TOA-2*stdPtoa), (xh1TOA+2*stdPtoa),label="95% CI ($\pm 2\sqrt{\hat{p}^H_t}$) of OHCA state $\^{H }_{t}$", color=colorstate)
    #TOA_meas_artif=z[30:,2].T[0]
    plt.plot(dates, TOA_meas_artif,'o',markersize=2,color=colorgrey)
    plt.plot(dates, TOA_swout2,'--',color='darkgoldenrod')
            #plt.plot(dates, TOA_swout,'o',markersize=2,color='red')
    plt.fill_between(dates, TOA_meas_artif-2*np.sqrt(TOA_meas_artif_var), TOA_meas_artif+2*np.sqrt(TOA_meas_artif_var),label="associated 95% uncertainty of $\Psi_{t}$", color="grey",zorder=5,alpha=0.2,lw=0)

    



    
    plt.show()



