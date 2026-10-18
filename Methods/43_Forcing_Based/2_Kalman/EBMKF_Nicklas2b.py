from pathlib import Path as _V11Path
import sys as _V11Sys
_V11Sys.path.insert(0, str(_V11Path(__file__).resolve().parent))
from ebmkf_v11_support import (load_inputs as _v11_load_inputs,
    precompute as _v11_precompute, complete_hidden_jacobian as _v11_hidden_jacobian,
    nonnegative_sqrt as _v11_sqrt, check_covariances as _v11_check_covariances,
    covariance_inverse as _v11_covariance_inverse,
    _apply_historical_aerosol_percentile as _v11_apply_aerosol,
    _apply_historical_eco2_percentile as _v11_apply_eco2,
    parameter_choice_density as _v11_choice_density,
    reference_preserving_volcanic_denominator as _v11_volcanic_denominator,
    inverse_volcanic_denominator as _v11_inverse_volcanic)
import numpy as np
# import pykalman
import matplotlib.pyplot as plt
import matplotlib        as mpl
import scipy.stats as stats
from sklearn.metrics import r2_score
#plt.rcParams['text.usetex'] = True
#plt.rcParams["font.family"] = "Times New Roman"
#mpl.rcParams['mathtext.fontset'] = "cm"
import matplotlib.transforms as mtransforms
import pdb

colorekf=(26./255, 44./255, 105./255)
colorstate=(62./255, 207./255, 117./255)
coloruncert=(52./255, 235./255, 235./255)
colorgrey=(0.5,0.5,0.5)




n_iters = len(_v11_load_inputs(__file__))
sz = (n_iters,2) # size of array
sz2d=(n_iters,2,2)

zJ_from_W = 5.1006447*3.154 # global-Earth flux conversion; no extra ocean-area factor
zJtomm=0.121
pcolor='deeppink'
colorrts='goldenrod'
nbins=24
sdate=1850
data = _v11_load_inputs(__file__)
dates=data[:,0]
n_iters = len(data)
dates[0]=sdate


#temps=data[:, 1]+287
lCo2=np.log10(data[:,2])
opt_depth=data[:,3]*0.001 #*0.053) #/1000?0.000100025
aerosol_reference = float(np.mean(data[(1999-1850):(2006-1850),7]))
anthro_clouds = data[:,7] - aerosol_reference

Roc_tvar=np.square(data[:,6])/zJ_from_W/zJ_from_W
tsi=data[:,8]
ocean_heat_measured = data[:,5]
#critexp = np.ones(ln(transmd))

#data2 = np.genfromtxt(open("HadCRUT5.global.annual.csv", "rb"),dtype=float, delimiter=',')
temps=data[:, 1]
R_tvar=np.square(data[:,4]) #still in temperature units



JonesOffset=13.85
offset= -np.mean(temps[(1960-sdate):(1990-sdate)]) +JonesOffset + 273.15 #Jones2013 13.7 to 14 (Jones1999)
#print(offset)
temps=temps+offset
#pindavg= np.mean(temps[(1850-sdate):(1930-sdate)])

# heatCp is derived by precompute_coeffs below (Cs + non-ocean capacity).
# EBMKF_Nicklas2b explicitly retains its separate fitted-capacity override.

sig=5.6704e-8

a_refl= 0.834 #constant clearsky albedo
g_refl=0.909
sw_in=340.2

#Tconst=286.64 #86.7



T02=temps[2002-1850] #~287.55 #in 2002 #(JMN rev 2024)
Teq1850=np.mean(temps[0:25])
pindavg=Teq1850 #286.7 #preindustrial avg
#print("Teq1850",Teq1850)

shaldepth=86 *0.5
deepdepth = 1141 +700 #+1800 # 700

fdbkS = 0.35
fdbkA = 0.40 # starting value copied from author-tuned 2a
fdbkW = 1.3

# Starting forcing settings copied from author-tuned 2a on 2026-10-17.
# Aerosol percentile selects historical AND SSP aerosol curves.
# The eCO2 percentile currently selects observational inputs.
AEROSOL_FORCING_PERCENTILE = 85.0
eCO2_FORCING_PERCENTILE = 50.0
use_observational_forcing_percentiles = True

v11_fitted_heat_capacity = 8.5 # preserve the MPI variant's deliberately fitted capacity
# AR6-reference volcanic sensitivity (magnitude): 20 +/- 5 W/m² per raw AOD.
# The +/-5 documents the supplied uncertainty; this deterministic run uses 20.
volcanic_target_efficiency = 20.0
volcanic_target_efficiency_uncertainty = 5.0

def volcanic_denominator(raw_aod):
    return _v11_volcanic_denominator(globals(), raw_aod)

def raw_aod_from_inverse_denominator(inverse_denominator):
    return _v11_inverse_volcanic(globals(), inverse_denominator)

def precompute_coeffs(printThem):
    if use_observational_forcing_percentiles:
        _v11_apply_aerosol(globals())
        _v11_apply_eco2(globals())
    _v11_precompute(globals(), printThem)
    if printThem and use_observational_forcing_percentiles:
        _v11_choice_density(globals())

        
precompute_coeffs(False)


gad= 1.1 #0.67
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

Q=(covM/30);

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
    inboundd=tsik*inbndf/volcanic_denominator(optdn)*((dfaS+dfaA)+2*(dfaS*dfaA)*(x-T02)+dfaS*anthro_cloud/reference_absorbed_sw)
    outgoingd= outbndf*(4-powp1)*np.power(x, 3-powp1)*(1-B0*lCo2n) #0.0655
# Jacobian of the complete next state: d(ocean temperature)/dT = -Cs/Cd.
# Heat-content increment contains Cs*(T_next-T), so J_HT has no standalone +Cs.
    return np.array([[(1 + (inboundd-outgoingd - gad*(1+Cs/Cd))/heatCp ) , gad/Cd/heatCp], \
                     [ (inboundd-outgoingd-gad*(1+Cs/Cd))*Cs/heatCp +gad*(1+Cs/Cd)  ,  ( 1 - (1-Cs/heatCp)*gad/Cd)]])

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
    inbound=tsik*inbndf*(1+dfaA*(x-T02)+anthro_cloud/reference_absorbed_sw)*(1+dfaS*(x-T02))/volcanic_denominator(optdn)
    outgoing= outbndf*np.power(x, 4-powp1 )*(1-B0*lCo2n)
    oc=(H - Cs*(x-Teq1850))/Cd+oc1850
    if (x)<280:
        raise Exception("Too Cold "+str(k+1850))
    #by convention addition of prior x, oc is not included
    Tchange=(inbound - outgoing - gad*(x-Teq1850 - oc + oc1850) )/heatCp
    return np.array([Tchange, Tchange*Cs + gad*(x-Teq1850 - oc + oc1850) ])


# allocate space for arrays

xhat=np.zeros(sz)      # a posteri estimate of x
P=np.zeros(sz2d)         # a posteri error estimate
F=np.zeros(sz2d)         # state transitions
xhatminus=np.zeros(sz) # a priori estimate of x
Pminus=np.zeros(sz2d)    # a priori error estimate
K=np.zeros(sz2d)         # gain or blending factor

xhathat=np.zeros(sz)   # smoothed a priori estimate of x
Phat=np.zeros(sz2d)      # smoothed posteri error estimate 
Khat=np.zeros(sz2d)      # smoothed gain or blending factor
Shat=np.zeros(sz2d)
S=np.zeros(sz2d)
xblind=np.zeros(sz)


lml=np.zeros(sz)
lsml=0
y=np.zeros(sz)
qqy=np.zeros(sz)

qqyh=[]
qqyk=[]

xnorm = np.linspace(-5.5, 5.5, 200)
qqyker=np.zeros((2,200))
qqykerall=[]
ynorm  = stats.norm.pdf(xnorm)
involcavg = np.mean(1/volcanic_denominator(opt_depth))

def ekf_future(startyr,Pstart,xhatstart,endyr,case,caseA,noVolcs=True,trail_avg_volcs=False): #pass VolcanoFit module as last parameter
    nfiter=endyr-startyr
    Phatf=np.zeros((nfiter,2,2))
    Phatf[0]=Pstart
    xhatf=np.zeros((nfiter,2))
    xhatf[0]=xhatstart
    ##gen volcanoes here
    #
    past_opt_depth=data[:,3]*0.001 #reload past data in case we lost it
    if(noVolcs==True):
        volcs = np.ones(nfiter)*raw_aod_from_inverse_denominator(involcavg) #np.mean(opt_depth) #barely makes a difference
        
    else:
        volcs = noVolcs.genEruptions(nfiter+30)*0.001
        if(trail_avg_volcs):
            allvolcsrec=np.concatenate((past_opt_depth,volcs))
            
            wt_opt_depths = 1/volcanic_denominator(allvolcsrec)
            N = 30 ; cN=int(np.ceil(N/2)) ; fN=int(np.floor(N/2))
            nwt_opt_depths=np.empty(nfiter+30); nwt_opt_depths[:]=involcavg
            for i in range(n_iters,n_iters+nfiter+30-1):
                firsta=i-fN
                nwt_opt_depths[i-n_iters] = (np.sum(wt_opt_depths[(firsta):i+1])+ involcavg*(cN-1))/N
                
                #computing half-average - future is assumed to be the average 
            volcs=raw_aod_from_inverse_denominator(nwt_opt_depths)
               
            
    #print(2*np.sqrt(Phatf[0][0,0]))
    for k in range(1,nfiter):
        xhatf[k]= xhatf[k-1] + compute_update(xhatf[k-1],0,volcs[k-1],case[startyr-2015+k-1],caseA[startyr-2015+k-1]-aerosol_reference)
        #print(xhatf[k], volcs[k-1],case[startyr-2015+k-1],caseA[startyr-2015+k-1])
        Fnow=compute_slope(xhatf[k-1],0,volcs[k-1],case[startyr-2015+k-1],caseA[startyr-2015+k-1]-aerosol_reference)         
        Phatf[k]= np.matmul(np.matmul(Fnow ,Phatf[k-1]), np.transpose(Fnow)) +Q
        #if(k==nfiter-1):
        #    print(Fnow)
        #    print(2*np.sqrt(Phatf[k][0,0]))
    return (xhatf,Phatf)
        

def ekf_run(z,n_iter,retPs=False):
    # intial guesses
    process_cov_history = np.zeros_like(Pminus[:n_iter])
    qqyker.fill(0)
    qqykerall.clear()
    xhat[0] = [Teq1850,oc_meas[0]]
    xblind[0]= xhat[0]
    P[0,:,:] = 1.0
    P[0,1,1] = 20
    Pminus[0,:,:] = P[0,:,:]


    F[0]=compute_slope(xhat[0,:],0) #necessary for last step of RTS smoother

    for k in range(1,n_iter):
        # time update
        F[k]=compute_slope(xhat[k-1],k)
        

        xhatminus[k] = xhat[k-1] + compute_update(xhat[k-1],k)
        xblind[k]= xblind[k-1] + compute_update(xblind[k-1],k)
        
        process_cov_history[k] = Q
        Pminus[k] = F[k] @ P[k-1] @ F[k].T + process_cov_history[k]

        # measurement update if(Rvary):

        S[k]= Pminus[k] + Q*30 + np.matrix([[R_tvar[k],0],[0,Roc_tvar[k]]])

        K[k] = np.matmul(Pminus[k],_v11_covariance_inverse(S[k]))
        y[k]=z[k]-(xhatminus[k])
        xhat[k] = xhatminus[k]+np.matmul(K[k],y[k])
        P[k] = np.matmul((np.eye(2)- K[k] ),Pminus[k])
        stdevS=_v11_sqrt(np.diag(S[k]))
        qqy[k]=y[k]/stdevS
        gmstpdf=stats.norm.pdf(xnorm, loc=qqy[k,0], scale=np.sqrt(R_tvar[k])/stdevS[0])
        ohcapdf=stats.norm.pdf(xnorm, loc=qqy[k,1], scale=np.sqrt(Roc_tvar[k])/stdevS[1])
        qqykerall.append([gmstpdf,ohcapdf])
        qqyker[0] += gmstpdf/(n_iter-1)
        qqyker[1] += ohcapdf/(n_iter-1)
        #lml[k]= -0.5* (np.transpose(y[k])/S[k]*y[k] + np.log(S[k]) + np.log(2*np.pi)) #need to sort this out later
        


    xhathat[n_iter-1]=xhat[n_iter-1]
    Phat[n_iter-1]=P[n_iter-1]
    xhathat[0]=xhat[0]
    Phat[0]=P[0]

    lsml=0
    ybark=0
    

##    #compute moving averages
##    N = 30
##    cumsum, moving_aves = [0], []
##    for i, x in enumerate(temps, 1):
##        cumsum.append(cumsum[i-1] + x)
##        if i<N/2:
##            moving_aves.append(xhat[0][0])
##        if i>=N:
##            moving_ave = (cumsum[i] - cumsum[i-N])/N
##            #can do stuff with moving_ave here
##            moving_aves.append(moving_ave)
    if(True):
        for ik in range(2,n_iter+1):
        # RTS Smoother
            k=n_iter-ik
        # measurement update
            # Symmetric eigensolver; do not reuse a neighbouring year's gain.
            Khat[k] = (P[k] @ F[k+1].T) @ _v11_covariance_inverse(Pminus[k+1])
            xhathat[k] = xhat[k] + np.matmul(Khat[k], xhathat[k+1]-xhatminus[k+1])
            # Stable Joseph form of the RTS conditional covariance.
            conditional_map = np.eye(P.shape[-1]) - Khat[k] @ F[k+1]
            Phat[k] = (conditional_map @ P[k] @ conditional_map.T
                       + Khat[k] @ process_cov_history[k+1] @ Khat[k].T
                       + Khat[k] @ Phat[k+1] @ Khat[k].T)
            Phat[k] = (Phat[k] + Phat[k].T) / 2
            yrts=z[k]-H *xhathat[k]
#        qqyh.append(float(yrts/np.sqrt(H*Phat[k]*np.transpose(H) +np.matrix([[R_tvar[k],0],[0,Roc_tvar[k]]]))))
            Shat[k]= Phat[k]+ Q*30  +np.matrix([[R_tvar[k],0],[0,Roc_tvar[k]]])
####        if (k<len(moving_aves) and k>N/2):
####            ybark= xhathat[k] -moving_aves[k]
####            lsml=lsml - 0.5*(np.log(np.abs(Phat[k])) + np.log(2*np.pi) + np.transpose(ybark)/np.abs(Phat[k])*ybark)
####            qqyh2.append(float(ybark/np.sqrt(np.abs(Phat[k]))))
##    print(sum(lml))
##    print(lsml)
    _v11_check_covariances(globals(), n_iter)
    if (retPs==True):
        return (xhat[0:n_iter]).copy(), (P[0:n_iter]).copy()
    elif (retPs==2):
        return (xhat[0:n_iter]).copy(), (P[0:n_iter]).copy(), (S[0:n_iter]).copy(), (xhatminus[0:n_iter,0]).copy()
    elif (retPs==3):
       # breakpoint()
        return (xhat[0:n_iter,0]).copy(), (P[0:n_iter,0,0]).copy(), (xhathat[0:n_iter,0]).copy(), (Phat[0:n_iter,0,0]).copy()
    elif (retPs==4):
        return (xhat[0:n_iter,0]).copy(), (P[0:n_iter,0,0]).copy(), (xblind[0:n_iter,0]).copy(), (Phat[0:n_iter,0,0]).copy()
    else:
        return (xhat[0:n_iter]).copy()


observ = np.transpose(np.array([temps,oc_meas]))
#print(observ)
this_xhat=ekf_run(observ,n_iters)
xh1s=this_xhat[:,0]
xh0s=xhatminus[:,0]
stdS=np.sqrt(np.abs(S))[:,0,0]
stdP=np.sqrt(np.abs(P))[:,0,0]
Plastretain=np.abs(P)[-1,:,:]
xhh1s=xhathat[:,0]
xlastretain=this_xhat[-1,:]
stdPh=np.sqrt(np.abs(Phat))[:,0,0]
stdSh=np.sqrt(np.abs(Shat))[:,0,0]

xh1d=this_xhat[:,1]
xh0d=xhatminus[:,1]
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
    parameter_choice_density_report = _v11_choice_density(globals())
    r2 = r2_score(observ[:,0], xblind[:,0])
    print('r2 score for GMST to blind is', r2)
    r2 = r2_score(observ[:,1], xblind[:,1])
    print('r2 score for OCHA to blind is', r2)

    r2 = r2_score(moving_aves[15:-16], xh1s[15:-16]) #
    print('r2 score for 30-mean GMST to EBM-KF is', r2)
    r2 = r2_score(ocean_aves[15:-16], xh1d[15:-16])
    print('r2 score for 30-mean OCHA to EBM-KF is', r2)

    plt.rcParams['figure.figsize'] = (7, 6)
   # plt.figure(1)
    fig, (ax1,ax4)= plt.subplots(2, 1, figsize=(7,10), gridspec_kw={ "hspace": 0.3})
    plot_boilerplate(ax1)
    plot_boilerplate(ax4)
    ax1.plot(dates,temps,'o',label='noisy HadCRUT5 GMST measurements $Y_{t}$',markersize=2,color=colorgrey)
    ax1.fill_between(dates, temps-2*data[:,4], temps+2*data[:,4],label="associated 95% uncertainty of $Y_{t}$", color="lightgrey")
   # ax1.plot(dates,this_xhat[:,0],'-',label='posterior GMST EBM-KF-uf state estimate $\^{T }_{t}}$', color=colorekf)
   # ax1.plot(dates,xblind[:,0],'--',label='blind model $\~{T}_{t+1} = \mathbf{F}( \~{T}_{t}, \~{H}_{t} ; [eCO_2]_t ,AOD_t,AC_{t},(\\frac{1}{4}G_{SC})_t )$',color='darkgoldenrod')
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
    ax4.fill_between(dates, ocean_heat_measured-2*data[:,6], ocean_heat_measured+2*data[:,6],label="associated 95% uncertainty of $ \Psi_{t}$", color="lightgrey")
    #ax4.plot(dates,this_xhat[:,1]*zJ_from_W,'-',label='posterior OHCA EBM-KF-uf state estimate $\^{H }_{t}}$', color=colorekf)
   # ax4.plot(dates,xblind[:,1]*zJ_from_W,'--',label='blind model $ \~{H}_{t+1} = \mathbf{F}( \~{T}_{t}, \~{H}_{t} ; [eCO_2]_t ,AOD_t,AC_{t},(\\frac{1}{4}G_{SC})_t )$',color='darkgoldenrod')
    #ax4.legend(loc="upper left", fontsize="10.5")
    handles, labels = ax4.get_legend_handles_labels()
   # order = [1,2,0,3]
    ax4.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc="upper left",prop={'size': 10.5})
    ax4.set_title('Ocean Total Heat State')
    ax4.set_xlabel('Year')
    ax4.set_ylabel('Heat (ZJ)')
    ax4.set_ylim([-150,600])
    ax4b = ax4.twinx()
    mn, mx = ax4.get_ylim()
    ax4b.set_yticks(np.arange(-1,8,1))
    ax4b.set_ylim(mn*zJtomm/10, mx*zJtomm/10)
    
    ax4b.set_ylabel('Thermosteric Sea Level Rise (cm)')

    axes=(ax1,ax4)
    axlabels=['a)','b)','c)','d)']
    for i in range(2):
        ax=axes[i]
        label=axlabels[i]
          # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-(40)/72, 15/72, fig.dpi_scale_trans) #-20/72, 7/72
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, fontsize='large', va='bottom')

    fig.savefig("weighted_average_pres.pdf",format="pdf")
    fig.savefig("weighted_average_pres.png", dpi=400,format="png")


    fig = plt.figure(figsize=(8.25,6))
    #(H - Cs*(x-Teq1850))/Cd+oc1850
    plt.plot(dates,(oc_meas - Cs*(temps-Teq1850))/Cd+oc1850,'o',markersize=2,label='inferred from both Zanna and HADCRUT5 measurements',color=colorgrey)
    #plt.plot(dates,(this_xhat[:,1]- Cs*(this_xhat[:,0]-Teq1850))/Cd+oc1850,'-', label="posterior ocean cons. temp. EBM-KF state estimate $ \^{\\theta }_{t} $",color=colorekf)
                                                                                                                
    #plt.plot(dates,(xblind[:,1]- Cs*(xblind[:,0]-Teq1850))/Cd+oc1850,'--',label='blind model prediction $\~{\\theta }_{t}$',color='darkgoldenrod')
    plt.legend()
    ax=plt.gca()
    ax.set_ylabel('Temperature (K)')
    mn, mx = ax.get_ylim()
    ax1b = ax.twinx()
    ax1b.set_ylim(mn-273.15, mx-273.15)
    ax1b.set_ylabel('Temperature (°C)')
    ax.set_title("Deep Ocean Conservative Temperature $\\theta$")
    #fig.savefig("DOPT.pdf",format="pdf")
    #fig.savefig("DOPT.png", dpi=400,format="png")   

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
    ax_dict["a"].fill_between(dates, ocean_heat_measured-2*data[:,6], ocean_heat_measured+2*data[:,6],label="associated 95% uncertainty of $\Psi_{t}$", color="grey",zorder=5,alpha=0.2,lw=0)
    #plt.fill_between(dates, xhh1s-std, xhh1s+std, color=(171./255, 245./255, 206./255))
    #ax_dict["a"].legend(fontsize="10.5")
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [0,3,2,1,4]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],prop={'size': 10.5})
    ##plt.savefig("Figure_2.pdf", format = "pdf") #dpi=300,format="png")
    ax_dict["a"].set_ylim([-150,600])
    ax4b = ax_dict["a"].twinx()
    mn, mx = ax_dict["a"].get_ylim()
    ax4b.set_ylim(mn*zJtomm/10, mx*zJtomm/10)
    ax4b.set_yticks(np.arange(-1,8,1))
    ax4b.set_ylabel('Thermosteric Sea Level Rise (cm)')

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

    plt.savefig("mainOHCA.pdf", format = "pdf")
    plt.savefig("mainOHCA.png", dpi=400,format="png")


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
    plot_boilerplate()
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
    plt.savefig("RTS_smoother_changes.pdf", format = "pdf")
    plt.savefig("RTS_smoother_changes.png", dpi=400,format="png")




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
        
    
    plt.show()



