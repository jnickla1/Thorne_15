import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import netCDF4

def colorgen(c,step,const):
    y=(((c+1)%3-1)*step[0], ((c//3+1)%3-1)*step[1], ((c//9+1)%3-1)*step[2])
    color=tuple( max(min(sum(c), 255),0) / 255 for c in zip(y,const))
    return color #


fig, ax = plt.subplots(1,1, figsize=(7, 5))


data = pd.read_csv("../Common_Data/HadCRUT5.csv")
temps_obs = data.loc[:,"Anomaly"].to_numpy()
preind_base = np.mean(temps_obs[0:50]) #preindustrial baseline 1850-1899
temps_obs = temps_obs - preind_base #remove this baseline
years=data.loc[:,"Time"].to_numpy()
#years[0]=1850
nyrs = len(years)

temperature = temps_obs


ax.scatter(years, temperature, color='darkgrey', label='HadCRUT5 observations',s=10, zorder=3)
ax.set_xlabel('Year')
ax.set_ylabel('GMST Anomaly (°C)\nRelative to 1850-1900')
ax.set_xlim([1990,2035])
ax.set_ylim([0.55,2.1])


sdate=1850
#data = np.genfromtxt(open("HadCRUT5.global.annual.csv", "rb"),dtype=float, delimiter=',')
dates=np.arange(sdate,2100)



avg_len_l=9
avg_len_u=1 #actually 9 yrs below, that year, 0 yrs after
empser  = np.full(np.shape(years),np.nan)


WMOoffset = 0.88 # for the WMO data 
#forec = pd.read_csv("GlobalT_WMOLC-ADCPforecast_1991-2020.csv")
filein = netCDF4.Dataset("../Methods/44_EarthModel_CGWL//tasAnom_rcp45_land-prob_global_glb_sample_b8100_1y_ann_18591201-20991130.nc",'r')
comput_temps = np.array(filein.variables['tasAnom'])
    #(240, 3000), 1860 is first year
nsamps = np.shape(comput_temps)[1]
cutoff_n = int(nsamps/30)
comput_temps_baselines = np.mean(comput_temps[0:50,:],axis=0)
comput_temps_align = comput_temps - np.tile(comput_temps_baselines, (np.shape(comput_temps)[0], 1))  
samp_cur = np.full((np.shape(years)[0],cutoff_n) ,np.nan)
kcolors=['firebrick','k','k']
percentiles=[50,20,80]

for i in range(1990-1850, len(years) ):
    chunk=temperature[i-avg_len_l:i+avg_len_u]
 #   chunk_uncert=temps_1std[i-avg_len_l:i+avg_len_u]
    chunk_avg = np.mean(chunk)
    hindc_samps_mean = np.mean(comput_temps_align[i-avg_len_l-10:i+avg_len_u-10,:], axis = 0)
    hindc_samps = comput_temps_align[i-10,:] #-10 because it starts in 1860 not 1850
    forec_samps0 = np.mean(comput_temps_align[(i+1-10):(i+1),:], axis = 0) #10 computed into the future

    #sort_indices = np.argsort(np.abs(hindc_samps-chunk_avg))
    sort_indices = np.argsort(np.abs(hindc_samps-temperature[i]) + np.abs(hindc_samps_mean-chunk_avg))
    forec_samps = forec_samps0[sort_indices[:cutoff_n]] #taking the 1/30th of samples with past 10yr mean closest to observed, so 100
    forec_samps_all = comput_temps_align[(i-10):(i+1),sort_indices[:cutoff_n]] #allyears in this future
    #print(i)
    #print(forec_samps)
    means = np.mean(chunk)/2 + forec_samps/2
    
    #tot_uncert = np.var(chunk) + np.mean(chunk_uncert**2) #this is going into a standard error
    #don't need to increase uncertainty further - taking prior 10 years as having no uncertainty
 #   tot_uncert0 =  np.nanvar(forec_samps)
 #   ses[i] = np.sqrt(tot_uncert0 )/2 #+tot_uncert/ len(chunk) 
 #   samp_cur[i,:] = ( np.mean(chunk)/2 + forec_samps/2 )
    mmedian = [np.percentile(means,50),np.percentile(means,20),np.percentile(means,80)]
    smedian = [np.percentile(forec_samps,50),np.percentile(forec_samps,20),np.percentile(forec_samps,80)]

    for r in range(1,4):
        color=colorgen(r*5,(40,20,40),(52,205,205)) 
        ax.plot(years[i]+.5, mmedian[r-1], '|', color=color, markersize=4)
        if (i ) % ((avg_len_l+avg_len_u)*2) == 14:
            ax.plot([dates[i-avg_len_l], dates[i]+.3], [mmedian[r-1], mmedian[r-1]], color=kcolors[r-1])
            ax.plot([dates[i+1]-.3, dates[i+10]], [mmedian[r-1], mmedian[r-1]], color=color)
            ax.plot([dates[i+1]-.3, dates[i+10]], [mmedian[r-1], mmedian[r-1]],linestyle=(0,(2,5)), color=kcolors[r-1])
            forc_ind = (np.abs(forec_samps-smedian[r-1])).argmin()
            if dates[i]<2015:
                ax.plot(dates[i:i+11],forec_samps_all[:,forc_ind], color=color,
                    linewidth=2,zorder=1, label = f"Selected UKCIP, {percentiles[r-1]:1d}th %ile")
                
            else:
                ax.plot(dates[i:i+11],forec_samps_all[:,forc_ind], color=color,
                    linewidth=2,zorder=1) 
            #ax.plot(dates[i].5, lag_mean, '|', color=color, markersize=6)
        #else:
            #ax.plot(dates[i]-.5, lag_mean, '|', color=color, markersize=2)
ax.plot(0,15, color=kcolors[1],label = "10-year trailing obs mean")
ax.plot(0,15, color=kcolors[1],linestyle=(0,(2,5)),label = "10-y mean of UKCIP member")
ax.plot(0,15, color=kcolors[0],linestyle=(0,(2,5)),label = "10-y mean, median selected")
ax.legend(prop={'size': 11})
ax.set_title("Betts (2023) CGWL Method")
plt.show()
