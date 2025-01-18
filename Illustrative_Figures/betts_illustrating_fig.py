import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def colorgen(c,step,const):
    y=(((c+1)%3-1)*step[0], ((c//3+1)%3-1)*step[1], ((c//9+1)%3-1)*step[2])
    color=tuple( max(min(sum(t), 255),0) / 255 for t in zip(y,const))
    return color #


fig, ax = plt.subplots(1,1, figsize=(7, 5))

# Generate fictitious temperature data (for demonstration)
data = np.genfromtxt(open("../Common_Data/toyKFmodelData8c.csv", "rb"),dtype=float, delimiter=',')
years=data[:,0]
years[0]=1850
nyrs = len(years)
temperature = data[:, 1]
offset_t = np.mean(temperature[0:51]) #1850-1900
temperature = temperature - offset_t


ax.scatter(years, temperature, color='darkgrey', label='HadCRUT5',s=10, zorder=3)
ax.set_xlabel('Year')
ax.set_ylabel('GMST Anomaly (°C)\nRelative to 1850-1900')
ax.set_xlim([1996-0.2,2040+0.2])
ax.set_ylim([0.6,2.1])

data_simulation = np.genfromtxt(open("../Common_Data/all_nl_cesm2_ts_0000-0089.csv", "rb"),dtype=float, delimiter=',')

twTR=np.where(data_simulation == 0, float("nan"), data_simulation)

[sms, smyrs]=np.shape(twTR)
#dif=0.21
twTRmean=np.nanmean(twTR,axis=0)
offset_s = np.mean(twTRmean[0:51]) #1850-1900
sdate=1850
#data = np.genfromtxt(open("HadCRUT5.global.annual.csv", "rb"),dtype=float, delimiter=',')
dates=np.arange(sdate,sdate+smyrs)
twTR = twTR-offset_s
dists_23 = np.abs(twTR[:,2022-1850] - temperature[2022-1850])
order  = dists_23.argsort()

for r in range(1,4):
    ax.plot(dates[nyrs-2:],twTR[order[r-1],nyrs-2:],'-', color=colorgen(r*5,(80,80,80),(52,205,205)),
            linewidth=2,zorder=1, label = f"ESM Sim. #{r:1d}")


avg_len_l=10
avg_len_u=9 #actually 15 yrs below, that year, 14 yrs after
labeled=False
for i in range(1996-1850, len(years) - avg_len_u):
    lag_mean = np.mean(temperature[i-avg_len_l:i+avg_len_u])
    color = 'black' 
    if (i ) % (avg_len_l) == ((2024-1850) % avg_len_l):
        ax.plot([years[i-avg_len_l], years[i+avg_len_u]], [lag_mean, lag_mean], color=color)
        if not(labeled):
            ax.plot(years[i]-.5, lag_mean, '|', color=color, markersize=6,label='20-yr mean')
            labeled=True
        else:
            ax.plot(years[i]-.5, lag_mean, '|', color=color, markersize=6)
        
    else:
        ax.plot(years[i]-.5, lag_mean, '|', color=color, markersize=2)

for r in range(1,4):
    twTR_cur = twTR[order[r-1],nyrs:]
    color=colorgen(r*5,(80,80,80),(52,205,205))
    proj_temp = np.concatenate((temperature,twTR_cur)) 
    for i in range(len(years) - avg_len_u, 2024-1850):
        lag_mean = np.mean(proj_temp[i-avg_len_l:i+avg_len_u])
        if (i ) % (avg_len_l+avg_len_u-4) == 8:
            ax.plot([dates[i-avg_len_l], dates[i+avg_len_u]], [lag_mean, lag_mean], color=color)
            ax.plot(dates[i]-.5, lag_mean, '|', color=color, markersize=6)
        else:
            ax.plot(dates[i]-.5, lag_mean, '|', color=color, markersize=2)

ax.legend(prop={'size': 14})
ax.set_title("Betts (2023) CGWL")
plt.show()
