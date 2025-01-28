import numpy as np
#import pdb;

def get_trend(yr):
    trendstarts= np.array([1850,1880,1960,1980,2020])
    trend_all = np.array([1.10,1.10,1.04,0.76,0])
    
    trend_part = -np.diff(trend_all)
    chunk_lens = np.diff(trendstarts)
    trend_yr = trend_part/chunk_lens
   # 0.89 to 1.32
   # 0.93 to 1.14
   # 0.65 to 0.87
    trend_high = np.array([1.32, 1.32, 1.14, 0.87])
    trend_low = np.array([0.89,  0.89,0.93,0.65]) #repeat first interval - big uncertainty there even though trend is 0
    trend_std  = (trend_high - trend_low)/2/chunk_lens
    one_hot = np.logical_and( yr >= trendstarts[0:-1], yr < trendstarts[1:])
    if(yr<2020):
        i=np.where(one_hot)[0][0]
        return [trend_yr[i],trend_std[i]]
    else:
        return [trend_yr[-1],trend_std[-1]]
    
def run_method(years, temperature, uncert, model_run, experiment_type):
    avg_len_l=9
    avg_len_u=1 #actually 9 yrs below, that year, 0 yrs after
    empser = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    if (experiment_type!="historical"):
        return means, ses, means, ses #forecasts not valude for future runds, return blanks
    
    for i in range(avg_len_l, len(years) - avg_len_u+1):
        chunk=temperature[i-avg_len_l:i+avg_len_u]
        chunk_uncert=temps_1std[i-avg_len_l:i+avg_len_u]
        g = get_trend(i+1850)
        #print(g)
        means[i] = np.mean(chunk) + g[0]*5
        tot_uncert = np.var(chunk)/ len(chunk) + np.mean(chunk_uncert**2) + 10*(g[1])**2 
        ses[i] = np.sqrt(tot_uncert) 
    return means, ses, empser.copy(), empser.copy()

