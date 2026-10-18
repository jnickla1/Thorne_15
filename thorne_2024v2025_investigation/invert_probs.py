"""Recover per-member current mean (c) and sd (s) from saved exceedance probabilities at two
thresholds:  P_j = Phi((c - t_j)/s)  =>  s = (t2-t1)/(z1-z2),  c = t1 + s*z1.
Validated on the 2025 vintage against the actual pickled current outputs."""
import numpy as np, pandas as pd, pickle
from scipy.stats import norm
BK='/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-rerun-2026-10-13/before/Results/'
REAL=['lowess1dt36wnc','GAM_AR1','Kal_flexLin','FaIR_comb_unB','EBMKF_ta4','GWI_tot_CGWL','CGWL10y_sfUKCP']; ANTH=['FaIR_nonat_unB','GWI_anthro_CGWL']
W={'rpi_real':np.array([7.5,7.4,7.9,24.6,16.4,25.7,10.5]),'rpi_anthro':np.array([69.,31.])}
def invert(P, yr, sth):
    i=yr-2005  # windows: thr1 (1.0) uses 1990.., thr2 (1.5) uses 2005..
    p1=P[:, (yr-1990), 1]; p2=P[:, i, 2]
    z1=norm.ppf(np.clip(p1,1e-300,1-1e-16)); z2=norm.ppf(np.clip(p2,1e-300,1-1e-16))
    s=0.5/(z1-z2); c=(1.0+sth)+s*z1
    return c,s,p1,p2
res=[]
for exp,meths in [('rpi_real',REAL),('rpi_anthro',ANTH)]:
    nm=pd.read_csv(f'snap25/st_{exp}0.csv')['method_name'].tolist()
    for vint,yrs in [('snap24',[2024]),('snap25',[2024,2025])]:
        P=np.array([np.load(f'{vint}/fp_{exp}{m}.npy') for m in range(100)])
        for yr in yrs:
            for name in meths+['PIVW_comb','equal_weight_comb']:
                c,s,p1,p2=invert(P[:,nm.index(name)],yr,0.0)
                ok=np.isfinite(c)&(p1<1-1e-15)&(p2>1e-12)
                truth=np.nan
                if vint=='snap25':
                    tv=[]
                    for m in range(100):
                        d=pickle.load(open(BK+f'resultshistens_{exp}{m}.pickle','rb'))
                        if name in d: tv.append(np.asarray(d[name]['LT_trend'][0])[yr-1850])
                        else:
                            st=np.array([np.asarray(d[k]['LT_trend'][0])[yr-1850] for k in meths]); tv.append(np.sum(st*(W[exp] if name=='PIVW_comb' else 1))/np.sum(W[exp] if name=='PIVW_comb' else np.ones(len(meths))))
                    truth=np.mean(tv)
                res.append(dict(exp=exp,vintage=vint,year=yr,method=name,n_ok=int(ok.sum()),c_inv_mean=np.nanmean(np.where(ok,c,np.nan)),s_inv_median=np.nanmedian(np.where(ok,s,np.nan)),true_mean=truth))
df=pd.DataFrame(res); pd.set_option('display.width',200); print(df.round(4).to_string(index=False)); df.to_csv('inverted_levels.csv',index=False)
