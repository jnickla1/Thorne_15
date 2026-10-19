import numpy as np, pandas as pd
d=pd.read_csv('Results2/historical_names_var_scale.csv',index_col=0)
S=np.load('Results2/historical_newsamples.npy'); C=np.load('Results2/historical_central_est.npy')
st=np.load('Results2/historical_standard.npy'); n=st.shape[0]
sl=slice(n-100,n-10)
for j,r in d.iterrows():
    s=S[len(d.index[:list(d.index).index(j)])][sl]; c=C[list(d.index).index(j)][sl]; y=st[sl]
    sd=np.nanstd(s,axis=1); ok=~np.isnan(sd)&(sd>0)
    if ok.sum()<10: continue
    z=(y-c)[ok]/sd[ok]
    lo,hi=np.nanpercentile(s[ok],5,axis=1),np.nanpercentile(s[ok],95,axis=1)
    cov=np.mean((y[ok]>=lo)&(y[ok]<=hi))
    print(f"{r.method_name:24s} scale={r.best_alter_scale:6.2f} rms_z={np.sqrt(np.mean(z**2)):5.2f} cov90={cov:4.2f}")
