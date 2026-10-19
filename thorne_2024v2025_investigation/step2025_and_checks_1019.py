"""19 Oct vintage of step2025_and_checks.py: 2025 one-year step per method (heads/tails members), cent20y monotonicity
since 1966, and weight sensitivity of the 2025 PIVW level."""
import numpy as np, pandas as pd, pickle
from cfg_1019 import RES, REAL, ANTH, W, W_PAPER, HISTSTATS
rows=[]; c20=[]; L={}
for exp,g,ms in [('rpi_real','real',REAL),('rpi_anthro','anthro',ANTH)]:
    for k in range(100):
        d=pickle.load(open(RES+f'resultshistens_{exp}{k}.pickle','rb'))
        v={m:np.asarray(d[m]['LT_trend'][0],float) for m in ms}
        v['PIVW_'+g]=sum(W[g][i]*v[m] for i,m in enumerate(ms))/W[g].sum()
        for m,a in v.items(): rows.append(dict(method=m,member=k,lvl2024=a[174],lvl2025=a[175],step=a[175]-a[174]))
        if g=='real':
            b=np.asarray(d['cent20y']['LT_trend'][2],float); c20.append(b)
            for m in ms: L.setdefault(m,[]).append(v[m][175])
df=pd.DataFrame(rows); df.to_csv('step_2024to2025_members_1019.csv',index=False)
s=df.groupby('method').agg(mean2024=('lvl2024','mean'),mean2025=('lvl2025','mean'),mean_step=('step','mean'),
    frac_members_neg=('step',lambda x:(x<0).mean()),p05_step=('step',lambda x:x.quantile(.05)),min_step=('step','min'))
pd.set_option('display.width',200); print(s.round(4)); s.to_csv('step_2024to2025_headstails_1019.csv')
c20=np.array(c20); yrs=np.arange(1850,2026)
db=np.diff(c20,axis=1); sel=(yrs[1:]>=1966)&np.isfinite(db).all(0)
print('\ncent20y (heads/tails members): last finite year', yrs[np.isfinite(c20).all(0)][-1])
print('fraction of member-years 1966-2015 with cent20y decrease:', (db[:,sel]<0).mean().round(4), ' members with any decrease since 1966:', (db[:,sel]<0).any(1).sum(), ' min step', db[:,sel].min().round(4))
sel2=(yrs[1:]>=1900)&np.isfinite(db).all(0); print('since 1900: frac decreasing', (db[:,sel2]<0).mean().round(3))
h=pd.read_csv('../Common_Data/HadCRUT.5.1.csv').set_index('Time')['Anomaly']
c=np.array([np.mean([h.loc[y-9:y+9].mean(),h.loc[y-10:y+10].mean()]) for y in range(1860,2016)]); dc=np.diff(c); yy=np.arange(1861,2016)
print('HadCRUT5.1 cent20y decreasing years since 1966:', list(yy[(dc<0)&(yy>=1966)]), ' since 1900:', int(((dc<0)&(yy>=1900)).sum()), 'of', int((yy>=1900).sum()))
lv=np.array([np.mean(L[m]) for m in REAL]); out=[]
for lab,f in [('19 Oct var_scale (1/err_var100)','W'),('17 Oct paper weights','P'),('1/RMS^2 Dec2025 HadCRUT5.0','../current_methods_statistics_251221False.csv'),
              ('1/RMS^2 Oct10 HadCRUT5.1 (False)','../current_methods_statistics_261010False.csv'),('1/RMS^2 Oct19 HadCRUT5.1 var_adjusted (True)',HISTSTATS),('equal',0)]:
    if f=='W': w=W['real']
    elif f=='P': w=W_PAPER['real']
    elif f==0: w=np.ones(7)
    else: st=pd.read_csv(f).drop_duplicates('method_name').set_index('method_name'); w=1/st.loc[REAL,'RMS'].to_numpy()**2
    w=w/w.sum(); print(f'{lab:46s} PIVW2025={np.sum(w*lv):.4f}  weights%={np.round(100*w,1)}'); out.append(dict(weights=lab,PIVW2025=np.sum(w*lv),**{m:100*x for m,x in zip(REAL,w)}))
pd.DataFrame(out).to_csv('weight_sensitivity_1019.csv',index=False)
