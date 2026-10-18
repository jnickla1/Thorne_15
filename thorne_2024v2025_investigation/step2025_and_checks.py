"""2025 one-year step per method (heads/tails members), cent20y monotonicity since 1965, and
weight sensitivity of the 2025 PIVW level."""
import numpy as np, pandas as pd, pickle
BK='/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-rerun-2026-10-13/before/'
REAL=['lowess1dt36wnc','GAM_AR1','Kal_flexLin','FaIR_comb_unB','EBMKF_ta4','GWI_tot_CGWL','CGWL10y_sfUKCP']; ANTH=['FaIR_nonat_unB','GWI_anthro_CGWL']
W={'real':np.array([7.5,7.4,7.9,24.6,16.4,25.7,10.5]),'anthro':np.array([69.,31.])}
rows=[]; c20=[]
L={}
for exp,g,ms in [('rpi_real','real',REAL),('rpi_anthro','anthro',ANTH)]:
    for k in range(100):
        d=pickle.load(open(BK+f'Results/resultshistens_{exp}{k}.pickle','rb'))
        v={m:np.asarray(d[m]['LT_trend'][0],float) for m in ms}
        v['PIVW_'+g]=sum(W[g][i]*v[m] for i,m in enumerate(ms))/W[g].sum()
        for m,a in v.items(): rows.append(dict(method=m,member=k,lvl2024=a[174],lvl2025=a[175],step=a[175]-a[174]))
        if g=='real':
            b=np.asarray(d['cent20y']['LT_trend'][2],float); c20.append(b)
            for m in ms: L.setdefault(m,[]).append(v[m][175])
df=pd.DataFrame(rows)
s=df.groupby('method').agg(mean2024=('lvl2024','mean'),mean2025=('lvl2025','mean'),mean_step=('step','mean'),
    frac_members_neg=('step',lambda x:(x<0).mean()),p05_step=('step',lambda x:x.quantile(.05)),min_step=('step','min'))
pd.set_option('display.width',200); print(s.round(4)); s.to_csv('step_2024to2025_headstails.csv')
c20=np.array(c20); yrs=np.arange(1850,2026)
db=np.diff(c20,axis=1); sel=(yrs[1:]>=1966)&np.isfinite(db).all(0)
print('\ncent20y (heads/tails members): last finite year', yrs[np.isfinite(c20).all(0)][-1])
print('fraction of member-years 1966-2015 with cent20y decrease:', (db[:,sel]<0).mean().round(4), ' members with any decrease since 1966:', (db[:,sel]<0).any(1).sum(), ' min step', db[:,sel].min().round(4))
sel2=(yrs[1:]>=1900)&np.isfinite(db).all(0); print('since 1900: frac decreasing', (db[:,sel2]<0).mean().round(3))
# HadCRUT5.1 cent20y
h=pd.read_csv('../Common_Data/HadCRUT.5.1.csv').set_index('Time')['Anomaly']
c=np.array([np.mean([h.loc[y-9:y+9].mean(),h.loc[y-10:y+10].mean()]) for y in range(1860,2016)]); dc=np.diff(c); yy=np.arange(1861,2016)
print('HadCRUT5.1 cent20y decreasing years since 1966:', list(yy[(dc<0)&(yy>=1966)]), ' since 1900:', int(((dc<0)&(yy>=1900)).sum()), 'of', int((yy>=1900).sum()))
# weight sensitivity: paper weights vs 1/RMS^2 from Dec-2025 (HadCRUT5.0) and Oct-2026 (HadCRUT5.1) statistics
lv=np.array([np.mean(L[m]) for m in REAL])
for lab,f in [('paper (Oct-10 var_scale)',None),('1/RMS^2 Dec2025 HadCRUT5.0','../current_methods_statistics_251221False.csv'),('1/RMS^2 Oct2026 HadCRUT5.1','../current_methods_statistics_261010False.csv'),('equal',0)]:
    if f is None: w=W['real']
    elif f==0: w=np.ones(7)
    else: st=pd.read_csv(f).drop_duplicates('method_name').set_index('method_name'); w=1/st.loc[REAL,'RMS'].to_numpy()**2
    w=w/w.sum(); print(f'{lab:32s} PIVW2025={np.sum(w*lv):.4f}  weights%={np.round(100*w,1)}')
