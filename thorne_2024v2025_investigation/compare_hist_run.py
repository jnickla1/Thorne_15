"""(a)/(c) Same-year comparison of the single HadCRUT observational run:
HadCRUT5.0.2 run to 2024 (git 93793643) vs HadCRUT5.1 run to 2025 (git f36ba84b, paper state)."""
import pickle, numpy as np, pandas as pd
REAL=['lowess1dt36wnc','GAM_AR1','Kal_flexLin','FaIR_comb_unB','EBMKF_ta4','GWI_tot_CGWL','CGWL10y_sfUKCP']; ANTH=['FaIR_nonat_unB','GWI_anthro_CGWL']
W=dict(real=np.array([7.5,7.4,7.9,24.6,16.4,25.7,10.5])/100, anth=np.array([69,31])/100)
h24=pickle.load(open('hist_results_93793643.pickle','rb')); h25=pickle.load(open('hist_results_f36ba84b.pickle','rb'))
def c(h,m,y):
    t=h[m]['LT_trend']; a=np.asarray(t[0],float); return a[y-1850] if y-1850<len(a) else np.nan
def r(h,m,y):
    t=h[m]['LT_trend']; a=np.asarray(t[2],float); return a[y-1850] if y-1850<len(a) else np.nan
rows=[]
for m in REAL+ANTH+['cent20y']:
    rows.append(dict(method=m, v24_cur2024=c(h24,m,2024), v25_cur2024=c(h25,m,2024), v25_cur2025=c(h25,m,2025),
                     v24_retro2024=r(h24,m,2024), v25_cur2015=c(h25,m,2015), v24_cur2015=c(h24,m,2015)))
d=pd.DataFrame(rows).set_index('method')
d['vintage_eff_2024']=d.v25_cur2024-d.v24_cur2024; d['year_step_25']=d.v25_cur2025-d.v25_cur2024; d['net']=d.v25_cur2025-d.v24_cur2024
d['vintage_eff_2015']=d.v25_cur2015-d.v24_cur2015
for g,ms in [('real',REAL),('anth',ANTH)]:
    X=d.loc[ms]; w=pd.DataFrame(np.isfinite(X.to_numpy())*W[g][:,None],index=X.index,columns=X.columns)
    d.loc['PIVW_'+g]=(X.fillna(0)*w).sum()/w.sum()   # renormalise over methods available (2024-run EBMKF_ta4 ends 2023)
pd.set_option('display.width',220); print(d.round(4)); d.to_csv('hadcrut_run_2024v2025.csv')
