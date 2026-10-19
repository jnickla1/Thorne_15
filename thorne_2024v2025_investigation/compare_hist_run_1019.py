"""19 Oct vintage of compare_hist_run.py: same-year comparison of the single observational run,
HadCRUT5.0.2 run to 2024 (git 93793643) vs HadCRUT5.1 run to 2025 (hist_results.pickle, 19 Oct 11:29).
PIVW rows use the 19-Oct 1/err_var100 weights; *_paperW rows use the 17-Oct paper weights for comparison."""
import pickle, numpy as np, pandas as pd
from cfg_1019 import REAL, ANTH, W, W_PAPER, HIST
h24=pickle.load(open('hist_results_93793643.pickle','rb')); h25=pickle.load(open(HIST,'rb'))
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
base=d.copy()
for lab,WW in [('',W),('_paperW',W_PAPER)]:
    for g,ms in [('real',REAL),('anthro',ANTH)]:
        X=base.loc[ms]; w=pd.DataFrame(np.isfinite(X.to_numpy())*(WW[g]/100)[:,None],index=X.index,columns=X.columns)
        d.loc['PIVW_'+g+lab]=(X.fillna(0)*w).sum()/w.sum()   # renormalise over methods available (2024-run EBMKF_ta4 ends 2023)
pd.set_option('display.width',220); print(d.round(4)); d.to_csv('hadcrut_run_2024v2025_1019.csv')
