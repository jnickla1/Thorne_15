"""Heads/tails ensemble: uniform level shift delta that maps the 2025-vintage 2024 members onto the
2024-vintage mean P(>1.5C) in 2024.  sigma_m = saved SE * best_alter_scale (combination: fixed scale).
First checks delta=0 reproduces the saved 2025-vintage P (validates the Gaussian approximation)."""
import numpy as np, pandas as pd, pickle
from scipy.stats import norm; from scipy.optimize import brentq
BK='/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-rerun-2026-10-13/before/'
vs=pd.read_csv(BK+'Results2/historical_names_var_scale.csv').set_index('method_name')
sc={'real':np.load(BK+'Results2/comb_real_scales.npy'),'anthro':np.load(BK+'Results2/comb_anthro_scales.npy')}
REAL=['lowess1dt36wnc','GAM_AR1','Kal_flexLin','FaIR_comb_unB','EBMKF_ta4','GWI_tot_CGWL','CGWL10y_sfUKCP']; ANTH=['FaIR_nonat_unB','GWI_anthro_CGWL']
W={'real':np.array([7.5,7.4,7.9,24.6,16.4,25.7,10.5]),'anthro':np.array([69.,31.])}
p15=pd.read_csv('p15_by_method_vintage.csv',header=[0,1],index_col=[0,1])
out=[]
for exp,g,ms in [('rpi_real','real',REAL),('rpi_anthro','anthro',ANTH)]:
    C={m:[] for m in ms}; S={m:[] for m in ms}
    for k in range(100):
        d=pickle.load(open(BK+f'Results/resultshistens_{exp}{k}.pickle','rb'))
        for m in ms:
            t=d[m]['LT_trend']; C[m].append(np.asarray(t[0])[2024-1850]); S[m].append(np.asarray(t[1])[2024-1850]*vs.loc[m,'best_alter_scale'])
    C={m:np.array(v) for m,v in C.items()}; S={m:np.array(v) for m,v in S.items()}
    st=np.array([C[m] for m in ms]); C['PIVW_comb']=(st*W[g][:,None]).sum(0)/W[g].sum(); S['PIVW_comb']=np.full(100,sc[g][1])
    for m in ms+['PIVW_comb']:
        P=lambda dl: norm.sf(1.5,loc=C[m]+dl,scale=S[m]).mean()
        tgt=p15.loc[(exp,m),('2024','2024')]; saved25=p15.loc[(exp,m),('2025','2024')]
        dl=brentq(lambda x:P(x)-tgt,-0.5,0.5) if tgt>0 else np.nan
        out.append(dict(exp=exp,method=m,mean2024_v25=C[m].mean(),P_gauss_d0=P(0)*100,P_saved_v25=saved25*100,P_saved_v24=tgt*100,equiv_shift_v24_minus_v25=dl))
df=pd.DataFrame(out); print(df.round(4).to_string(index=False)); df.to_csv('equivalent_level_shift.csv',index=False)
