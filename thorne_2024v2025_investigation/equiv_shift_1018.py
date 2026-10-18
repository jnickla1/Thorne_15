"""18 Oct vintage of equiv_shift.py: uniform level shift delta mapping the 18-Oct-vintage 2024 members onto the
2024-vintage mean P(>1.5C) in 2024. sigma_m = saved SE * best_alter_scale (combination: fixed scale).
Checks delta=0 reproduces the saved 18-Oct P (Gaussian approximation)."""
import numpy as np, pandas as pd, pickle
from scipy.stats import norm; from scipy.optimize import brentq
from cfg_1018 import RES, VS as vs, SC as sc, REAL, ANTH, W
p15=pd.read_csv('p15_by_method_vintage_1018.csv',header=[0,1],index_col=[0,1])
out=[]
for exp,g,ms in [('rpi_real','real',REAL),('rpi_anthro','anthro',ANTH)]:
    C={m:[] for m in ms}; S={m:[] for m in ms}
    for k in range(100):
        d=pickle.load(open(RES+f'resultshistens_{exp}{k}.pickle','rb'))
        for m in ms:
            t=d[m]['LT_trend']; C[m].append(np.asarray(t[0])[2024-1850]); S[m].append(np.asarray(t[1])[2024-1850]*vs.loc[m,'best_alter_scale'])
    C={m:np.array(v) for m,v in C.items()}; S={m:np.array(v) for m,v in S.items()}
    st=np.array([C[m] for m in ms]); C['PIVW_comb']=(st*W[g][:,None]).sum(0)/W[g].sum(); S['PIVW_comb']=np.full(100,sc[g][1])
    for m in ms+['PIVW_comb']:
        P=lambda dl: norm.sf(1.5,loc=C[m]+dl,scale=S[m]).mean()
        tgt=p15.loc[(exp,m),('2024','2024')]; saved25=p15.loc[(exp,m),('2025','2024')]
        dl=brentq(lambda x:P(x)-tgt,-0.5,0.5) if tgt>0 else np.nan
        out.append(dict(exp=exp,method=m,mean2024_v25=C[m].mean(),P_gauss_d0=P(0)*100,P_saved_v25=saved25*100,P_saved_v24=tgt*100,equiv_shift_v24_minus_v25=dl))
df=pd.DataFrame(out); print(df.round(4).to_string(index=False)); df.to_csv('equivalent_level_shift_1018.csv',index=False)
