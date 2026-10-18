"""Sensitivity of equiv_shift.py: before solving for delta, widen member deviations from the
ensemble mean by the 2024/2025 input-ensemble SD ratio in 2024 (0.0927/0.0838), so delta measures
a level shift only. Also the direct (retrospective-excess) estimate of the 2024-vintage current PIVW."""
import numpy as np, pandas as pd, pickle
from scipy.stats import norm; from scipy.optimize import brentq
BK='/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/kalman-rerun-2026-10-13/before/'
e24=np.loadtxt('../Hist_heads_tails/to_2024/sst_pseudo_rpi.csv',delimiter=','); e25=np.loadtxt('../Hist_heads_tails/sst_pseudo.csv',delimiter=',')
ratio=e24[e24[:,0]==2024,1:].std()/e25[e25[:,0]==2024,1:].std(); print('SD ratio 2024 (old/new):',round(ratio,4))
vs=pd.read_csv(BK+'Results2/historical_names_var_scale.csv').set_index('method_name')
sc={'real':np.load(BK+'Results2/comb_real_scales.npy'),'anthro':np.load(BK+'Results2/comb_anthro_scales.npy')}
REAL=['lowess1dt36wnc','GAM_AR1','Kal_flexLin','FaIR_comb_unB','EBMKF_ta4','GWI_tot_CGWL','CGWL10y_sfUKCP']; ANTH=['FaIR_nonat_unB','GWI_anthro_CGWL']
W={'real':np.array([7.5,7.4,7.9,24.6,16.4,25.7,10.5]),'anthro':np.array([69.,31.])}
p15=pd.read_csv('p15_by_method_vintage.csv',header=[0,1],index_col=[0,1]); out=[]
for exp,g,ms in [('rpi_real','real',REAL),('rpi_anthro','anthro',ANTH)]:
    C={m:[] for m in ms}; S={m:[] for m in ms}
    for k in range(100):
        d=pickle.load(open(BK+f'Results/resultshistens_{exp}{k}.pickle','rb'))
        for m in ms:
            t=d[m]['LT_trend']; C[m].append(np.asarray(t[0])[174]); S[m].append(np.asarray(t[1])[174]*vs.loc[m,'best_alter_scale'])
    C={m:np.array(v) for m,v in C.items()}; S={m:np.array(v) for m,v in S.items()}
    st=np.array([C[m] for m in ms]); C['PIVW_comb']=(st*W[g][:,None]).sum(0)/W[g].sum(); S['PIVW_comb']=np.full(100,sc[g][1])
    for m in ms+['PIVW_comb']:
        cw=C[m].mean()+(C[m]-C[m].mean())*ratio
        tgt=p15.loc[(exp,m),('2024','2024')]
        f=lambda x: norm.sf(1.5,loc=cw+x,scale=S[m]).mean()-tgt
        out.append(dict(exp=exp,method=m,shift_spread_adjusted=brentq(f,-0.5,0.5)))
df=pd.DataFrame(out); old=pd.read_csv('equivalent_level_shift.csv')
df['shift_unadjusted']=old.equiv_shift_v24_minus_v25.values; print(df.round(4).to_string(index=False)); df.to_csv('equivalent_level_shift_spread_adjusted.csv',index=False)
# direct estimate: 2024 rate-CSV (retrospective endpoints) minus retrospective excess measured in the 2025 vintage
r24=dict(lowess1dt36wnc=1.3625355,GAM_AR1=1.3827720,Kal_flexLin=1.3839089,FaIR_comb_unB=1.5663461,EBMKF_ta4=1.3497012,GWI_tot_CGWL=1.3195769,CGWL10y_sfUKCP=1.3811250)
excess=dict(FaIR_comb_unB=1.5517415-1.3676957, CGWL10y_sfUKCP=1.3863918-1.3978902)  # legacy retro minus current, 2025 vintage (git 51de5173 vs f36ba84b)
w=W['real']/W['real'].sum(); cur24=np.array([r24[m]-excess.get(m,0) for m in REAL])
print('\nPIVW realised 2024 vintage: reported(retro)',round(np.sum(w*np.array([r24[m] for m in REAL])),4),' direct current estimate',round(np.sum(w*cur24),4))
an24=dict(FaIR_nonat_unB=1.3465887,GWI_anthro_CGWL=1.3239456); exn=1.315439-1.337859  # changelog retro-before minus current (2025 vintage)
wa=W['anthro']/W['anthro'].sum(); print('PIVW anthro 2024 vintage: reported',round(wa[0]*an24['FaIR_nonat_unB']+wa[1]*an24['GWI_anthro_CGWL'],4),' direct current estimate',round(wa[0]*(an24['FaIR_nonat_unB']-exn)+wa[1]*an24['GWI_anthro_CGWL'],4))
