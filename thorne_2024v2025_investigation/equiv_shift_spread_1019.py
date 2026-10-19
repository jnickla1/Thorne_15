"""19 Oct vintage of equiv_shift_spread.py: widen member deviations by the 2024/2025 input-ensemble SD ratio in 2024
before solving for delta. Also the direct (retrospective-excess) estimate of the 2024-vintage current PIVW.
The direct estimate reconstructs the 2024 rate CSV, so it keeps the weights that reproduce its reported 1.406/1.340."""
import numpy as np, pandas as pd, pickle
from scipy.stats import norm; from scipy.optimize import brentq
from cfg_1019 import RES, VS as vs, SC as sc, REAL, ANTH, W, W_PAPER
e24=np.loadtxt('../Hist_heads_tails/to_2024/sst_pseudo_rpi.csv',delimiter=','); e25=np.loadtxt('../Hist_heads_tails/sst_pseudo.csv',delimiter=',')
ratio=e24[e24[:,0]==2024,1:].std()/e25[e25[:,0]==2024,1:].std(); print('SD ratio 2024 (old/new):',round(ratio,4))
p15=pd.read_csv('p15_by_method_vintage_1019.csv',header=[0,1],index_col=[0,1]); out=[]
for exp,g,ms in [('rpi_real','real',REAL),('rpi_anthro','anthro',ANTH)]:
    C={m:[] for m in ms}; S={m:[] for m in ms}
    for k in range(100):
        d=pickle.load(open(RES+f'resultshistens_{exp}{k}.pickle','rb'))
        for m in ms:
            t=d[m]['LT_trend']; C[m].append(np.asarray(t[0])[174]); S[m].append(np.asarray(t[1])[174]*vs.loc[m,'best_alter_scale'])
    C={m:np.array(v) for m,v in C.items()}; S={m:np.array(v) for m,v in S.items()}
    st=np.array([C[m] for m in ms]); C['PIVW_comb']=(st*W[g][:,None]).sum(0)/W[g].sum(); S['PIVW_comb']=np.full(100,sc[g][1])
    for m in ms+['PIVW_comb']:
        cw=C[m].mean()+(C[m]-C[m].mean())*ratio
        tgt=p15.loc[(exp,m),('2024','2024')]
        f=lambda x: norm.sf(1.5,loc=cw+x,scale=S[m]).mean()-tgt
        out.append(dict(exp=exp,method=m,mean2024_v25=C[m].mean(),shift_spread_adjusted=brentq(f,-0.5,0.5),implied_v24_level=C[m].mean()+brentq(f,-0.5,0.5)))
df=pd.DataFrame(out); old=pd.read_csv('equivalent_level_shift_1019.csv')
df['shift_unadjusted']=old.equiv_shift_v24_minus_v25.values; print(df.round(4).to_string(index=False)); df.to_csv('equivalent_level_shift_spread_adjusted_1019.csv',index=False)
# direct estimate (2024 vintage; unchanged inputs): 2024 rate-CSV retrospective endpoints minus the retrospective excess
r24=dict(lowess1dt36wnc=1.3625355,GAM_AR1=1.3827720,Kal_flexLin=1.3839089,FaIR_comb_unB=1.5663461,EBMKF_ta4=1.3497012,GWI_tot_CGWL=1.3195769,CGWL10y_sfUKCP=1.3811250)
excess=dict(FaIR_comb_unB=1.5517415-1.3676957, CGWL10y_sfUKCP=1.3863918-1.3978902)  # FaIR_comb and CGWL current values are bit-identical in the 19-Oct vintage
an24=dict(FaIR_nonat_unB=1.3465887,GWI_anthro_CGWL=1.3239456); exn=1.315439-1.337859
for lab,WW in [('paper (17 Oct)',W_PAPER),('19 Oct var_scale',W)]:
    w=WW['real']/WW['real'].sum(); cur24=np.array([r24[m]-excess.get(m,0) for m in REAL]); wa=WW['anthro']/WW['anthro'].sum()
    print(f'{lab:18s} PIVW realised 2024 vintage: reported(retro) {np.sum(w*np.array([r24[m] for m in REAL])):.4f}  direct current {np.sum(w*cur24):.4f}'
          f' | anthro reported {wa[0]*an24["FaIR_nonat_unB"]+wa[1]*an24["GWI_anthro_CGWL"]:.4f}  direct current {wa[0]*(an24["FaIR_nonat_unB"]-exn)+wa[1]*an24["GWI_anthro_CGWL"]:.4f}')
