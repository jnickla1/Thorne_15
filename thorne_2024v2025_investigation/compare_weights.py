"""(b) PIVW weights: reconstruct from historical statistics of each vintage and compare."""
import pandas as pd, numpy as np
R='../'; B='/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/thorne-full-rerun-2026-10-17/before/'
REAL=['lowess1dt36wnc','GAM_AR1','Kal_flexLin','FaIR_comb_unB','EBMKF_ta4','GWI_tot_CGWL','CGWL10y_sfUKCP']; ANTH=['FaIR_nonat_unB','GWI_anthro_CGWL']
out={}
vs=pd.read_csv(B+'Results2/historical_names_var_scale.csv').set_index('method_name')
for grp in [REAL,ANTH]:
    iv=1/vs.loc[grp,'err_var100']; out[('2025 var_scale(before)','w')]=out.get(('2025 var_scale(before)','w'),[])+list(iv/iv.sum())
for f,lab in [('current_methods_statistics_251221False.csv','Dec2025 (HadCRUT5.0, to 2024)'),('current_methods_statistics_261010False.csv','Oct2026 (HadCRUT5.1, to 2025)'),('current_methods_statistics_261013True.csv','Oct13 2026 True')]:
    d=pd.read_csv(R+f); d=d.drop_duplicates('method_name').set_index('method_name')
    col='100RMS' if '100RMS' in d.columns else 'RMS'
    w=[]
    for grp in [REAL,ANTH]:
        iv=1/d.loc[grp,col]**2; w+=list(iv/iv.sum())
    out[(lab,'w['+col+']')]=w
    out[(lab,'RMS')]=list(d.loc[REAL+ANTH,'RMS'])
    out[(lab,'best_scale/avg_unc')]=list(d.loc[REAL+ANTH,'avg_unc.(1se)'])
t=pd.DataFrame(out,index=REAL+ANTH); pd.set_option('display.width',250)
print(t.round(4)); t.to_csv('weights_by_vintage.csv')
print('\nvar_scale 2025 (before) rows:'); print(vs.loc[REAL+ANTH,['err_var','err_var100','err_var75','best_alter_scale']])
