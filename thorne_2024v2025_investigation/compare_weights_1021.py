"""21 Oct vintage of compare_weights.py: PIVW weights from each saved calibration, incl. the 17 Oct 23:58 Results2 var_scale."""
import pandas as pd, numpy as np
from cfg_1021 import REAL, ANTH, VS, HISTSTATS
R='../'; B='/Users/JohnMatthew/Documents/Codex/2026-10-11/hi-chat-this-is-a-more/outputs/thorne-full-rerun-2026-10-17/before/'
out={}
for lab,vs in [('Oct17 before-rerun var_scale',pd.read_csv(B+'Results2/historical_names_var_scale.csv').set_index('method_name')),('Results2 var_scale 17Oct23:58 (21 Oct vintage)',VS)]:
    w=[]
    for grp in [REAL,ANTH]:
        iv=1/vs.loc[grp,'err_var100']; w+=list(100*iv/iv.sum())
    out[(lab,'w%')]=w
for f,lab in [('current_methods_statistics_251221False.csv','Dec2025 (HadCRUT5.0)'),('current_methods_statistics_261010False.csv','Oct10 2026 (HadCRUT5.1)'),(HISTSTATS[3:],'Oct21 2026 (1850-1900)')]:
    d=pd.read_csv(R+f).drop_duplicates('method_name').set_index('method_name'); w=[]
    for grp in [REAL,ANTH]:
        iv=1/d.loc[grp,'RMS']**2; w+=list(100*iv/iv.sum())
    out[(lab,'w%[1/RMS2]')]=w; out[(lab,'RMS')]=list(d.loc[REAL+ANTH,'RMS'])
t=pd.DataFrame(out,index=REAL+ANTH); pd.set_option('display.width',250)
print(t.round(3).to_string()); t.to_csv('weights_by_vintage_1021.csv')
