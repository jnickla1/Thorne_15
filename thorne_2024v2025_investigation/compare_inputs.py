"""(a) Compare observational inputs: HadCRUT5.0.2 vs 5.1 and the 2024 vs 2025 heads/tails ensembles."""
import numpy as np, pandas as pd
R='../'
def pi(df): return df.set_index('Time')['Anomaly']-df[(df.Time>=1850)&(df.Time<=1900)].Anomaly.mean()
h50=pi(pd.read_csv(R+'Common_Data/HadCRUT5.csv')); h51=pi(pd.read_csv(R+'Common_Data/HadCRUT.5.1.csv'))
def c20(s,y): return np.mean([s.loc[y-9:y+9].mean(), s.loc[y-10:y+10].mean()])
print('HadCRUT (rel 1850-1900)   5.0.2   5.1   diff')
for lab,f in [('1981-2010 mean',lambda s:s.loc[1981:2010].mean()),('2015-2024 mean',lambda s:s.loc[2015:2024].mean()),
              ('2024 annual',lambda s:s.loc[2024]),('2025 annual',lambda s:s.get(2025,np.nan)),
              ('cent20y 2005',lambda s:c20(s,2005)),('cent20y 2010',lambda s:c20(s,2010)),
              ('lag10 2015-2024',lambda s:s.loc[2015:2024].mean()),('lag10 2016-2025',lambda s:s.loc[2016:2025].mean())]:
    a,b=f(h50),f(h51); print(f'{lab:22s} {a:7.3f} {b:7.3f} {b-a:+7.3f}')
def ens(p):
    d=pd.read_csv(p,header=None,index_col=0); return d
for lab,p24,p25 in [('rpi',R+'Hist_heads_tails/to_2024/sst_pseudo_rpi.csv',R+'Hist_heads_tails/sst_pseudo.csv'),
                    ('satcal',R+'Hist_heads_tails/to_2024/sst_pseudo_satcal.csv',R+'Hist_heads_tails/sst_pseudo_1981-2010.csv')]:
    e24,e25=ens(p24),ens(p25); print(f'\nheads/tails {lab}: shapes {e24.shape} {e25.shape}; years {e24.index[0]}-{e24.index[-1]} vs {e25.index[0]}-{e25.index[-1]}')
    print('quantity                 2024ens(mean,sd)      2025ens(mean,sd)     diff-of-means')
    for q,f in [('1850-1900 mean',lambda e:e.loc[1850:1900].mean()),('1981-2010 mean',lambda e:e.loc[1981:2010].mean()),
                ('2005-2014 mean',lambda e:e.loc[2005:2014].mean()),('2015-2024 mean',lambda e:e.loc[2015:2024].mean()),
                ('2023',lambda e:e.loc[2023]),('2024',lambda e:e.loc[2024]),('2025',lambda e:e.loc[2025] if 2025 in e.index else e.loc[2024]*np.nan),
                ('2016-2025 mean',lambda e:e.loc[2016:2025].mean())]:
        a,b=f(e24),f(e25); print(f'{q:22s} {a.mean():7.3f} ({a.std():.3f})     {b.mean():7.3f} ({b.std():.3f})     {b.mean()-a.mean():+.3f}')
