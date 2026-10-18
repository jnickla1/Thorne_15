"""Compare per-method P(current warming > threshold) between the 2024 and 2025 heads/tails vintages."""
import numpy as np, pandas as pd, sys
names = pd.read_csv('snap25/st_rpi_real0.csv')['method_name'].tolist()
rows=[]
for exp in ['rpi_real','rpi_anthro','satcal_real','satcal_anthro']:
    nm = pd.read_csv(f'snap25/st_{exp}0.csv')['method_name'].tolist()
    for vint,last in [('snap24',2024),('snap25',2025)]:
        P = np.array([np.load(f'{vint}/fp_{exp}{m}.npy') for m in range(100)])  # (100, nmeth, nyr, 3)
        for yr in ([2024] if last==2024 else [2024,2025]):
            idx = yr-2005
            for i,n in enumerate(nm):
                rows.append(dict(exp=exp,vintage=last,year=yr,method=n,P15=np.nanmean(P[:,i,idx,2])))
df=pd.DataFrame(rows)
t=df.pivot_table(index=['exp','method'],columns=['vintage','year'],values='P15')
pd.set_option('display.width',200)
print((t*100).round(1))
t.to_csv('p15_by_method_vintage.csv')
