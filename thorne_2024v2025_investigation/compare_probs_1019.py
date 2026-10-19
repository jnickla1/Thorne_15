"""19 Oct vintage: compare per-method P(current warming > 1.5) between the 2024 (snap24) and 19-Oct (snap25_1019) heads/tails vintages."""
import numpy as np, pandas as pd, sys
names = pd.read_csv('snap25_1019/st_rpi_real0.csv')['method_name'].tolist()
rows=[]
for exp in ['rpi_real','rpi_anthro','satcal_real','satcal_anthro']:
    nm = pd.read_csv(f'snap25_1019/st_{exp}0.csv')['method_name'].tolist()
    for vint,last in [('snap24',2024),('snap25_1019',2025)]:
        P = np.array([np.load(f'{vint}/fp_{exp}{m}.npy') for m in range(100)])  # (100, nmeth, nyr, 3)
        for yr in ([2024] if last==2024 else [2024,2025]):
            idx = yr-2005
            for i,n in enumerate(nm):
                rows.append(dict(exp=exp,vintage=last,year=yr,method=n,P15=np.nanmean(P[:,i,idx,2])))
df=pd.DataFrame(rows)
t=df.pivot_table(index=['exp','method'],columns=['vintage','year'],values='P15')
pd.set_option('display.width',200)
print((t*100).round(1))
t.to_csv('p15_by_method_vintage_1019.csv')
