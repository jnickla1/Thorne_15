"""Diagnostic plot of the seven user-requested quartic crossing estimates."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]
d=json.loads((root/'data/rpi_real.json').read_text())
years=np.array(d['years']);c10=np.array(d['centered10'],float);c20=np.array(d['centered20'],float)
s=d['definitive'][1]
fig,axes=plt.subplots(2,4,figsize=(15,7.4),sharex=True,sharey=True)
for ax,member in zip(axes.ravel(),s['estimated_members']):
    ix=np.flatnonzero(np.isfinite(c10[member]))[-25:]
    q=np.polynomial.Polynomial.fit(years[ix],c10[member,ix],4)
    xx=np.linspace(1996,2023,250)
    keep=years>=1996
    ax.plot(years[keep],c10[member,keep],color='#616975',lw=1.4,label='Centered 10-year mean')
    ax.plot(years[keep],c20[member,keep],color='#1731bb',lw=2.4,label='Centered 20-year mean')
    ax.plot(xx,q(xx),color='#78b6e9',lw=1.8,label='Quartic fit to 1996–2020')
    cr=s['crossings'][member]
    ax.axhline(1,color='#a8acb2',ls=':',lw=1)
    ax.axvline(2015,color='#d1d4da',ls='--',lw=.8)
    ax.scatter([cr],[1],s=30,color='#4a93d0',zorder=6)
    ax.text(cr,1.035,f'{cr:.2f}',color='#287ebd',ha='center',fontsize=10)
    ax.set_title(f'Member {member}',loc='left',fontsize=12)
    ax.set_xlim(1995,2023);ax.set_ylim(.38,1.35);ax.grid(alpha=.14)
    ax.spines[['top','right']].set_visible(False)
    ax.set_xticks([2000,2010,2020]);ax.tick_params(labelsize=9)
    ax.set_xlabel('Year',fontsize=10)
axes[0,0].set_ylabel('Warming above 1850–1900 (°C)',fontsize=10)
axes[1,0].set_ylabel('Warming above 1850–1900 (°C)',fontsize=10)
axes.ravel()[-1].axis('off')
axes.ravel()[-1].text(.02,.88,'1.0°C crossing benchmark',transform=axes.ravel()[-1].transAxes,fontsize=13,weight='bold')
axes.ravel()[-1].text(.02,.7,'93 observed crossings\n7 quartic estimates\n\nMean over all 100 members: 2011.22\n95% member range: 2004.02–2016.58',transform=axes.ravel()[-1].transAxes,fontsize=11,linespacing=1.6)
fig.suptitle('Quartic estimates for the seven members without an observed centered 20-year crossing',x=.055,ha='left',fontsize=17,weight='bold')
handles,labels=axes[0,0].get_legend_handles_labels()
fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False,fontsize=10,bbox_to_anchor=(.5,.028))
fig.text(.055,.013,'Fits use 25 available centered 10-year values (1996–2020). Dashed line marks the last available centered 20-year year (2015). Member IDs are zero-based.',fontsize=9,color='#5a626c')
fig.subplots_adjust(left=.06,right=.98,top=.89,bottom=.14,hspace=.28,wspace=.18)
fig.savefig(root/'review/quartic_crossing_review.png',dpi=160)
plt.close(fig)
