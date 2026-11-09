# Fig. S1 (English labels): municipal population shares in three prefectures, national model vs detailed model with buildings
import sys, pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
M,OUT=sys.argv[1],sys.argv[2]
import os; os.makedirs(f'{OUT}/eps',exist_ok=True)
plt.rcParams.update({'ps.fonttype':3,'font.family':'Liberation Sans','font.size':9,'axes.edgecolor':'#52514e'})
m=pd.read_csv(M,dtype={'muni':str,'pref':str})
P=[('27','Osaka','#2a78d6','o'),('28','Hyogo','#d9822b','s'),('30','Wakayama','#3c9a5f','^')]
fig,axs=plt.subplots(1,2,figsize=(7.2,3.7),dpi=200)
for ax,(c,t) in zip(axs,[('60','Lnight ≥ 60 dB'),('55','Lnight ≥ 55 dB')]):
    x,y=m[f'q{c}']*100,m[f's{c}']*100; r=np.corrcoef(x,y)[0,1]; rho=x.rank().corr(y.rank()); lim=max(x.max(),y.max())*1.05
    ax.plot([0,lim],[0,lim],color='#b8b7b2',lw=.8,zorder=1)
    for pc,nm,col,mk in P:
        g=m[m.pref==pc]; ax.scatter(g[f'q{c}']*100,g[f's{c}']*100,s=16,color=col,marker=mk,edgecolor='white',lw=.5,zorder=2,label=f'{nm} (n = {len(g)})')
    ax.set_xlim(0,lim); ax.set_ylim(0,lim); ax.set_aspect('equal'); ax.grid(color='#e5e4e0',lw=.5)
    ax.set_xlabel('Detailed model with buildings (%)'); ax.set_ylabel('National model (%)')
    ax.set_title(f'{t}\nSpearman\u2019s ρ = {rho:.2f}, r = {r:.2f} ({len(m)} municipalities)',loc='left',fontsize=9)
    for s in ['top','right']: ax.spines[s].set_visible(False)
axs[0].legend(frameon=False,fontsize=8,loc='upper left')
fig.tight_layout(); fig.savefig(f'{OUT}/figS1_3pref_simple_vs_detailed.png'); fig.savefig(f'{OUT}/eps/figS1_3pref_simple_vs_detailed.eps'); print('done')
