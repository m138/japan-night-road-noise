"""Out-of-sample check of the nighttime profile (Section 3.3): profile from one part of the Osaka 24-h sections, evaluated on the rest.
Error per section = 10log10(estimated / observed 22-06 sound energy), with vehicle weights from ASJ coefficients (speed cancels)."""
import pandas as pd, numpy as np, sys
h=pd.read_csv(sys.argv[1],dtype=str); hc=list(h.columns); hrs=hc[11:35]; hnum=[(7+i)%24 for i in range(24)]
NI=[hrs[hnum.index(x)] for x in [22,23,0,1,2,3,4,5]]; N12=[hrs[i] for i in range(24) if hnum[i]>=19 or hnum[i]<7]
for x in hrs: h[x]=pd.to_numeric(h[x],errors='coerce')
o=h[(h[hc[4]]=='1')&(h[hc[6]].str.strip()=='2')].copy(); o['grp']=np.where(o[hc[2]].isin(['1','2']),'mw','gen')
o['sec']=o[hc[0]]+'_'+o[hc[1]]; o['veh']=o[hc[10]]; o['n12']=o[N12].sum(1); o['night']=o[NI].sum(1)
o['frac']=o[NI].div(o[N12].sum(1),axis=0).mean(1)          # per-row mean hourly fraction (as in night_profile)
W={('gen','1'):1,('gen','2'):10**((88.8-82.3)/10),('mw','1'):1,('mw','2'):10**((53.2-45.8)/10)}
o['w']=[W[(g,v)] for g,v in zip(o.grp,o.veh)]
def prof(d): return {(g,v):(x[NI].to_numpy().sum(0)/x[N12].to_numpy().sum()).mean() for (g,v),x in d.groupby(['grp','veh'])}
def err(d,p):
    est=(d.w*d.n12*[p[(g,v)] for g,v in zip(d.grp,d.veh)]).groupby(d.sec).sum(); obs=(d.w*d.night/8).groupby(d.sec).sum()
    return 10*np.log10(est/obs)
secs=o.sec.unique(); print('sections',len(secs))
e_in=err(o,prof(o)); print('in-sample: MAE %.2f median %+.2f'%(e_in.abs().mean(),e_in.median()))
rng=np.random.default_rng(20261109); res=[]
for rep in range(100):
    f=dict(zip(secs,rng.integers(0,2,len(secs)))); k=o.sec.map(f)
    e=pd.concat([err(o[k==i],prof(o[k!=i])) for i in (0,1)]); res.append((e.abs().mean(),e.median()))
r=np.array(res); print('2-fold split by section, 100 repeats: MAE %.2f (range %.2f-%.2f), median %+.2f'%(r[:,0].mean(),r[:,0].min(),r[:,0].max(),r[:,1].mean()))
e=pd.concat([err(o[o.sec==s],prof(o[o.sec!=s])) for s in secs]); print('leave-one-section-out: MAE %.2f median %+.2f'%(e.abs().mean(),e.median()))
