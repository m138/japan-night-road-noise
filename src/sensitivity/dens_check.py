"""Does the stability of municipal rankings only reflect population density?
Municipal gross density and population-weighted density (from census small areas); Spearman correlation of mean Lnight with density;
partial Spearman (ranks residualised on log density, loess-free linear) between the main analysis and each sensitivity scenario;
median within-density-quintile Spearman. usage: dens_check.py <sens dir> <estat shp dir> <out dir>"""
import pandas as pd, numpy as np, glob, sys, geopandas as gpd
sys.argv=[sys.argv[0]]+sys.argv[1:]; D,SHP,OUT=sys.argv[1:4]
sys.path.insert(0,D); import importlib.util
spec=importlib.util.spec_from_file_location('ag',f'{D}/agg_sens.py')
# reuse functions from agg_sens without running it: copy the needed code
def load(s):
    X=pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{D}/{s}/rec_*.parquet'))],ignore_index=True)
    X['pref']=X.PREF.astype(str).str.zfill(2); X['muni']=X.pref+X.CITY.astype(str).str.zfill(3)
    allon=X.groupby('KEY_CODE').oncw.transform('all'); tot=X.groupby('KEY_CODE').w.sum(); X=X[~X.oncw|allon].copy()
    X['w']=X.KEY_CODE.map(tot)/X.KEY_CODE.map(X.groupby('KEY_CODE').size()); return X.reset_index(drop=True)
dB=lambda E: np.where(E>0,10*np.log10(np.maximum(E,1e-30)),np.nan)
def muni(X,L,floor=40):
    Lf=np.clip(np.nan_to_num(L,nan=floor),floor,None)
    g=pd.DataFrame({'m':X.muni,'w':X.w,'wl':X.w*Lf,'w60':X.w*(L>=60)}).groupby('m').sum(); return g.wl/g.w, g.w60/g.w*100
X=load('s0'); L=dB(X.E200); base,base60=muni(X,L)
sc={'Floor 30 dB':muni(X,L,30)[0],'Floor 35 dB':muni(X,L,35)[0],'Floor 45 dB':muni(X,L,45)[0],'Cutoff 500 m':muni(X,dB(X.E500))[0]}
for s in ('s1','s2','s3'):
    Y=load(s); sc[f'Seed {s[1:]}']=muni(Y,dB(Y.E200))[0]
# density
G=pd.concat([gpd.read_file(f,ignore_geometry=True)[['KEY_CODE','PREF','CITY','AREA','JINKO']] for f in sorted(glob.glob(f'{SHP}/r2ka*.shp'))])
G['muni']=G.PREF.astype(str).str.zfill(2)+G.CITY.astype(str).str.zfill(3); G=G[G.AREA>0]
G['d']=G.JINKO/(G.AREA/1e6)
M=G.groupby('muni').agg(pop=('JINKO','sum'),area=('AREA','sum')); M['gross']=M['pop']/(M.area/1e6)
M['pwd']=(G.JINKO*G.d).groupby(G.muni).sum()/M['pop']
idx=base.index.intersection(M.index[M['pop']>0]); M=M.loc[idx]; b=base.loc[idx]
r=lambda a,c: a.rank().corr(c.rank())
def resid(y,x):
    yr,xr=y.rank(),np.log(x).rank(); beta=np.polyfit(xr,yr,1); return yr-np.polyval(beta,xr)
rows=[{'item':'n_municipalities','value':len(idx)},
      {'item':'rho mean Lnight vs gross density','value':r(b,M.gross)},
      {'item':'rho mean Lnight vs population-weighted density','value':r(b,M.pwd)},
      {'item':'rho share>=60 dB vs population-weighted density','value':r(base60.loc[idx],M.pwd)}]
q=pd.qcut(M.pwd.rank(method='first'),5,labels=False)
for k,v in sc.items():
    v=v.loc[idx]
    rows.append({'item':f'{k}: rho with main','value':r(b,v)})
    rows.append({'item':f'{k}: partial rho with main, adjusted for population-weighted density','value':resid(b,M.pwd).corr(resid(v,M.pwd))})
    rows.append({'item':f'{k}: median rho within density quintiles','value':np.median([r(b[q==i],v[q==i]) for i in range(5)])})
    rows.append({'item':f'{k}: min rho within density quintiles','value':min([r(b[q==i],v[q==i]) for i in range(5)])})
R=pd.DataFrame(rows); import os; os.makedirs(OUT,exist_ok=True); R.to_csv(f'{OUT}/density_check.csv',index=False); print(R.round(3).to_string())
