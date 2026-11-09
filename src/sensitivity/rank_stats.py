"""Median and maximum absolute change in municipal rank (mean Lnight, descending) and number of the top-100 municipalities
that stay in the top 100, for each sensitivity scenario against the v4 main analysis. usage: rank_stats.py <sens dir> <out csv>"""
import pandas as pd, numpy as np, glob, sys
D,OUT=sys.argv[1:3]
def load(s,v4=True):
    X=pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{D}/{s}/rec_*.parquet'))],ignore_index=True)
    X['muni']=X.PREF.astype(str).str.zfill(2)+X.CITY.astype(str).str.zfill(3)
    if v4:
        tot=X.groupby('KEY_CODE').w.sum(); allon=X.groupby('KEY_CODE').oncw.transform('all'); X=X[~X.oncw|allon].copy()
        X['w']=X.KEY_CODE.map(tot)/X.KEY_CODE.map(X.groupby('KEY_CODE').size())
    return X.reset_index(drop=True)
dB=lambda E: np.where(E>0,10*np.log10(np.maximum(E,1e-30)),np.nan)
def muni(X,L,floor=40):
    Lf=np.clip(np.nan_to_num(L,nan=floor),floor,None)
    g=pd.DataFrame({'m':X.muni,'w':X.w,'wl':X.w*Lf}).groupby('m').sum(); return g.wl/g.w
X=load('s0'); L=dB(X.E200); base=muni(X,L)
sc={'Floor 30 dB':muni(X,L,30),'Floor 35 dB':muni(X,L,35),'Floor 45 dB':muni(X,L,45),'Cutoff 500 m':muni(X,dB(X.E500)),
    'ASJ RTN-Model 2023':muni(X,dB(X.E200_23))}
for s in ('s1','s2','s3'):
    Y=load(s); sc[f'Seed {s[1:]}']=muni(Y,dB(Y.E200))
Y=load('s0',v4=False); sc['Carriageway receptors kept']=muni(Y,dB(Y.E200))
rb=base.rank(ascending=False); top=set(rb[rb<=100].index); rows=[]
for k,v in sc.items():
    v=v.reindex(base.index); rv=v.rank(ascending=False); d=(rb-rv).abs()
    rows.append({'scenario':k,'n':len(base),'rho':rb.corr(rv),'median_abs_rank_change':d.median(),'max_abs_rank_change':d.max(),
                 'top100_kept':len(top & set(rv[rv<=100].index))})
R=pd.DataFrame(rows); R.to_csv(OUT,index=False); print(R.round(3).to_string())
