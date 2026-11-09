"""Summarize the sensitivity runs: receptors, no-road and 40 dB floor shares, floor/cutoff/ASJ2023/on-carriageway/seed effects,
and the effect on municipal rankings of the population-weighted mean Lnight. usage: agg_sens.py <sens dir> <out dir> [v4]
With v4, receptors on carriageways are removed in every scenario (population re-divided within the small area), as in the v4 national model."""
import pandas as pd, numpy as np, glob, sys, os
D,OUT=sys.argv[1:3]; V4=len(sys.argv)>3 and sys.argv[3]=='v4'; os.makedirs(OUT,exist_ok=True)
def load(s): 
    X=pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{D}/{s}/rec_*.parquet'))],ignore_index=True)
    X['pref']=X.PREF.astype(str).str.zfill(2); X['muni']=X.pref+X.CITY.astype(str).str.zfill(3)
    if V4:
        allon=X.groupby('KEY_CODE').oncw.transform('all'); X=X[~X.oncw|allon].copy(); pass
        tot=pd.concat([pd.read_parquet(f,columns=['KEY_CODE','w']) for f in sorted(glob.glob(f'{D}/{s}/rec_*.parquet'))]).groupby('KEY_CODE').w.sum()
        X['w']=X.KEY_CODE.map(tot)/X.KEY_CODE.map(X.groupby('KEY_CODE').size()); X=X.reset_index(drop=True)
    return X
dB=lambda E: np.where(E>0,10*np.log10(np.maximum(E,1e-30)),np.nan)
def ind(X,L,floor=40,w=None):
    w=X.w if w is None else w; W=w.sum(); Lf=np.clip(np.nan_to_num(L,nan=floor),floor,None)
    return {'p55':(w*(L>=55)).sum()/W*100,'p60':(w*(L>=60)).sum()/W*100,'p65':(w*(L>=65)).sum()/W*100,'Lmean':(w*Lf).sum()/W}
def muni(X,L,floor=40,w=None):
    w=X.w if w is None else w; Lf=np.clip(np.nan_to_num(L,nan=floor),floor,None)
    g=pd.DataFrame({'m':X.muni,'w':w,'wl':w*Lf,'w60':w*(L>=60)}).groupby('m').sum(); return g.wl/g.w, g.w60/g.w*100
rows=[]; X=load('s0'); L=dB(X.E200); base_m,base_p=muni(X,L)
def add(name,ind_,mm=None):
    r={'scenario':name,**ind_}
    if mm is not None:
        m,p=mm; r['rho_Lmean_muni']=base_m.rank().corr(m.reindex(base_m.index).rank()); r['max_rank_change_Lmean']=int((base_m.rank(ascending=False)-m.reindex(base_m.index).rank(ascending=False)).abs().max())
        r['rho_p60_muni']=base_p.rank().corr(p.reindex(base_p.index).rank()); r['mean_abs_diff_Lmean_dB']=(m.reindex(base_m.index)-base_m).abs().mean()
    rows.append(r)
add('Main (200 m, floor 40 dB, ASJ 2018, seed 20261106'+(', carriageway receptors removed)' if V4 else ')'),ind(X,L),(base_m,base_p))
for f in (30,35,45): add(f'Floor {f} dB',ind(X,L,f),muni(X,L,f))
L5=dB(X.E500); add('Cutoff 500 m',ind(X,L5),muni(X,L5))
L23=dB(X.E200_23); add('ASJ RTN-Model 2023 (with its own bias correction)',ind(X,L23),muni(X,L23))
if not V4:
    # receptors on a carriageway removed; the population of each small area re-divided among the remaining receptors
    keep=~X.oncw; k=X[keep].groupby('KEY_CODE').size(); tot=X.groupby('KEY_CODE').w.sum()
    w2=np.where(keep&X.KEY_CODE.isin(k.index),X.KEY_CODE.map(tot)/X.KEY_CODE.map(k).fillna(1),0.0); w2=pd.Series(w2,index=X.index)
    w2[~X.KEY_CODE.isin(k.index)]=X.w[~X.KEY_CODE.isin(k.index)]   # small areas where every receptor is on a carriageway: keep as is
    add('Receptors on carriageways removed',ind(X,L,w=w2),muni(X,L,w=w2))
for s in ('s1','s2','s3'):
    if glob.glob(f'{D}/{s}/rec_*.parquet'):
        Y=load(s); LY=dB(Y.E200); add(f'Random placement, seed {s[1:]}',ind(Y,LY),muni(Y,LY))
T=pd.DataFrame(rows); T.round(3).to_csv(f'{OUT}/sensitivity_summary.csv',index=False); print(T.round(3).to_string())
# descriptive facts on receptors
W=X.w.sum(); P=X.groupby('pref').apply(lambda d: pd.Series({'noroad200':d.w[d.E200==0].sum()/d.w.sum()*100,'set40':d.w[~(L[d.index]>=40)].sum()/d.w.sum()*100,
     'noroad500':d.w[d.E500==0].sum()/d.w.sum()*100,'oncw':d.w[d.oncw].sum()/d.w.sum()*100}))
F={'receptors':len(X),'small_areas':X.KEY_CODE.nunique(),'pop':W,'noroad200_%':X.w[X.E200==0].sum()/W*100,'below40_or_noroad_%':X.w[~(L>=40)].sum()/W*100,
   'below40_with_road_%':X.w[(L<40)].sum()/W*100,'noroad500_%':X.w[X.E500==0].sum()/W*100,'oncw_%':X.w[X.oncw].sum()/W*100,'dmin_lt5m_%':X.w[X.dmin<5].sum()/W*100,
   'dmin_lt10m_%':X.w[X.dmin<10].sum()/W*100,'pref_noroad200_min':P.noroad200.min(),'pref_noroad200_max':P.noroad200.max(),'pref_set40_min':P.set40.min(),'pref_set40_max':P.set40.max()}
pd.Series(F).round(3).to_csv(f'{OUT}/receptor_facts.csv',header=['value']); print(pd.Series(F).round(3).to_string()); P.round(2).to_csv(f'{OUT}/receptor_facts_by_pref.csv')
# prefecture ranking of Lmean under floors
pr=lambda f: X.assign(l=np.clip(np.nan_to_num(L,nan=f),f,None)*X.w).groupby('pref').l.sum()/X.groupby('pref').w.sum()
for f in (30,35,45): print('pref rank rho floor',f, round(pr(40).rank().corr(pr(f).rank()),3))
