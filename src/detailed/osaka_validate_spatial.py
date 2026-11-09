"""Osaka: match continuous-monitoring sites (2019-2023) to census sections by route name and distance, and flag sites that do not face a
census road (match == 'none'). Only the match flag, coordinates, height and night LAeq are used later (osaka_local_eval.py); the estimate
columns written here use an older light-vehicle coefficient and are not used in the paper.
Inputs in <data dir>: census27.csv, kasyo27.csv (Road Traffic Census section table for Osaka), osaka_souon_2019_2023.csv (Osaka records
of the continuous monitoring, fiscal years 2019-2023, from data/tn), geo/traffic_census_2021_converted.parquet.
usage: osaka_validate_spatial.py <data dir>  -> <data dir>/osaka_validation_points_spatial.csv"""
import pandas as pd, numpy as np, sys, re, unicodedata as ud
S=sys.argv[1]
# --- night-hour profile within 19-7 (same as lnight.py) ---
h=pd.read_csv(f'{S}/census27.csv',dtype=str); hc=list(h.columns); hrs=hc[11:35]; hnum=[(7+i)%24 for i in range(24)]
NI=[hnum.index(x) for x in [22,23,0,1,2,3,4,5]]; N12=[x for x in range(24) if x>=19 or x<7]
for x in hrs: h[x]=pd.to_numeric(h[x],errors='coerce')
o=h[(h[hc[4]]=='1')&(h[hc[6]].str.strip()=='2')].copy(); o['grp']=np.where(o[hc[2]].isin(['1','2']),'mw','gen')
prof={(g,v):x[[hrs[i] for i in NI]].to_numpy().sum(0)/x[[hrs[hnum.index(t)] for t in N12]].to_numpy().sum() for (g,v),x in o.groupby(['grp',hc[10]])}
def lwa(V,veh,st): return ((46.7 if veh=='1' else 53.2)+30*np.log10(V)) if st else ((82.3 if veh=='1' else 88.8)+10*np.log10(V))
def en(Q,V,veh,st,r): return 10**((lwa(V,veh,st)-8)/10)*np.pi*Q/(1000*V*r)
# --- sections ---
k=pd.read_csv(f'{S}/kasyo27.csv',dtype=str); kc=list(k.columns); num=lambda i: pd.to_numeric(k[kc[i]],errors='coerce')
K=pd.DataFrame({'id':k[kc[0]],'rt':k[kc[3]].str.strip(),'route':k[kc[5]],'muni':k[kc[18]],'motorway':k[kc[19]].str.strip(),'oneway':k[kc[26]].str.strip(),
 'u12s':num(38),'u12l':num(39),'u24s':num(41),'u24l':num(42),'d12s':num(50),'d12l':num(51),'d24s':num(53),'d24l':num(54),
 'roadw':num(112),'carw':num(113),'medw':num(115),'vlim':num(140)})
def city(m): m=str(m); return '27100' if m.startswith('271') and '27102'<=m<='27128' else ('27140' if '27141'<=m<='27147' else m)
K['muni2']=K.muni.map(city)
norm=lambda s: ud.normalize('NFKC',str(s)).replace(' ','')
K['rn']=K.route.map(norm); K['rn0']=K.rn.str.replace(r'[（(].*?[）)]','',regex=True)
def model(r, dist_edge, hgt, f=(1,1)):
    st = r.motorway=='1' or r.rt in ('1','2')
    V = r.vlim if r.vlim and r.vlim>0 else (80 if st else 40); V=float(np.clip(V,40,120) if st else np.clip(V,10,60))
    med=r.medw if r.medw==r.medw else 0
    if r.oneway=='1': ds={'1':dist_edge+r.carw/2,'2':dist_edge+r.carw/2}
    else: w=max((r.carw-med)/2,1.5); ds={'1':dist_edge+w/2,'2':dist_edge+w+med+w/2}
    E=np.zeros(8); g='mw' if st else 'gen'
    for dr,p in (('1','u'),('2','d')):
        rr=np.hypot(ds[dr],hgt)
        for veh,ff in (('1',f[0]),('2',f[1])):
            n12=r[f'{p}24{"s" if veh=="1" else "l"}']-r[f'{p}12{"s" if veh=="1" else "l"}']
            E+=en(max(n12,0) if n12==n12 else 0, V,veh,st,rr)*ff*prof[(g,veh)] if False else en((max(n12,0) if n12==n12 else 0)*ff*prof[(g,veh)],V,veh,st,rr)
    return 10*np.log10(E.mean()) if E.mean()>0 else np.nan

import geopandas as gpd
from shapely.geometry import Point
M=pd.read_csv(f'{S}/osaka_souon_2019_2023.csv',dtype=str)
for c in ['distance_from_roadway','measure_height','night_leq','longitude','latitude']: M[c]=pd.to_numeric(M[c],errors='coerce')
M['rn']=M.road1_name.map(norm); M['rn0']=M.rn.str.replace(r'[（(].*?[）)]','',regex=True)
g=gpd.read_parquet(f'{S}/geo/traffic_census_2021_converted.parquet',columns=['census','geometry'])
g=g[g.census.isin(set(K.id))].to_crs(6674)
Kg=K.set_index('id')
pts=gpd.GeoDataFrame(M,geometry=gpd.points_from_xy(M.longitude,M.latitude),crs=4326).to_crs(6674)
out=[]
for i,m in pts.iterrows():
    near=g[g.distance(m.geometry)<=150].copy()
    if near.empty: out.append({**m.drop('geometry').to_dict(),'match':'none'}); continue
    near['dist']=near.distance(m.geometry); near['rn']=Kg.loc[near.census,'rn'].values; near['rn0']=Kg.loc[near.census,'rn0'].values
    same=near[(near.rn==m.rn)|(near.rn0==m.rn0)]
    lvl='name+near' if not same.empty else 'nearest_only'
    best=(same if not same.empty else near).sort_values('dist').iloc[0]
    r=Kg.loc[best.census]
    if isinstance(r,pd.DataFrame): r=r.iloc[0]
    r=r.copy(); r['id']=best.census
    if pd.isna(r.carw) or pd.isna(r.u24s): out.append({**m.drop('geometry').to_dict(),'match':'nodata'}); continue
    out.append({**m.drop('geometry').to_dict(),'match':lvl,'census':best.census,'gap_m':best.dist,'rt':r.rt,'route_census':r.route,
      'est':model(r,m.distance_from_roadway,m.measure_height),'est_c':model(r,m.distance_from_roadway,m.measure_height,(1.11,1.11))})
R=pd.DataFrame(out); R.to_csv(f'{S}/osaka_validation_points_spatial.csv',index=False)
print(R.match.value_counts().to_dict())
V=R[R.match=='name+near'].copy(); V['err']=V.est-V.night_leq
names={'1':'高速自動車国道','2':'都市高速','3':'一般国道','4':'主要地方道(府道)','5':'主要地方道(市道)','6':'一般府道','7':'指定市一般市道'}
V['種別']=V.rt.map(names)
f=lambda x: pd.Series({'n':len(x),'実測中央値':x.night_leq.median(),'推定中央値':x.est.median(),'誤差_平均':x.err.mean(),'誤差_中央':x.err.median(),'誤差SD':x.err.std(),'MAE':x.err.abs().mean(),'r':np.corrcoef(x.est,x.night_leq)[0,1] if len(x)>2 else np.nan,'±3dB以内':(x.err.abs()<=3).mean()})
pd.set_option('display.width',250)
print('gap_m quantiles', np.percentile(V.gap_m,[50,90,99]).round(1))
print(V.groupby('種別').apply(f).round(2).to_string())
G=V[~V.rt.isin(['1','2'])]; print('一般道'); print(f(G).round(2).to_string())
G2=G[(G.gap_m<=30)&(G.measure_height<=2)]; print('一般道・区間まで30m以内・高さ2m以下'); print(f(G2).round(2).to_string())
T=G[G.err.abs()<=10]; print('|err|<=10:',len(T),round(T.err.mean(),2),round(T.err.std(),2),round(np.corrcoef(T.est,T.night_leq)[0,1],2))
b=np.polyfit(G.est,G.night_leq,1); print('calib measured = %.2f*est + %.2f'%tuple(b))
# within-point uniqueness: same census used by multiple years -> repeated; per unique point
U=G.groupby(['jis_code','noise_measure_point_code']).agg(est=('est','first'),obs=('night_leq','mean')).reset_index()
print('unique points',len(U),'bias',round((U.est-U.obs).mean(),2),'SD',round((U.est-U.obs).std(),2),'r',round(np.corrcoef(U.est,U.obs)[0,1],2))
