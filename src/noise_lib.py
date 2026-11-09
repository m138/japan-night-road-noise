"""Shared functions for national nighttime road-noise estimation (census sections, ASJ RTN-Model 2018 simplified).
Set the environment variable ASJ_RTN=2023 to use the ASJ RTN-Model 2023 coefficient for light vehicles, non-steady flow (81.4 instead of 82.3)."""
import pandas as pd, numpy as np, geopandas as gpd, shapely, os
LS_NS=81.4 if os.environ.get('ASJ_RTN')=='2023' else 82.3   # light vehicles, non-steady flow
def night_profile(census_hourly_csv):
    h=pd.read_csv(census_hourly_csv,dtype=str); hc=list(h.columns); hrs=hc[11:35]; hnum=[(7+i)%24 for i in range(24)]
    NI=[hnum.index(x) for x in [22,23,0,1,2,3,4,5]]; N12=[x for x in range(24) if x>=19 or x<7]
    for x in hrs: h[x]=pd.to_numeric(h[x],errors='coerce')
    o=h[(h[hc[4]]=='1')&(h[hc[6]].str.strip()=='2')].copy(); o['grp']=np.where(o[hc[2]].isin(['1','2']),'mw','gen')
    return {(g,v):x[[hrs[i] for i in NI]].to_numpy().sum(0)/x[[hrs[hnum.index(t)] for t in N12]].to_numpy().sum() for (g,v),x in o.groupby(['grp',hc[10]])}
def lwa(V,veh,st): return ((45.8 if veh=='1' else 53.2)+30*np.log10(V)) if st else ((LS_NS if veh=='1' else 88.8)+10*np.log10(V))
def load_sections(gp):
    G=gpd.read_parquet(gp)
    n=lambda c: pd.to_numeric(G[c],errors='coerce')
    K=pd.DataFrame({'id':G.census.astype(str),'rt':G['道路種別'].astype(str).str.strip(),'route':G['路線名'].astype(str),'muni':G['市区町村コード'].astype(str).str.zfill(5),
      'motorway':G['自動車専用道路の別'].astype(str).str.strip(),'oneway':G['一方通行フラグ'].astype(str).str.strip(),
      'u12s':n('上り／昼間１２時間自動車類交通量／小型車（台）'),'u12l':n('上り／昼間１２時間自動車類交通量／大型車（台）'),
      'u24s':n('上り／２４時間自動車類交通量／小型車（台）'),'u24l':n('上り／２４時間自動車類交通量／大型車（台）'),
      'd12s':n('下り／昼間１２時間自動車類交通量／小型車（台）'),'d12l':n('下り／昼間１２時間自動車類交通量／大型車（台）'),
      'd24s':n('下り／２４時間自動車類交通量／小型車（台）'),'d24l':n('下り／２４時間自動車類交通量／大型車（台）'),
      'roadw':n('幅員構成／道路部幅員（ｍ）'),'carw':n('幅員構成／車道部幅員（ｍ）'),'medw':n('幅員構成／中央帯幅員（ｍ）'),'vlim':n('指定最高速度（ｋｍ／ｈ）')})
    K['geometry']=G.geometry.values
    return gpd.GeoDataFrame(K,geometry='geometry',crs=G.crs)
def source_table(K,prof,bias=None):
    """Hour-averaged night energy coefficient per section (S such that E = S*theta/D for each half), lane offset."""
    st=(K.motorway=='1')|K.rt.isin(['1','2']); g=np.where(st,'mw','gen')
    V=K.vlim.where(K.vlim>0); V=np.where(st,np.clip(V.fillna(80),40,120),np.clip(V.fillna(40),10,60))
    S=np.zeros(len(K))
    for p in ('u','d'):
        for veh,s in (('1','s'),('2','l')):
            n12=(K[f'{p}24{s}']-K[f'{p}12{s}']).clip(lower=0).fillna(0).values
            pf=np.where(st,prof[('mw',veh)].mean(),prof[('gen',veh)].mean())   # mean hourly fraction of 19-07 volume within 22-06
            L=np.where(st,(45.8 if veh=='1' else 53.2)+30*np.log10(V),(LS_NS if veh=='1' else 88.8)+10*np.log10(V))
            S+=10**((L-8)/10)*n12*pf/(1000*V)
    if bias is not None: S=S*10**(-np.where(st,bias['mw'],bias['gen'])/10)
    med=K.medw.fillna(0).values; carw=K.carw.values
    off=np.where(K.oneway.values=='1',0.0,med/2+np.maximum((carw-med)/2,1.5)/2)
    ok=~np.isnan(carw)&(K.u24s.notna()|K.d24s.notna()).values
    return pd.DataFrame({'id':K.id.values,'grp':g,'S':S,'off':off,'ok':ok})
def segments(geoms,ids):
    rows=[]
    for cid,geom in zip(ids,geoms):
        for ln in getattr(geom,'geoms',[geom]):
            c=shapely.get_coordinates(ln)
            if len(c)<2: continue
            a,b=c[:-1],c[1:]; m=np.hypot(*(b-a).T)>0.01
            for (ax,ay),(bx,by) in zip(a[m],b[m]): rows.append((cid,ax,ay,bx,by))
    return pd.DataFrame(rows,columns=['id','ax','ay','bx','by'])
def seg_energy(ax,ay,bx,by,px,py,S,off,h):
    ux,uy=bx-ax,by-ay; L=np.hypot(ux,uy); ux/=L; uy/=L
    qx,qy=px-ax,py-ay; tA=-(qx*ux+qy*uy); tB=L+tA; d=np.maximum(np.abs(qx*uy-qy*ux),0.5)
    E=np.zeros(len(ax))
    for dd in (np.maximum(d-off,1.0),d+off):
        D3=np.hypot(dd,h); E+=0.5*S*(np.arctan(tB/D3)-np.arctan(tA/D3))/D3
    return E
def utm_epsg(lon): return 6688+int((lon+180)//6)+1-51
