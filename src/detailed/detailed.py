"""Detailed (building-level) nighttime road noise for one prefecture, consistent with national v3:
same census sections, ASJ RTN-Model 2018 coefficients (45.8), same region/road-type bias correction, 4 m height.
Receptors = residential OSM footprints (most exposed facade = nearest point of footprint), population allocated by
floor area within each 2020 small area; small areas lacking buildings (none, or < 0.05 buildings per household) use a 25 m grid.
Shielding: a receptor-section path crossing another residential footprint is attenuated by A dB (main 10; 5/15 sensitivity).
usage: detailed.py <v3 run dir S> <noise_lib dir> <pref code> <buildings parquet> <out dir>"""
import pandas as pd, numpy as np, geopandas as gpd, shapely, sys, time, os
S,LIB,PC,BLD,OUT=sys.argv[1:6]; sys.path.insert(0,LIB); import noise_lib as nl
os.makedirs(OUT,exist_ok=True); H=4.0; CUT=200.0; t0=time.time()
REG={**{p:'北海道・東北' for p in ['01','02','03','04','05','06','07']},**{p:'関東' for p in ['08','09','10','11','12','13','14']},**{p:'中部' for p in ['15','16','17','18','19','20','21','22','23']},
     **{p:'近畿' for p in ['24','25','26','27','28','29','30']},**{p:'中国・四国' for p in ['31','32','33','34','35','36','37','38','39']},**{p:'九州・沖縄' for p in ['40','41','42','43','44','45','46','47']}}
# ---- sources: identical to noise_national.py (v3) ----
V=pd.read_csv(f'{S}/national_validation_points.csv',dtype={'pref':str})
bias_gen=V[V.grp=='gen'].groupby(V.pref.map(REG)).err.mean().to_dict(); V['rtc']=V.rt.astype(str).str[0].where(V.rt.astype(str).str[0].isin(['1','2']),'other')
bias_mw=V[V.grp=='mw'].groupby('rtc').err.mean().to_dict()
prof=nl.night_profile(f'{S}/census27.csv')
K=nl.load_sections(f'{S}/geo/traffic_census_2021_converted.parquet')
SRC=nl.source_table(K.drop_duplicates('id'),prof).set_index('id')
K=K[K.id.isin(SRC.index[SRC.ok])].dissolve(by='id').reset_index(); K['pref']=K.muni.str[:2]
rtc=K.rt.astype(str).str[0].where(K.rt.astype(str).str[0].isin(['1','2']),'other')
b=np.where(SRC.loc[K.id,'grp'].values=='mw',rtc.map(bias_mw).values,K.pref.map(REG).map(bias_gen).fillna(np.mean(list(bias_gen.values()))).values)
K['S']=SRC.loc[K.id,'S'].values*10**(-b/10); K['off']=SRC.loc[K.id,'off'].values
# ---- small areas ----
A=gpd.read_file(f'{S}/estat_all/shp/r2ka{PC}.shp'); A=A[A.HCODE==8101]
A=A.dissolve(by='KEY_CODE',aggfunc={'CITY':'first','CITY_NAME':'first','JINKO':'sum','SETAI':'sum'}).reset_index()
z=nl.utm_epsg(A.to_crs(4326).union_all().centroid.x); A=A.to_crs(z)
tree=shapely.STRtree(A.geometry.values)
# ---- buildings -> residential receptors ----
nonres={'industrial','warehouse','commercial','retail','office','school','public','roof','garage','garages','greenhouse','church','temple','shrine','hospital',
 'train_station','hangar','parking','factory','kindergarten','university','college','hotel','service','construction','shed','carport','transportation','government',
 'civic','supermarket','sports_hall','stadium','barn','farm_auxiliary','toilets','bridge','ruins','storage_tank','grandstand','fire_station','gatehouse','kiosk','cowshed','power','water_tower','silo'}
Bd=gpd.read_parquet(BLD).to_crs(z)
rp=shapely.point_on_surface(Bd.geometry.values); idx=tree.query(rp,predicate='within')
keep=np.full(len(Bd),-1); keep[idx[0]]=idx[1]; m=keep>=0
Bd=Bd[m].copy(); Bd['KEY_CODE']=A.KEY_CODE.values[keep[m]]; Bd['area']=Bd.geometry.area
Bd['lv']=pd.to_numeric(Bd.levels,errors='coerce').clip(1,60)
res=(~Bd.building.isin(nonres))&Bd.amenity.isna()&Bd.shop.isna()&(Bd.area>=20)&~((Bd.building=='yes')&(Bd.area>5000))
B=Bd[res].reset_index(drop=True)
A['nres']=A.KEY_CODE.map(B.groupby('KEY_CODE').size()).fillna(0); A['ratio']=A.nres/A.SETAI.where(A.SETAI>0)
flag=set(A[(A.JINKO>0)&((A.nres==0)|(A.ratio<0.05))].KEY_CODE)
Bk=B[~B.KEY_CODE.isin(flag)].copy(); Bk['w']=Bk.area*Bk.lv.fillna(1)
Bk['pop']=Bk.w/Bk.KEY_CODE.map(Bk.groupby('KEY_CODE').w.sum())*Bk.KEY_CODE.map(A.set_index('KEY_CODE').JINKO)
recg=list(Bk.geometry.values); rec=pd.DataFrame({'KEY_CODE':Bk.KEY_CODE.values,'pop':Bk['pop'].values,'kind':'bldg'})
grid=[]
for a in A[A.KEY_CODE.isin(flag)].itertuples():
    x0,y0,x1,y1=a.geometry.bounds; xs,ys=np.meshgrid(np.arange(x0+12.5,x1,25),np.arange(y0+12.5,y1,25))
    p=shapely.points(xs.ravel(),ys.ravel()); p=p[shapely.contains(a.geometry,p)]
    if len(p)==0: p=np.array([a.geometry.representative_point()])
    grid+=[(a.KEY_CODE,a.JINKO/len(p),q) for q in p]
if grid:
    gk,gp,gg=zip(*grid); rec=pd.concat([rec,pd.DataFrame({'KEY_CODE':gk,'pop':gp,'kind':'grid'})],ignore_index=True); recg+=list(gg)
recg=np.array(recg,dtype=object); cen=shapely.centroid(recg); cx=shapely.get_x(cen); cy=shapely.get_y(cen)
print(PC,'osm buildings in pref',int(m.sum()),'residential',len(B),'receptors',len(rec),'grid',len(grid),
      'grid pop share',round(sum(g[1] for g in grid)/A.JINKO.sum(),3),round(time.time()-t0),'s',flush=True)
obst=B.geometry.values; otree=shapely.STRtree(obst)
self_idx=np.full(len(rec),-1); self_idx[:len(Bk)]=Bk.index.values      # B index == position (reset_index)
rtree=shapely.STRtree(recg)
# ---- sections near the prefecture ----
Kz=K.to_crs(z); Kz=Kz[Kz.intersects(shapely.box(*A.total_bounds).buffer(CUT+100))].copy()
Kz['geometry']=shapely.simplify(Kz.geometry.values,2.0)
SG=nl.segments(Kz.geometry.values,Kz.id.values); src=Kz.set_index('id'); secgeom=src.geometry
E_open=np.zeros(len(rec)); E_blk=np.zeros(len(rec)); ids=SG.id.unique()
for bi in range(0,len(ids),150):
    sg=SG[SG.id.isin(ids[bi:bi+150])].reset_index(drop=True)
    lines=shapely.linestrings(np.stack([np.c_[sg.ax,sg.ay],np.c_[sg.bx,sg.by]],axis=1))
    si,ri=rtree.query(lines,predicate='dwithin',distance=CUT)
    if len(si)==0: continue
    ax,ay,bx,by=sg.ax.values[si],sg.ay.values[si],sg.bx.values[si],sg.by.values[si]
    ux,uy=bx-ax,by-ay; L=np.hypot(ux,uy); ux/=L; uy/=L
    qx,qy=cx[ri]-ax,cy[ri]-ay; tA=-(qx*ux+qy*uy); tB=L+tA; dc=np.abs(qx*uy-qy*ux)
    d=np.maximum(dc-(shapely.distance(lines[si],cen[ri])-shapely.distance(lines[si],recg[ri])),0.5)   # nearest facade
    Sx=src.S.loc[sg.id.values[si]].values; off=src.off.loc[sg.id.values[si]].values; E=np.zeros(len(si))
    for dd in (np.maximum(d-off,1.0),d+off):
        D3=np.hypot(dd,H); E+=0.5*Sx*(np.arctan(tB/D3)-np.arctan(tA/D3))/D3
    df=pd.DataFrame({'r':ri,'s':sg.id.values[si],'E':E}).groupby(['r','s'],as_index=False).E.sum()
    sl=shapely.shortest_line(secgeom.loc[df.s.values].values,recg[df.r.values])
    li,oi=otree.query(sl,predicate='crosses'); ok=oi!=self_idx[df.r.values[li]]
    blk=np.zeros(len(df),bool); blk[np.unique(li[ok])]=True
    np.add.at(E_open,df.r.values[~blk],df.E.values[~blk]); np.add.at(E_blk,df.r.values[blk],df.E.values[blk])
    if bi%3000==0: print(bi,len(ids),round(time.time()-t0),'s',flush=True)
rec['E_open']=E_open; rec['E_blk']=E_blk
rec=rec.merge(A[['KEY_CODE','CITY','CITY_NAME']],on='KEY_CODE',how='left'); rec['muni_code']=PC+rec.CITY.astype(str).str.zfill(3)
for a in (5,10,15):
    E=rec.E_open+rec.E_blk*10**(-a/10); rec[f'Ln_A{a}']=np.where(E>0,10*np.log10(np.maximum(E,1e-30)),np.nan)
rec.drop(columns='CITY').to_parquet(f'{OUT}/receptors_{PC}.parquet')
def agg(x,col='Ln_A10'):
    p=x['pop'].sum(); l=x[col]
    return pd.Series({'pop':p,**{f'p{t}':x['pop'][l>=t].sum()/p for t in (50,55,60,65)},'Lmean':np.average(l.fillna(40).clip(lower=40),weights=x['pop']),
                      'grid_share':x['pop'][x.kind=='grid'].sum()/p,'p55_A5':x['pop'][x.Ln_A5>=55].sum()/p,'p55_A15':x['pop'][x.Ln_A15>=55].sum()/p})
M=rec.groupby(['muni_code','CITY_NAME']).apply(agg).reset_index(); M.round(4).to_csv(f'{OUT}/detailed_by_municipality_{PC}.csv',index=False,encoding='utf-8-sig')
T=agg(rec); print(PC,'pref totals',T.round(3).to_dict(),'done',round(time.time()-t0),'s',flush=True)
