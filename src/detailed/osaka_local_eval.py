"""Osaka: (1) add non-census OSM roads (local roads) to the v3-consistent detailed version; (2) compare with the official
'面的評価' (roadside-area assessment). Same assumptions as the earlier version (noise_local.py, aggregate_local.py, check_eval.py),
except: v3 coefficients/receptors (detailed.py), Kinki regional general-road bias from v3, UTM coordinates.
usage: osaka_local_eval.py <v3 run dir S> <noise_lib dir> <buildings parquet> <detailed out dir> <pbf> <validation_points_spatial.csv> <out dir>"""
import pandas as pd, numpy as np, geopandas as gpd, shapely, sys, time, os, glob, osmium, shapely.wkb as wkb
S,LIB,BLD,DOUT,PBF,VP,OUT=sys.argv[1:8]; sys.path.insert(0,LIB); import noise_lib as nl
os.makedirs(OUT,exist_ok=True); t0=time.time(); PC='27'; H=4.0
# ---- rebuild receptor geometries exactly as detailed.py ----
A=gpd.read_file(f'{S}/estat_all/shp/r2ka{PC}.shp'); A=A[A.HCODE==8101]
A=A.dissolve(by='KEY_CODE',aggfunc={'CITY':'first','CITY_NAME':'first','JINKO':'sum','SETAI':'sum'}).reset_index()
z=nl.utm_epsg(A.to_crs(4326).union_all().centroid.x); A=A.to_crs(z); tree=shapely.STRtree(A.geometry.values)
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
Bk=B[~B.KEY_CODE.isin(flag)]; recg=list(Bk.geometry.values)
for a in A[A.KEY_CODE.isin(flag)].itertuples():
    x0,y0,x1,y1=a.geometry.bounds; xs,ys=np.meshgrid(np.arange(x0+12.5,x1,25),np.arange(y0+12.5,y1,25))
    p=shapely.points(xs.ravel(),ys.ravel()); p=p[shapely.contains(a.geometry,p)]
    if len(p)==0: p=np.array([a.geometry.representative_point()])
    recg+=list(p)
recg=np.array(recg,dtype=object); cen=shapely.centroid(recg); cx=shapely.get_x(cen); cy=shapely.get_y(cen)
R=pd.read_parquet(f'{DOUT}/receptors_{PC}.parquet'); assert len(R)==len(recg) and (R.kind.values[:len(Bk)]=='bldg').all()
self_idx=np.full(len(R),-1); self_idx[:len(Bk)]=Bk.index.values; otree=shapely.STRtree(B.geometry.values); rtree=shapely.STRtree(recg)
print('receptors rebuilt',len(R),round(time.time()-t0),'s',flush=True)
# ---- OSM roads in the prefecture minus census roads ----
x0,y0,x1,y1=A.to_crs(4326).total_bounds; wf=osmium.geom.WKBFactory(); rows=[]
keepc={'motorway','trunk','primary','secondary','tertiary','unclassified','residential','living_street','motorway_link','trunk_link','primary_link','secondary_link','tertiary_link'}
for o in osmium.FileProcessor(PBF).with_locations().with_filter(osmium.filter.KeyFilter('highway')):
    if not isinstance(o,osmium.osm.Way) or o.tags.get('highway') not in keepc: continue
    try: g=wkb.loads(wf.create_linestring(o),hex=True)
    except Exception: continue
    c=g.centroid
    if x0<c.x<x1 and y0<c.y<y1: rows.append((o.tags.get('highway'),g))
O=gpd.GeoDataFrame(pd.DataFrame(rows,columns=['highway','geometry']),geometry='geometry',crs=4326).to_crs(z)
ci,_=tree.query(shapely.centroid(O.geometry.values),predicate='dwithin',distance=100); O=O.iloc[np.unique(ci)].reset_index(drop=True)
K=nl.load_sections(f'{S}/geo/traffic_census_2021_converted.parquet'); prof=nl.night_profile(f'{S}/census27.csv')
SRC=nl.source_table(K.drop_duplicates('id'),prof).set_index('id'); K=K[K.id.isin(SRC.index[SRC.ok])].dissolve(by='id').reset_index()
Kz=K.to_crs(z); Kz=Kz[Kz.intersects(shapely.box(*A.total_bounds).buffer(300))].copy()
cbuf=shapely.buffer(shapely.simplify(Kz.geometry.values,2.0),20,quad_segs=2); ct=shapely.STRtree(cbuf)
li,bi=ct.query(O.geometry.values,predicate='intersects'); geo=O.geometry.values.copy()
for l,grp in pd.Series(bi).groupby(li): geo[l]=shapely.difference(geo[l],shapely.union_all(cbuf[grp.values]))
O['geometry']=geo; O=O[~O.geometry.is_empty].explode(index_parts=False); O=O[O.geometry.geom_type=='LineString']; O=O[O.geometry.length>=5]
cls={h:'L1' for h in ['motorway','trunk','primary','secondary','tertiary','motorway_link','trunk_link','primary_link','secondary_link','tertiary_link']}
cls.update({'unclassified':'L2','residential':'L3','living_street':'L3'})
O['cls']=O.highway.map(cls); O['geometry']=shapely.simplify(O.geometry.values,2.0); O=O.reset_index(drop=True); O['wid']=np.arange(len(O))
print('local roads km',{c:round(g.geometry.length.sum()/1000) for c,g in O.groupby('cls')},round(time.time()-t0),'s',flush=True)
# ---- unit-flow energies (same assumptions as noise_local.py; bias = v3 Kinki general-road bias) ----
V3=pd.read_csv(f'{S}/national_validation_points.csv',dtype={'pref':str})
KINKI=['24','25','26','27','28','29','30']; BIAS=float(V3[(V3.grp=='gen')&V3.pref.isin(KINKI)].err.mean())
CUT={'L1':150.0,'L2':60.0,'L3':60.0}; HV={'L1':0.08,'L2':0.03,'L3':0.02}; VV={'L1':40.0,'L2':30.0,'L3':30.0}; OFF={'L1':3.0,'L2':1.5,'L3':1.0}
lw=lambda V,veh: (82.3 if veh=='s' else 88.8)+10*np.log10(V)
UNIT={c:((1-HV[c])*10**((lw(VV[c],'s')-8)/10)+HV[c]*10**((lw(VV[c],'l')-8)/10))*(0.08/8)/(1000*VV[c])*10**(-BIAS/10) for c in CUT}
segs=[(w,c,a[0],a[1],b[0],b[1]) for w,c,g in zip(O.wid,O.cls,O.geometry) for xy in [shapely.get_coordinates(g)] for a,b in zip(xy[:-1],xy[1:]) if np.hypot(*(b-a))>0.01]
SG=pd.DataFrame(segs,columns=['wid','cls','ax','ay','bx','by'])
def energy(sg,si,ri,px_,py_,rg,cn,h):
    ax,ay,bx,by=sg.ax.values[si],sg.ay.values[si],sg.bx.values[si],sg.by.values[si]
    ux,uy=bx-ax,by-ay; L=np.hypot(ux,uy); ux/=L; uy/=L
    qx,qy=px_[ri]-ax,py_[ri]-ay; tA=-(qx*ux+qy*uy); tB=L+tA; dc=np.abs(qx*uy-qy*ux)
    if rg is not None:
        lines=shapely.linestrings(np.stack([np.c_[ax,ay],np.c_[bx,by]],axis=1)); dc=dc-(shapely.distance(lines,cn[ri])-shapely.distance(lines,rg[ri]))
    d=np.maximum(dc,0.5); off=sg.cls.map(OFF).values[si]; u=sg.cls.map(UNIT).values[si]; E=np.zeros(len(si))
    for dd in (np.maximum(d-off,1.0),d+off):
        D3=np.hypot(dd,h); E+=0.5*u*(np.arctan(tB/D3)-np.arctan(tA/D3))/D3
    return E
M=pd.read_csv(VP,dtype=str)
for c in ['longitude','latitude','measure_height','night_leq']: M[c]=pd.to_numeric(M[c],errors='coerce')
mp=gpd.GeoSeries(gpd.points_from_xy(M.longitude,M.latitude),crs=4326).to_crs(z); mx,my=mp.x.values,mp.y.values
mres=np.zeros((len(M),3))
for j,c in enumerate(['L1','L2','L3']):
    sg=SG[SG.cls==c].reset_index(drop=True); lines=shapely.linestrings(np.stack([np.c_[sg.ax,sg.ay],np.c_[sg.bx,sg.by]],axis=1))
    ri,si=shapely.STRtree(lines).query(mp.values,predicate='dwithin',distance=CUT[c])
    np.add.at(mres[:,j],ri,energy(sg,si,ri,mx,my,None,None,M.measure_height.fillna(1.2).values[ri]))
M[['uL1','uL2','uL3']]=mres
wgeom=O.set_index('wid').geometry; U={}
for c in ['L1','L2','L3']:
    Eo=np.zeros(len(R)); Eb=np.zeros(len(R)); wids=SG[SG.cls==c].wid.unique()
    for k0 in range(0,len(wids),4000):
        sg=SG[SG.wid.isin(wids[k0:k0+4000])].reset_index(drop=True)
        lines=shapely.linestrings(np.stack([np.c_[sg.ax,sg.ay],np.c_[sg.bx,sg.by]],axis=1))
        si,ri=rtree.query(lines,predicate='dwithin',distance=CUT[c])
        if len(si)==0: continue
        df=pd.DataFrame({'r':ri,'w':sg.wid.values[si],'E':energy(sg,si,ri,cx,cy,recg,cen,H)}).groupby(['r','w'],as_index=False).E.sum()
        sl=shapely.shortest_line(wgeom.loc[df.w.values].values,recg[df.r.values])
        li,oi=otree.query(sl,predicate='crosses'); ok=oi!=self_idx[df.r.values[li]]
        blk=np.zeros(len(df),bool); blk[np.unique(li[ok])]=True
        np.add.at(Eo,df.r.values[~blk],df.E.values[~blk]); np.add.at(Eb,df.r.values[blk],df.E.values[blk])
    U[f'u{c}_open']=Eo; U[f'u{c}_blk']=Eb; print(c,'done',round(time.time()-t0),'s',flush=True)
U=pd.DataFrame(U); U.to_parquet(f'{OUT}/osaka_receptor_local_unit.parquet')
# ---- calibration of L1 flow and exposure (same as aggregate_local.py) ----
log=[]; P_=lambda *a: (print(*a,flush=True),log.append(' '.join(map(str,a))))
P_(f'general-road bias used for local roads (v3 Kinki): {BIAS:.2f} dB')
Q2,Q3=1000.0,300.0; c=M[(M.match=='none')&(M.uL1>0)].copy(); c['q1']=(10**(c.night_leq/10)-c.uL2*Q2-c.uL3*Q3)/c.uL1; c=c[c.q1>0]; Q1=float(np.median(c.q1))
c['est']=10*np.log10(c.uL1*Q1+c.uL2*Q2+c.uL3*Q3); c['err']=c.est-c.night_leq
P_(f'L1 calibration n={len(c)} Q1={Q1:.0f} veh/day IQR={np.percentile(c.q1,25):.0f}-{np.percentile(c.q1,75):.0f}; in-sample mean err {c.err.mean():.2f} SD {c.err.std():.2f}')
pts=c.noise_measure_point_code.astype(str)+'_'+c.jis_code.astype(str); errs=[]
for p in pts.unique():
    q=np.median(c.q1[pts!=p]); x=c[pts==p]; errs+=list(10*np.log10(x.uL1*q+x.uL2*Q2+x.uL3*Q3)-x.night_leq)
P_(f'leave-one-point-out mean err {np.mean(errs):.2f} SD {np.std(errs,ddof=1):.2f}')
P=R['pop'].sum(); res={}
for lab,(a,f2,f3,f1) in {'センサス道路のみ(A=10)':(10,0,0,0),'主:センサス＋生活道路(A=10)':(10,1,1,1),'生活道路交通量×0.5':(10,.5,.5,1),'生活道路交通量×2':(10,2,2,1),
                       'L1交通量×0.5':(10,1,1,.5),'L1交通量×2':(10,1,1,2),'遮蔽A=5':(5,1,1,1),'遮蔽A=15':(15,1,1,1)}.items():
    g=10**(-a/10); Ec=R.E_open+R.E_blk*g; El=sum(q*(U[f'u{k}_open']+U[f'u{k}_blk']*g) for k,q in (('L1',Q1*f1),('L2',Q2*f2),('L3',Q3*f3)))
    E=Ec+El; L=np.where(E>0,10*np.log10(np.maximum(E,1e-30)),np.nan)
    res[lab]={t: R['pop'][L>=t].sum()/P for t in (45,50,55,60,65)}|{'生活道路が主音源':R['pop'][El>Ec].sum()/P}
    if lab.startswith('主'): R['Ln_local']=L; R['local_dom']=El>Ec
D=pd.DataFrame(res).T; D.round(4).to_csv(f'{OUT}/osaka_local_roads_sensitivity_v3.csv',encoding='utf-8-sig'); P_((D*100).round(1).to_string())
# ---- comparison with official roadside-area assessment (same as check_eval.py) ----
def rd(f):
    for e in ('utf-8-sig','cp932'):
        try: return pd.read_csv(f,encoding=e,dtype=str)
        except UnicodeDecodeError: pass
E=pd.concat([rd(f) for f in glob.glob(f'{S}/tn/tn20*/eval*.csv')]); E.columns=[x.strip() for x in E.columns]; E=E[E.pref_code.str.lstrip('0')=='27']
n=lambda col: pd.to_numeric(E[col],errors='coerce').fillna(0)
E['yr']=E.report_year.astype(int); E=E.sort_values('yr').drop_duplicates(['jis_code','eval_section_code'],keep='last')
for pre,lab in (('prox_','近接空間'),('nonprox_','非近接空間')):
    tot=n(pre+'eval_house_num').sum(); nex=(n(pre+'day_below_std_house_num')+n(pre+'wholeday_exceed_std_house_num')).sum(); P_(f'面的評価 {lab}: 戸数 {tot:.0f} 夜間基準超過 {nex/tot:.3f}')
Rb=R[R.kind=='bldg'].copy(); pts=shapely.points(cx[:len(Bk)],cy[:len(Bk)])
G=Kz[(SRC.loc[Kz.id,'grp'].values=='gen')]; gt=shapely.STRtree(G.geometry.values)
ri,gi=gt.query_nearest(pts,max_distance=80,return_distance=False,all_matches=False)
dist=np.full(len(Rb),np.nan); dist[ri]=shapely.distance(pts[ri],G.geometry.values[gi]); half=np.full(len(Rb),np.nan); half[ri]=G.carw.values[gi]/2
Rb['edge']=dist-half
for lab,mk,std in (('近接空間(車道端〜20m)',Rb.edge<=20,65),('非近接空間(20〜50m)',(Rb.edge>20)&(Rb.edge<=50),60)):
    x=Rb[mk]
    for a in (5,10,15): P_(f'モデル {lab} A={a}: 建物数 {len(x)} 夜間>{std}dB 割合 {(x[f"Ln_A{a}"]>std).mean():.3f}')
city=lambda m: '27100' if '27102'<=m<='27128' else ('27140' if '27141'<=m<='27147' else m)
E['muni']=E.jis_code.str[:5].map(city)
e=E.assign(tot=n('prox_eval_house_num')+n('nonprox_eval_house_num'),nex=n('prox_day_below_std_house_num')+n('prox_wholeday_exceed_std_house_num')+n('nonprox_day_below_std_house_num')+n('nonprox_wholeday_exceed_std_house_num')).groupby('muni')[['tot','nex']].sum()
e['off']=e.nex/e.tot; Rb['muni']=Rb.muni_code.map(city); W=Rb[Rb.edge<=50].copy(); W['std']=np.where(W.edge<=20,65,60)
for a in (5,10,15):
    mm=W.assign(ex=W[f'Ln_A{a}']>W['std']).groupby('muni').ex.mean(); cc=e.join(mm.rename('model'),how='inner'); cc=cc[cc.tot>=500]
    P_(f'A={a}: 市町村数 {len(cc)} Spearman {cc.off.corr(cc.model,method="spearman"):.2f} Pearson {cc.off.corr(cc.model):.2f} 公式平均 {cc.off.mean():.3f} モデル平均 {cc.model.mean():.3f}')
    if a==10: cc.to_csv(f'{OUT}/osaka_eval_vs_model_by_municipality_v3.csv',encoding='utf-8-sig')
open(f'{OUT}/osaka_local_eval_v3_log.txt','w').write('\n'.join(log)); print('done',round(time.time()-t0),'s')
