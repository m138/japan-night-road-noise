"""National nighttime road noise (census roads only) with population-proportional random receptors (1 per ~10 residents)
uniformly placed in each 2020 census small area; no building shielding; region-specific bias correction from national validation.
v4: receptors that fall on a carriageway (distance to a census section center line < carriageway width/2) are removed and the
population of the small area is divided among the remaining receptors (if all fall on carriageways, they are kept).
usage: noise_national.py <v3 data dir> <noise_lib dir> <out dir> [prefs]"""
import pandas as pd, numpy as np, geopandas as gpd, shapely, sys, glob, os, time
sys.path.insert(0,sys.argv[2]); import noise_lib as nl
S=sys.argv[1]; NAT=sys.argv[3]; PREFS=sys.argv[4].split(',') if len(sys.argv)>4 else None
H=4.0; CUT=200.0; PER=10; rng=np.random.default_rng(20261106)
REG={**{p:'北海道・東北' for p in ['01','02','03','04','05','06','07']},**{p:'関東' for p in ['08','09','10','11','12','13','14']},**{p:'中部' for p in ['15','16','17','18','19','20','21','22','23']},
     **{p:'近畿' for p in ['24','25','26','27','28','29','30']},**{p:'中国・四国' for p in ['31','32','33','34','35','36','37','38','39']},**{p:'九州・沖縄' for p in ['40','41','42','43','44','45','46','47']}}
V=pd.read_csv(f'{S}/national_validation_points.csv',dtype={'pref':str})
bias_gen=V[V.grp=='gen'].groupby(V.pref.map(REG)).err.mean().to_dict(); V['rtc']=V.rt.astype(str).str[0].where(V.rt.astype(str).str[0].isin(['1','2']),'other')
bias_mw=V[V.grp=='mw'].groupby('rtc').err.mean().to_dict()   # v2: road-type-specific motorway correction
print('bias gen by region',{k:round(v,2) for k,v in bias_gen.items()},' mw',{k:round(v,2) for k,v in bias_mw.items()},flush=True)
prof=nl.night_profile(f'{S}/census27.csv')
K=nl.load_sections(f'{S}/geo/traffic_census_2021_converted.parquet')
SRC=nl.source_table(K.drop_duplicates('id'),prof).set_index('id')
K=K[K.id.isin(SRC.index[SRC.ok])].dissolve(by='id').reset_index()
K['pref']=K.muni.str[:2]
rtc=K.rt.astype(str).str[0].where(K.rt.astype(str).str[0].isin(['1','2']),'other')
b=np.where(SRC.loc[K.id,'grp'].values=='mw',rtc.map(bias_mw).values,K.pref.map(REG).map(bias_gen).fillna(np.mean(list(bias_gen.values()))).values)
K['S']=SRC.loc[K.id,'S'].values*10**(-b/10); K['off']=SRC.loc[K.id,'off'].values; K['hw']=K.carw.values/2; K['grp']=SRC.loc[K.id,'grp'].values
K['geometry']=K.geometry.simplify(0.00002)
print('sections',len(K),flush=True)
os.makedirs(NAT,exist_ok=True)
files=sorted(glob.glob(f'{S}/estat_all/shp/r2ka*.shp'))
for f in files:
    pc=os.path.basename(f)[4:6]
    if PREFS and pc not in PREFS: continue
    if os.path.exists(f'{NAT}/sa_{pc}.parquet'): continue
    t0=time.time()
    A=gpd.read_file(f); A=A[A.HCODE==8101][['KEY_CODE','PREF','CITY','JINKO','geometry']]
    A=A.dissolve(by='KEY_CODE',aggfunc={'PREF':'first','CITY':'first','JINKO':'sum'}).reset_index(); A=A[A.JINKO>0]
    z=nl.utm_epsg(A.to_crs(4326).union_all().centroid.x); A=A.to_crs(z)
    # receptors
    n=np.maximum(1,np.round(A.JINKO.values/PER).astype(int)); px=[];py=[];kc=[];pw=[]
    for (key,pop,geom),k in zip(A[['KEY_CODE','JINKO','geometry']].itertuples(index=False),n):
        x0,y0,x1,y1=geom.bounds; got=np.empty((0,2)); tries=0
        prep=shapely.prepared.prep(geom) if hasattr(shapely,'prepared') else None
        while len(got)<k and tries<50:
            m=max(4*(k-len(got)),16); c=np.c_[rng.uniform(x0,x1,m),rng.uniform(y0,y1,m)]
            ins=shapely.contains_xy(geom,c[:,0],c[:,1]); got=np.vstack([got,c[ins]]); tries+=1
        if len(got)==0: rp=geom.representative_point(); got=np.array([[rp.x,rp.y]])
        got=got[:k]; px.append(got[:,0]); py.append(got[:,1]); kc+= [key]*len(got); pw.append(np.full(len(got),pop/len(got)))
    px=np.concatenate(px); py=np.concatenate(py); pw=np.concatenate(pw); kc=np.array(kc)
    # sections near this prefecture
    box=shapely.box(*A.total_bounds).buffer(CUT+100)
    Kz=K.to_crs(z); Kz=Kz[Kz.intersects(box)]
    # v4: remove receptors on a carriageway and re-divide the small-area population among the rest
    pts=shapely.points(px,py); ktree=shapely.STRtree(Kz.geometry.values)
    pi,ki=ktree.query(pts,predicate='dwithin',distance=float(np.nanmax(Kz.hw.values)))
    on=np.zeros(len(pts),bool); hit=shapely.distance(pts[pi],Kz.geometry.values[ki])<Kz.hw.values[ki]; on[pi[hit]]=True
    R0=pd.DataFrame({'k':kc,'on':on}); allon=R0.groupby('k').on.transform('all').values; keep=~on|allon
    n_on=int(on.sum()); p_on=pw[on].sum()
    px,py,kc=px[keep],py[keep],kc[keep]; tot=A.set_index('KEY_CODE').JINKO; cnt=pd.Series(kc).value_counts()
    pw=(pd.Series(kc).map(tot)/pd.Series(kc).map(cnt)).values
    pts=shapely.points(px,py); rtree=shapely.STRtree(pts)
    SG=nl.segments(Kz.geometry.values,Kz.id.values); src=Kz.set_index('id')
    E=np.zeros(len(pts))
    for i0 in range(0,len(SG),20000):
        sg=SG.iloc[i0:i0+20000]
        lines=shapely.linestrings(np.stack([np.c_[sg.ax,sg.ay],np.c_[sg.bx,sg.by]],axis=1))
        si,ri=rtree.query(lines,predicate='dwithin',distance=CUT)
        if len(si)==0: continue
        e=nl.seg_energy(sg.ax.values[si],sg.ay.values[si],sg.bx.values[si],sg.by.values[si],px[ri],py[ri],
                        src.S.loc[sg.id.values[si]].values,src.off.loc[sg.id.values[si]].values,H)
        np.add.at(E,ri,e)
    Ln=np.where(E>0,10*np.log10(np.maximum(E,1e-30)),np.nan)
    D=pd.DataFrame({'KEY_CODE':kc,'w':pw,'Ln':Ln})
    D['w45']=D.w*(D.Ln>=45); D['w50']=D.w*(D.Ln>=50); D['w55']=D.w*(D.Ln>=55); D['w60']=D.w*(D.Ln>=60); D['w65']=D.w*(D.Ln>=65); D['wL']=D.w*np.clip(np.nan_to_num(D.Ln,nan=40),40,None)
    out=D.groupby('KEY_CODE')[['w','w45','w50','w55','w60','w65','wL']].sum().reset_index()
    out=out.merge(A[['KEY_CODE','PREF','CITY','JINKO']],on='KEY_CODE')
    out.to_parquet(f'{NAT}/sa_{pc}.parquet')
    print(pc,'smallareas',len(A),'receptors',len(pts),'removed_on_carriageway',n_on,'pop_share_removed',round(p_on/A.JINKO.sum(),4),'sections',len(Kz),'segments',len(SG),' share>=60',round(out.w60.sum()/out.w.sum(),3),round(time.time()-t0),'s',flush=True)
