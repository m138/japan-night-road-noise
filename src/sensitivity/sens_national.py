"""Sensitivity runs for the national model (receptor placement, 40 dB floor, 200 m cutoff, ASJ RTN-Model 2023).
Identical to noise_national.py (v3) except that per-receptor results are kept and the following are added:
 - E200, E500: energy from census sections within 200 m / 500 m (main model = E200)
 - E200_23: E200 with ASJ RTN-Model 2023 non-steady light-vehicle coefficient (81.4) and its own bias correction
 - dmin: distance to the nearest section centerline within 500 m; oncw: receptor lies on a carriageway (dmin < carriageway width/2)
usage: sens_national.py <v3 dir> <noise_lib dir> <validation csv ASJ2023> <seed> <cut: 200|500> <out dir> [prefs]"""
import pandas as pd, numpy as np, geopandas as gpd, shapely, sys, glob, os, time
V3,LIB,VAL23,SEED,CUTMAX,OUT=sys.argv[1:7]; PREFS=sys.argv[7].split(',') if len(sys.argv)>7 else None
sys.path.insert(0,LIB); import noise_lib as nl
H=4.0; PER=10; CUTMAX=float(CUTMAX); rng=np.random.default_rng(int(SEED)); os.makedirs(OUT,exist_ok=True)
REG={**{p:'北海道・東北' for p in ['01','02','03','04','05','06','07']},**{p:'関東' for p in ['08','09','10','11','12','13','14']},**{p:'中部' for p in ['15','16','17','18','19','20','21','22','23']},
     **{p:'近畿' for p in ['24','25','26','27','28','29','30']},**{p:'中国・四国' for p in ['31','32','33','34','35','36','37','38','39']},**{p:'九州・沖縄' for p in ['40','41','42','43','44','45','46','47']}}
def bias(f):
    V=pd.read_csv(f,dtype={'pref':str}); V['pref']=V.pref.str.zfill(2)
    bg=V[V.grp=='gen'].groupby(V.pref.map(REG)).err.mean().to_dict(); V['rtc']=V.rt.astype(str).str[0].where(V.rt.astype(str).str[0].isin(['1','2']),'other')
    return bg,V[V.grp=='mw'].groupby('rtc').err.mean().to_dict()
bg18,bm18=bias(f'{V3}/national_validation_points.csv'); bg23,bm23=bias(VAL23)
print('bias 2018',bg18,bm18,'\nbias 2023',bg23,bm23,flush=True)
prof=nl.night_profile(f'{V3}/census27.csv')
K=nl.load_sections(f'{V3}/geo/traffic_census_2021_converted.parquet')
SRC=nl.source_table(K.drop_duplicates('id'),prof).set_index('id')
# ASJ 2023: only the non-steady light-vehicle coefficient differs (82.3 -> 81.4); recompute S with that coefficient
KK=K.drop_duplicates('id'); st=(KK.motorway=='1')|KK.rt.isin(['1','2'])
Vk=KK.vlim.where(KK.vlim>0); Vk=np.where(st,np.clip(Vk.fillna(80),40,120),np.clip(Vk.fillna(40),10,60)); S23=np.zeros(len(KK))
for p in ('u','d'):
    for veh,s in (('1','s'),('2','l')):
        n12=(KK[f'{p}24{s}']-KK[f'{p}12{s}']).clip(lower=0).fillna(0).values
        pf=np.where(st,prof[('mw',veh)].mean(),prof[('gen',veh)].mean())
        L=np.where(st,(45.8 if veh=='1' else 53.2)+30*np.log10(Vk),(81.4 if veh=='1' else 88.8)+10*np.log10(Vk))
        S23+=10**((L-8)/10)*n12*pf/(1000*Vk)
SRC['S23']=pd.Series(S23,index=KK.id.values).reindex(SRC.index).values
K=K[K.id.isin(SRC.index[SRC.ok])].dissolve(by='id').reset_index(); K['pref']=K.muni.str[:2]
rtc=K.rt.astype(str).str[0].where(K.rt.astype(str).str[0].isin(['1','2']),'other'); mw=SRC.loc[K.id,'grp'].values=='mw'
b18=np.where(mw,rtc.map(bm18).values,K.pref.map(REG).map(bg18).values); b23=np.where(mw,rtc.map(bm23).values,K.pref.map(REG).map(bg23).values)
K['S']=SRC.loc[K.id,'S'].values*10**(-b18/10); K['S23']=SRC.loc[K.id,'S23'].values*10**(-b23/10)
K['off']=SRC.loc[K.id,'off'].values; K['hw']=K.carw.values/2
K['geometry']=K.geometry.simplify(0.00002)
for f in sorted(glob.glob(f'{V3}/estat_all/shp/r2ka*.shp')):
    pc=os.path.basename(f)[4:6]
    if PREFS and pc not in PREFS: continue
    t0=time.time()
    A=gpd.read_file(f); A=A[A.HCODE==8101][['KEY_CODE','PREF','CITY','JINKO','geometry']]
    A=A.dissolve(by='KEY_CODE',aggfunc={'PREF':'first','CITY':'first','JINKO':'sum'}).reset_index(); A=A[A.JINKO>0]
    z=nl.utm_epsg(A.to_crs(4326).union_all().centroid.x); A=A.to_crs(z)
    n=np.maximum(1,np.round(A.JINKO.values/PER).astype(int)); px=[];py=[];kc=[];pw=[]
    for (key,pop,geom),k in zip(A[['KEY_CODE','JINKO','geometry']].itertuples(index=False),n):
        x0,y0,x1,y1=geom.bounds; got=np.empty((0,2)); tries=0
        while len(got)<k and tries<50:
            m=max(4*(k-len(got)),16); c=np.c_[rng.uniform(x0,x1,m),rng.uniform(y0,y1,m)]
            ins=shapely.contains_xy(geom,c[:,0],c[:,1]); got=np.vstack([got,c[ins]]); tries+=1
        if len(got)==0: rp=geom.representative_point(); got=np.array([[rp.x,rp.y]])
        got=got[:k]; px.append(got[:,0]); py.append(got[:,1]); kc+=[key]*len(got); pw.append(np.full(len(got),pop/len(got)))
    px=np.concatenate(px); py=np.concatenate(py); pw=np.concatenate(pw); kc=np.array(kc)
    pts=shapely.points(px,py); rtree=shapely.STRtree(pts)
    Kz=K.to_crs(z); Kz=Kz[Kz.intersects(shapely.box(*A.total_bounds).buffer(CUTMAX+100))]
    SG=nl.segments(Kz.geometry.values,Kz.id.values); src=Kz.set_index('id')
    E200=np.zeros(len(pts)); E500=np.zeros(len(pts)); E23=np.zeros(len(pts)); dmin=np.full(len(pts),np.inf); oncw=np.zeros(len(pts),bool)
    for i0 in range(0,len(SG),20000):
        sg=SG.iloc[i0:i0+20000]
        lines=shapely.linestrings(np.stack([np.c_[sg.ax,sg.ay],np.c_[sg.bx,sg.by]],axis=1))
        si,ri=rtree.query(lines,predicate='dwithin',distance=CUTMAX)
        if len(si)==0: continue
        ids=sg.id.values[si]; dist=shapely.distance(lines[si],pts[ri])
        a=(sg.ax.values[si],sg.ay.values[si],sg.bx.values[si],sg.by.values[si],px[ri],py[ri])
        off=src.off.loc[ids].values
        e=nl.seg_energy(*a,src.S.loc[ids].values,off,H); in2=dist<=200
        np.add.at(E500,ri,e); np.add.at(E200,ri[in2],e[in2])
        e23=nl.seg_energy(*a,src.S23.loc[ids].values,off,H); np.add.at(E23,ri[in2],e23[in2])
        np.minimum.at(dmin,ri,dist); hw=src.hw.loc[ids].values; on=dist<hw; oncw[ri[on]]=True
    D=pd.DataFrame({'KEY_CODE':kc,'w':pw,'E200':E200,'E500':E500,'E200_23':E23,'dmin':dmin,'oncw':oncw})
    D=D.merge(A[['KEY_CODE','PREF','CITY']],on='KEY_CODE',how='left'); D.to_parquet(f'{OUT}/rec_{pc}.parquet')
    print(pc,'receptors',len(D),round(time.time()-t0),'s',flush=True)
