import pandas as pd, numpy as np, geopandas as gpd, shapely, glob, sys, unicodedata as ud
sys.path.insert(0,sys.argv[2]); import noise_lib as nl
S=sys.argv[1]
prof=nl.night_profile(f'{S}/census27.csv')
K=nl.load_sections(f'{S}/geo/traffic_census_2021_converted.parquet')
K=K.drop_duplicates('id')
SRC=nl.source_table(K,prof)
K=K.merge(SRC[['id','S','off','ok','grp']],on='id'); K=K[K.ok]
norm=lambda s: ud.normalize('NFKC',str(s)).replace(' ','').replace('　','')
K['rn']=K.route.map(norm); K['rn0']=K.rn.str.replace(r'[（(].*?[）)]','',regex=True)
def rd(f):
    for e in ('utf-8-sig','cp932'):
        try: return pd.read_csv(f,encoding=e,dtype=str)
        except UnicodeDecodeError: pass
M=pd.concat([rd(f) for f in sorted(glob.glob(f'{S}/tn/**/souon*.csv',recursive=True))]); M.columns=[c.strip() for c in M.columns]
for c in ['longitude','latitude','distance_from_roadway','measure_height','night_leq']: M[c]=pd.to_numeric(M[c],errors='coerce')
M=M.dropna(subset=['longitude','latitude','night_leq']); M['rn']=M.road1_name.map(norm); M['rn0']=M.rn.str.replace(r'[（(].*?[）)]','',regex=True)
M['pref']=M.pref_code.str.zfill(2)
out=[]
for z,mz in M.groupby(M.longitude.map(nl.utm_epsg)):
    Kz=K.to_crs(z); tr=shapely.STRtree(Kz.geometry.values)
    pts=gpd.GeoSeries(gpd.points_from_xy(mz.longitude,mz.latitude),crs=4326).to_crs(z).values
    pi,ki=tr.query(pts,predicate='dwithin',distance=150)
    df=pd.DataFrame({'pi':pi,'ki':ki}); df['dist']=shapely.distance(pts[pi],Kz.geometry.values[ki])
    df['same']=(Kz.rn.values[ki]==mz.rn.values[pi])|(Kz.rn0.values[ki]==mz.rn0.values[pi])
    df=df[df.same].sort_values('dist').drop_duplicates('pi')
    for _,r in df.iterrows():
        m=mz.iloc[int(r.pi)]; k=Kz.iloc[int(r.ki)]
        # receptor at measurement point: use segments of the matched section only (as in Osaka validation: edge distance from record)
        geom=k.geometry; segs=nl.segments([geom],[k.id])
        # place receptor on perpendicular at (centerline -> road edge + distance_from_roadway); approximate with actual point distance
        p=pts[int(r.pi)]; px,py=shapely.get_x(p),shapely.get_y(p)
        # use recorded distance from carriageway edge to define distance from centerline
        dcl=(k.carw/2 if k.oneway!='1' else k.carw/2)+ (m.distance_from_roadway if m.distance_from_roadway==m.distance_from_roadway else 0)
        # move receptor along the normal so that its distance to centerline equals dcl
        near=shapely.ops.nearest_points(geom,p)[0] if hasattr(shapely,'ops') else None
        import shapely.ops as so
        q=so.nearest_points(geom,p)[0]; vx,vy=px-q.x,py-q.y; vl=np.hypot(vx,vy) or 1.0
        rx,ry=q.x+vx/vl*dcl,q.y+vy/vl*dcl
        E=nl.seg_energy(segs.ax.values,segs.ay.values,segs.bx.values,segs.by.values,np.full(len(segs),rx),np.full(len(segs),ry),
                        np.full(len(segs),k.S),np.full(len(segs),k.off),m.measure_height if m.measure_height==m.measure_height else 1.2).sum()
        out.append({'pref':m.pref,'year':m.report_year,'point':m.noise_measure_point_code,'jis':m.jis_code,'grp':k.grp,'rt':k.rt,'gap':r.dist,'meas':m.night_leq,'est':10*np.log10(E) if E>0 else np.nan})
R=pd.DataFrame(out).dropna(subset=['est']); R['err']=R.est-R.meas
R.to_csv(f'{S}/national_validation_points.csv',index=False)
f=lambda x: pd.Series({'n':len(x),'meas_med':x.meas.median(),'bias':x.err.mean(),'sd':x.err.std(),'r':np.corrcoef(x.est,x.meas)[0,1],'within3':(x.err.abs()<=3).mean()})
pd.set_option('display.width',200)
print('matched points',len(R),' of measurement records',len(M))
print(R.groupby('grp').apply(f).round(2).to_string())
G=R[R.grp=='gen']
reg={**{p:'北海道・東北' for p in ['01','02','03','04','05','06','07']},**{p:'関東' for p in ['08','09','10','11','12','13','14']},**{p:'中部' for p in ['15','16','17','18','19','20','21','22','23']},
     **{p:'近畿' for p in ['24','25','26','27','28','29','30']},**{p:'中国・四国' for p in ['31','32','33','34','35','36','37','38','39']},**{p:'九州・沖縄' for p in ['40','41','42','43','44','45','46','47']}}
print(G.groupby(G.pref.map(reg)).apply(f).round(2).to_string())
print('Osaka only (check vs earlier +2.07):', round(G[G.pref=='27'].err.mean(),2), len(G[G.pref=='27']))
