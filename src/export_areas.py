"""Export the exposure indicators by municipality and by prefecture as CSV and as a GeoPackage of municipal polygons.
Municipalities follow the 2020 census codes (wards of designated cities separately). Names are in Japanese, as in e-Stat.
usage: export_areas.py <data dir> <out dir>
  -> municipalities.csv, prefectures.csv, municipalities.gpkg (layer 'municipalities', EPSG:6668 (JGD2011), simplified to about 10 m)"""
import pandas as pd, geopandas as gpd, glob, sys, os
D,OUT=sys.argv[1:3]; os.makedirs(OUT,exist_ok=True)
E=pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{D}/nat/sa_*.parquet'))])
E['PREF']=E.PREF.astype(str).str.zfill(2); E['CITY']=E.CITY.astype(str).str.zfill(3)
TH=(45,50,55,60,65)
def summ(g,keys):
    s=g.groupby(keys)[['w',*[f'w{t}' for t in TH],'wL']].sum()
    o=pd.DataFrame({'population':s.w.round(0)})
    for t in TH: o[f'share_ge{t}_pct']=(s[f'w{t}']/s.w*100).round(3)
    o['mean_Lnight_dB']=(s.wL/s.w).round(3); return o.reset_index()
A=[]
for f in sorted(glob.glob(f'{D}/estat_all/shp/r2ka*.shp')):
    if f[-6:-4] not in set(E.PREF): continue
    a=gpd.read_file(f,columns=['PREF','CITY','PREF_NAME','CITY_NAME','HCODE','geometry']); a=a[a.HCODE==8101].copy()
    a['PREF']=a.PREF.astype(str).str.zfill(2); a['CITY']=a.CITY.astype(str).str.zfill(3); a['geometry']=a.geometry.make_valid(); A.append(a)
A=gpd.GeoDataFrame(pd.concat(A),crs=A[0].crs)
names=A.drop_duplicates(['PREF','CITY'])[['PREF','CITY','PREF_NAME','CITY_NAME']]
M=summ(E,['PREF','CITY']).merge(names,on=['PREF','CITY'],how='left')
M.insert(0,'municipality_code',M.PREF+M.CITY)
M=M.rename(columns={'PREF_NAME':'prefecture_ja','CITY_NAME':'municipality_ja'})[['municipality_code','prefecture_ja','municipality_ja','population',*[f'share_ge{t}_pct' for t in TH],'mean_Lnight_dB']]
M.to_csv(f'{OUT}/municipalities.csv',index=False)
P=summ(E,['PREF']).merge(names.drop_duplicates('PREF')[['PREF','PREF_NAME']],on='PREF')
P=P.rename(columns={'PREF':'prefecture_code','PREF_NAME':'prefecture_ja'})[['prefecture_code','prefecture_ja','population',*[f'share_ge{t}_pct' for t in TH],'mean_Lnight_dB']]
P.to_csv(f'{OUT}/prefectures.csv',index=False)
A['municipality_code']=A.PREF+A.CITY
G=A[['municipality_code','geometry']].dissolve('municipality_code').reset_index().to_crs(6668)
G['geometry']=G.geometry.simplify(0.0001)
G=G.merge(M,on='municipality_code',how='inner'); G.to_file(f'{OUT}/municipalities.gpkg',layer='municipalities',driver='GPKG')
print('municipalities',len(M),'prefectures',len(P),'population',int(M.population.sum()))
