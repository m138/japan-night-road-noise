"""Three prefectures: Spearman and partial Spearman (rank residuals on rank of log population-weighted density, Pearson of residuals)
between national (s*) and detailed (q*) municipal values. Reproduces density_check_3pref.csv. usage: dens_check_3pref.py <estat shp dir> <municipal_3pref_v4.csv> <out csv>"""
import pandas as pd,numpy as np,sys,geopandas as gpd
SHP,MUN,OUT=sys.argv[1:4]
G=pd.concat([gpd.read_file(f,ignore_geometry=True)[['PREF','CITY','AREA','JINKO']] for f in [f'{SHP}/r2ka{p}.shp' for p in ('27','28','30')]])
G['muni']=G.PREF.astype(str).str.zfill(2)+G.CITY.astype(str).str.zfill(3); G=G[G.AREA>0]; G['d']=G.JINKO/(G.AREA/1e6)
pwd=(G.JINKO*G.d).groupby(G.muni).sum()/G.groupby('muni').JINKO.sum()
d=pd.read_csv(MUN,dtype={'muni':str}).set_index('muni'); d['pwd']=pwd.reindex(d.index); d=d.dropna(subset=['pwd'])
def resid(y,x):
    yr,xr=y.rank(),np.log(x).rank(); b=np.polyfit(xr,yr,1); return yr-np.polyval(b,xr)
rows=[{'indicator':a[1:],'rho':d[a].rank().corr(d[b].rank()),'partial_rho_density':resid(d[a],d.pwd).corr(resid(d[b],d.pwd)),
       'rho_national_vs_density':d[a].rank().corr(d.pwd.rank()),'rho_detailed_vs_density':d[b].rank().corr(d.pwd.rank())} for a,b in [('s60','q60'),('s55','q55'),('sL','qL')]]
pd.DataFrame(rows).to_csv(OUT,index=False)
