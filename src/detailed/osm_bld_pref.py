"""Extract OSM building footprints whose representative point falls in a prefecture's bounding box.
usage: osm_bld_pref.py <pbf> <estat shp dir> <pref code> <out parquet>"""
import osmium, sys, shapely.wkb as wkb, pandas as pd, geopandas as gpd
src,shpdir,pc,dst=sys.argv[1:5]
x0,y0,x1,y1=gpd.read_file(f'{shpdir}/r2ka{pc}.shp').to_crs(4326).total_bounds
wf=osmium.geom.WKBFactory(); rows=[]
fp=osmium.FileProcessor(src).with_areas(osmium.filter.KeyFilter('building')).with_filter(osmium.filter.KeyFilter('building'))
for o in fp:
    if not o.is_area(): continue
    try: g=wkb.loads(wf.create_multipolygon(o),hex=True)
    except Exception: continue
    c=g.representative_point()
    if not (x0<c.x<x1 and y0<c.y<y1): continue
    t=o.tags; rows.append((t.get('building'),t.get('building:levels'),t.get('amenity'),t.get('shop'),g))
d=gpd.GeoDataFrame(pd.DataFrame(rows,columns=['building','levels','amenity','shop','geometry']),geometry='geometry',crs=4326)
d.to_parquet(dst); print(pc,'buildings',len(d),'levels non-null',round(d.levels.notna().mean(),3))
