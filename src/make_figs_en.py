# Paper A figures, English labels (for submission). Tables 2 and 3 are written to OUT as CSV with English headers. Inputs: validation points, small-area exposure (data/nat/sa_*.parquet from noise_national.py), small-area boundaries.
import pandas as pd, numpy as np, geopandas as gpd, glob, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
plt.rcParams.update({'ps.fonttype':3,'font.family':'Liberation Sans','font.size':9,'axes.edgecolor':'#52514e','axes.labelcolor':'#0b0b0b',
                     'xtick.color':'#52514e','ytick.color':'#52514e','axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
RAW,OUT=sys.argv[1],sys.argv[2]
import os; os.makedirs(f'{OUT}/eps',exist_ok=True)
BLUE,ORANGE,INK2='#2a78d6','#eb6834','#52514e'
R=pd.read_csv(f'{RAW}/national_validation_points.csv',dtype={'pref':str}); R['pref']=R.pref.str.zfill(2)
BLK={'北海道・東北':range(1,8),'関東':range(8,15),'中部':range(15,24),'近畿':range(24,31),'中国・四国':range(31,40),'九州・沖縄':range(40,48)}
reg={f'{p:02d}':k for k,v in BLK.items() for p in v}; R['blk']=R.pref.map(reg)
R['rtc']=R.rt.str[0].where(R.rt.str[0].isin(['1','2']),'other')
def corr(tr,te):
    g=tr[tr.grp=='gen'].groupby('blk').err.mean(); m=tr[tr.grp=='mw'].groupby('rtc').err.mean()
    return np.where(te.grp=='gen',te.blk.map(g),te.rtc.map(m))
parts=[]
for p in R.pref.unique():
    te=R[R.pref==p].copy(); te['c']=corr(R[R.pref!=p],te); parts.append(te)
T=pd.concat(parts); T['est_cv']=T.est-T.c; T['err_cv']=T.est_cv-T.meas
# ---- Figure 1: predicted vs measured, before / after correction (LOPO) ----
fig,axs=plt.subplots(1,2,figsize=(7.2,3.6),sharex=True,sharey=True)
for ax,col,ttl in [(axs[0],'est','(a) Before correction'),(axs[1],'est_cv','(b) After correction (leave-one-prefecture-out)')]:
    for g,c,lab in [('gen',BLUE,'Ordinary roads'),('mw',ORANGE,'Motorways')]:
        d=T[T.grp==g]; ax.scatter(d.meas,d[col],s=2,alpha=.25 if g=='gen' else .6,color=c,linewidths=0,label=lab,rasterized=True)
    ax.plot([30,95],[30,95],color=INK2,lw=1,ls='--'); ax.set_xlim(35,85); ax.set_ylim(35,95)
    ax.set_title(ttl,loc='left',fontsize=9); ax.set_xlabel('Measured Lnight (dB)'); ax.grid(color='#e5e4e0',lw=.5)
    for g,y0,lab in [('gen',0,'Ordinary roads'),('mw',1,'Motorways')]:
        d=T[T.grp==g]; e=d[col]-d.meas
        ax.text(.03,.97-.08*y0,f'{lab}: mean error {("0.0" if abs(e.mean())<0.05 else f"{e.mean():+.1f}")} dB, RMSE {np.sqrt((e**2).mean()):.1f} dB, r = {np.corrcoef(d[col],d.meas)[0,1]:.2f}',transform=ax.transAxes,va='top',fontsize=7,color='#0b0b0b')
axs[0].set_ylabel('Estimated Lnight (dB)')
h=[plt.Line2D([],[],marker='o',ls='',color=c,ms=4) for c in (BLUE,ORANGE)]; axs[1].legend(h,['Ordinary roads (n = 11,978)','Motorways (n = 670)'],loc='lower right',frameon=False,fontsize=7)
fig.tight_layout(); fig.savefig(f'{OUT}/fig1_validation_scatter.png',dpi=300); fig.savefig(f'{OUT}/eps/fig1_validation_scatter.eps'); plt.close()
# ---- Figure 3: residual mean error by prefecture (LOPO, ordinary roads) ----
T['site']=T.pref+'_'+T.jis.astype(str)+'_'+T.point.astype(str)   # monitoring site = prefecture + municipality + point code
def cr_se(d):   # cluster-robust (CR1) standard error of the mean error, clustered by monitoring site (several fiscal years per site)
    u=d.err_cv-d.err_cv.mean(); s=u.groupby(d.site).sum(); G=len(s); return np.sqrt((s**2).sum()*G/(G-1))/len(d)
P=pd.DataFrame([{'pref':p,'n':len(d),'sites':d.site.nunique(),'bias':d.err_cv.mean(),'se':cr_se(d)} for p,d in T[T.grp=='gen'].groupby('pref')])
P.to_csv(f'{OUT}/fig3_data.csv',index=False); PN=pd.read_csv(f'{RAW}/pref_names.csv',dtype={'pref':str}) if False else None
names=dict(zip([f'{i:02d}' for i in range(1,48)],'Hokkaido Aomori Iwate Miyagi Akita Yamagata Fukushima Ibaraki Tochigi Gunma Saitama Chiba Tokyo Kanagawa Niigata Toyama Ishikawa Fukui Yamanashi Nagano Gifu Shizuoka Aichi Mie Shiga Kyoto Osaka Hyogo Nara Wakayama Tottori Shimane Okayama Hiroshima Yamaguchi Tokushima Kagawa Ehime Kochi Fukuoka Saga Nagasaki Kumamoto Oita Miyazaki Kagoshima Okinawa'.split()))
P=P.sort_values('bias'); P['nm']=P.pref.map(names)
fig,ax=plt.subplots(figsize=(3.6,7.2)); y=np.arange(len(P))
ax.axvspan(-2,2,color='#f0efec',zorder=0); ax.axvline(0,color=INK2,lw=.8)
ax.errorbar(P.bias,y,xerr=1.96*P.se,fmt='o',ms=4,color=BLUE,ecolor='#86b6ef',elinewidth=1.2,capsize=0)
ax.set_yticks(y); ax.set_yticklabels(P.nm,fontsize=7); ax.set_xlabel('Remaining mean error after correction (dB)\nordinary roads, leave-one-prefecture-out,\n95% CI clustered by monitoring site')
ax.set_ylim(-1,len(P)); ax.text(2.05,-0.6,'±2 dB',ha='left',fontsize=7,color=INK2)
fig.tight_layout(); fig.savefig(f'{OUT}/fig3_residual_bias_by_pref.png',dpi=300); fig.savefig(f'{OUT}/eps/fig3_residual_bias_by_pref.eps'); plt.close()
# ---- Table 2: accuracy by scheme ----
def met(e,p,m): return pd.Series({'n':len(e),'mean_error_dB':e.mean(),'sd_dB':e.std(),'rmse_dB':np.sqrt((e**2).mean()),'mae_dB':e.abs().mean(),'within_3dB_pct':(e.abs()<=3).mean()*100,'r':np.corrcoef(p,m)[0,1],'rho':pd.Series(np.asarray(p)).rank().corr(pd.Series(np.asarray(m)).rank())})
rows=[]
for g,lab in [('gen','Ordinary roads'),('mw','Motorways')]:
    d=R[R.grp==g]; rows.append(('Before correction',lab,met(d.err,d.est,d.meas)))
    ins=R.copy(); ins['c']=corr(R,R); d=ins[ins.grp==g]; rows.append(('After correction, same data',lab,met(d.err-d.c,d.est-d.c,d.meas)))
    d=T[T.grp==g]; rows.append(('Leave-one-prefecture-out',lab,met(d.err_cv,d.est_cv,d.meas)))
    BL=[]; 
    for b in R.blk.unique():
        tr,te=R[R.blk!=b],R[R.blk==b].copy()
        gm=tr[tr.grp=='gen'].err.mean(); mm=tr[tr.grp=='mw'].groupby('rtc').err.mean(); te['c']=np.where(te.grp=='gen',gm,te.rtc.map(mm)); BL.append(te)
    d=pd.concat(BL); d=d[d.grp==g]; rows.append(('Leave-one-block-out',lab,met(d.err-d.c,d.est-d.c,d.meas)))
    tr,te=R[R.year<=2021],R[R.year>=2022].copy(); te['c']=corr(tr,te); d=te[te.grp==g]; rows.append(('Split by year (2019-21 -> 2022-23)',lab,met(d.err-d.c,d.est-d.c,d.meas)))
    rng=np.random.default_rng(1); R['pid']=R.pref+'_'+R.jis.astype(str)+'_'+R.point.astype(str); pids=R.pid.unique(); fo=dict(zip(pids,rng.integers(0,10,len(pids)))); R['f']=R.pid.map(fo)
    FK=[]
    for k in range(10):
        te=R[R.f==k].copy(); te['c']=corr(R[R.f!=k],te); FK.append(te)
    d=pd.concat(FK); d=d[d.grp==g]; rows.append(('10-fold by site',lab,met(d.err-d.c,d.est-d.c,d.meas)))
T2=pd.DataFrame([{'scheme':a,'roads':b,**m} for a,b,m in rows]); T2.round(2).to_csv(f'{OUT}/table2_accuracy.csv',index=False); print(T2.round(2).to_string())
# ---- Table 3 & Figure 2: exposure ----
E=pd.concat([pd.read_parquet(f) for f in glob.glob(f'{RAW}/nat/sa_*.parquet')]); E['pref']=E.PREF.astype(str).str.zfill(2)
E['muni']=E.pref+E.CITY.astype(str).str.zfill(3)
P3=E.groupby('pref')[['w','w55','w60','w65','wL']].sum(); tot=P3.sum()
T3=pd.DataFrame({'prefecture':P3.index.map(names),'population_10k':P3.w/1e4,'share_ge55_pct':P3.w55/P3.w*100,'share_ge60_pct':P3.w60/P3.w*100,'share_ge65_pct':P3.w65/P3.w*100,'mean_Lnight_dB':P3.wL/P3.w})
T3.loc['00']=['Total',tot.w/1e4,tot.w55/tot.w*100,tot.w60/tot.w*100,tot.w65/tot.w*100,tot.wL/tot.w]
T3.round(1).to_csv(f'{OUT}/table3_exposure_by_pref.csv'); print(T3.round(1).tail(3).to_string())
M=E.groupby('muni')[['w','w60']].sum(); M['p60']=M.w60/M.w*100
G=[]
for f in sorted(glob.glob(f'{RAW}/estat_all/shp/r2ka*.shp')):
    A=gpd.read_file(f,columns=['PREF','CITY','HCODE','geometry']); A=A[A.HCODE==8101].copy(); A['geometry']=A.geometry.make_valid()
    A['muni']=A.PREF.astype(str).str.zfill(2)+A.CITY.astype(str).str.zfill(3); D_=A.dissolve('muni').reset_index()[['muni','geometry']]; D_['geometry']=D_.geometry.simplify(0.001); G.append(D_)
G=pd.concat(G); G=gpd.GeoDataFrame(G,crs=G.crs).to_crs(6677); G=G.merge(M[['p60']],left_on='muni',right_index=True,how='left')

bins=[0,1,3,5,10,20,100]; cols=['#cde2fb','#9ec5f4','#6da7ec','#3987e5','#256abf','#104281']
G['cls']=pd.cut(G.p60.fillna(0),bins,include_lowest=True,labels=False)
fig,ax=plt.subplots(figsize=(6.5,7))
isl=G.muni.str.startswith('47')|G.muni.isin([f'13{c}' for c in ['361','362','363','364','381','382','401','402','421']])
def draw(a,D):
    for i,c in enumerate(cols): D[D.cls==i].plot(ax=a,color=c,lw=0)
    a.set_axis_off(); a.set_aspect('equal')
draw(ax,G[~isl])
xmin,ymin,xmax,ymax=G[~isl].total_bounds; ax.set_xlim(xmin,xmax); ax.set_ylim(ymin,ymax)
ia=ax.inset_axes([0.04,0.66,0.2,0.26]); O=G[G.muni.str.startswith('47')]; draw(ia,O)
import shapely as _sh
bb=gpd.GeoSeries([_sh.box(127.6,26.05,128.35,26.9)],crs=4326).to_crs(6677).total_bounds
ia.set_xlim(bb[0],bb[2]); ia.set_ylim(bb[1],bb[3]); ia.set_axis_on(); ia.set_xticks([]); ia.set_yticks([]); ia.set_xlabel(''); ia.set_ylabel('')
for sp in ia.spines.values(): sp.set_color('#c3c2b7'); sp.set_visible(True)
ia.set_title('Okinawa Island',fontsize=7,loc='left')
h=[plt.Rectangle((0,0),1,1,color=c) for c in cols]; ax.legend(h,['<1%','1–3%','3–5%','5–10%','10–20%','≥20%'],title='Share of residents with\nLnight ≥ 60 dB (municipality)',loc='lower right',frameon=False,fontsize=7,title_fontsize=8)
fig.tight_layout(); fig.savefig(f'{OUT}/fig2_map_municipal_p60.png',dpi=300); fig.savefig(f'{OUT}/eps/fig2_map_municipal_p60.eps'); plt.close()
print('done')
