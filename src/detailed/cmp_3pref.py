"""National model (v4) vs detailed model with buildings in Osaka, Hyogo and Wakayama: prefecture summary (all small areas and small
areas with complete building data) and municipal values. Also the detailed receptors without screening (building-based placement).
usage: cmp_3pref_v4.py <national nat dir with w50> <detailed receptors dir> <out dir>"""
import pandas as pd, numpy as np, sys, os
NAT,REC,OUT=sys.argv[1:4]; os.makedirs(OUT,exist_ok=True)
S=[];M=[]
for pc in ['27','28','30']:
    n=pd.read_parquet(f'{NAT}/sa_{pc}.parquet'); n['muni']=pc+n.CITY.astype(str).str.zfill(3)
    r=pd.read_parquet(f'{REC}/receptors_{pc}.parquet')
    r['L0']=np.where(r.E_open+r.E_blk>0,10*np.log10(np.maximum(r.E_open+r.E_blk,1e-30)),np.nan)
    grid_sa=set(r.KEY_CODE[r.kind=='grid'])
    for st,keep_n,keep_r in [('all',np.ones(len(n),bool),np.ones(len(r),bool)),('bldg_only',~n.KEY_CODE.isin(grid_sa).values,~r.KEY_CODE.isin(grid_sa).values)]:
        a=n[keep_n]; b=r[keep_r]; W=a.w.sum(); P=b['pop'].sum(); L=b.Ln_A10
        row={'pref':pc,'set':st,'pop_M':round(W/1e6,2)}
        for t in (50,55,60,65): row[f'simple_p{t}']=a[f'w{t}'].sum()/W*100; row[f'detail_p{t}']=b['pop'][L>=t].sum()/P*100
        row['simple_L']=a.wL.sum()/W; row['detail_L']=np.average(L.fillna(40).clip(lower=40),weights=b['pop'])
        for t in (55,60): row[f'unscreened_p{t}']=b['pop'][b.L0>=t].sum()/P*100
        row['unscreened_L']=np.average(b.L0.fillna(40).clip(lower=40),weights=b['pop'])
        S.append(row)
    g=n.groupby('muni')[['w','w50','w55','w60','w65','wL']].sum()
    q=r.groupby('muni_code').apply(lambda x: pd.Series({'dpop':x['pop'].sum(),**{f'q{t}':x['pop'][x.Ln_A10>=t].sum()/x['pop'].sum() for t in (50,55,60,65)},
        'qL':np.average(x.Ln_A10.fillna(40).clip(lower=40),weights=x['pop'])}))
    m=pd.DataFrame({**{f's{t}':g[f'w{t}']/g.w for t in (50,55,60,65)},'sL':g.wL/g.w,'pop':g.w}).join(q,how='inner'); m['pref']=pc; M.append(m)
S=pd.DataFrame(S).round(1); S.to_csv(f'{OUT}/summary_3pref_v4.csv',index=False); print(S.to_string())
M=pd.concat(M); M.index.name='muni'; M.reset_index().to_csv(f'{OUT}/municipal_3pref_v4.csv',index=False)
for k,x in [('all',M)]+[(p,M[M.pref==p]) for p in ['27','28','30']]:
    print(k,'n',len(x),' r p55 %.2f p60 %.2f Lmean %.2f  spearman p60 %.2f'%(np.corrcoef(x.s55,x.q55)[0,1],np.corrcoef(x.s60,x.q60)[0,1],np.corrcoef(x.sL,x.qL)[0,1],x.s60.rank().corr(x.q60.rank())))
