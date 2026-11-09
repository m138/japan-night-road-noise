import sys, pandas as pd, numpy as np
R=pd.read_csv(sys.argv[1],dtype={'pref':str}); R['pref']=R.pref.str.zfill(2)
B={'北海道・東北':range(1,8),'関東':range(8,15),'中部':range(15,24),'近畿':range(24,31),'中国・四国':range(31,40),'九州・沖縄':range(40,48)}
reg={f'{p:02d}':k for k,v in B.items() for p in v}; R['blk']=R.pref.map(reg)
R['rtc']=R.rt.str[0]
def met(e,est,meas): return dict(n=len(e),bias=round(e.mean(),2),sd=round(e.std(),2),rmse=round(np.sqrt((e**2).mean()),2),mae=round(e.abs().mean(),2),w3=round((e.abs()<=3).mean()*100,1),r=round(np.corrcoef(est,meas)[0,1],3))
def fit(tr,te,mode):
    out=[]
    for g in ['gen','mw']:
        a=tr[tr.grp==g]; b=te[te.grp==g].copy()
        if len(b)==0: continue
        if mode=='block': c=b.blk.map(a.groupby('blk').err.mean()).fillna(a.err.mean()) if g=='gen' else a.err.mean()
        elif mode=='rt': c=b.rtc.map(a.groupby('rtc').err.mean()).fillna(a.err.mean())
        elif mode=='rt+block':
            m=a.groupby(['rtc','blk']).err.mean(); k=list(zip(b.rtc,b.blk))
            c=pd.Series([m.get(x,np.nan) for x in k],index=b.index).fillna(b.rtc.map(a.groupby('rtc').err.mean())).fillna(a.err.mean())
        elif mode=='linear':  # meas = a + b*est (+rt dummies)
            X=lambda d: np.column_stack([np.ones(len(d)),d.est]+[(d.rtc==r).astype(float) for r in sorted(a.rtc.unique())[1:]])
            beta=np.linalg.lstsq(X(a),a.meas,rcond=None)[0]; pred=X(b)@beta; c=b.est-pred
        b['p']=b.est-c; b['e']=b.p-b.meas; out.append(b)
    return pd.concat(out)
rows=[]
for mode in ['block','rt','rt+block','linear']:
    t=pd.concat([fit(R[R.pref!=p],R[R.pref==p],mode) for p in R.pref.unique()])
    for g in ['gen','mw']:
        D=t[t.grp==g]; rows.append({'補正方法':mode,'grp':g,**met(D.e,D.p,D.meas)})
O=pd.DataFrame(rows); print(O.to_string()); O.to_csv(''+(sys.argv[2] if len(sys.argv)>2 else '.')+'/cv_methods_lopo.csv',index=False,encoding='utf-8-sig')
G=R[R.grp=='gen']; print(G.groupby(pd.qcut(G.est,5)).apply(lambda x: pd.Series({'n':len(x),'est':x.est.mean(),'meas':x.meas.mean(),'err':x.err.mean()})).round(2).to_string())
print('slope meas~est gen', np.polyfit(G.est,G.meas,1).round(3))
