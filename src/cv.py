import sys, pandas as pd, numpy as np
R=pd.read_csv(sys.argv[1],dtype={'pref':str}); R['pref']=R.pref.str.zfill(2)
B={'北海道・東北':range(1,8),'関東':range(8,15),'中部':range(15,24),'近畿':range(24,31),'中国・四国':range(31,40),'九州・沖縄':range(40,48)}
reg={f'{p:02d}':k for k,v in B.items() for p in v}; R['blk']=R.pref.map(reg)
R['pid']=R.pref+'_'+R.point.astype(str)
def met(e,est,meas):
    return dict(n=len(e),bias=e.mean(),sd=e.std(),rmse=np.sqrt((e**2).mean()),mae=e.abs().mean(),w3=(e.abs()<=3).mean()*100,r=np.corrcoef(est,meas)[0,1])
def corr_insample(D):
    c=np.where(D.grp=='gen',D.groupby(['grp','blk']).err.transform('mean'),D.groupby('grp').err.transform('mean'))
    return c
def apply(train,test,mode):
    t=test.copy()
    if mode=='block':
        gb=train[train.grp=='gen'].groupby('blk').err.mean(); gn=train[train.grp=='gen'].err.mean()
        hw=train[train.grp!='gen'].err.mean()
        c=np.where(t.grp=='gen',t.blk.map(gb).fillna(gn),hw)
    elif mode=='national':
        m=train.groupby('grp').err.mean(); c=t.grp.map(m).values
    t['adj']=t.err-c; t['esta']=t.est-c; return t
rows=[]
print(R.grp.value_counts().to_dict(), 'unique points',R.pid.nunique(), 'years',sorted(R.year.unique()))
for g in ['gen','exp']:
    D=R[R.grp==g] if g=='gen' else R[R.grp!='gen']
    rows.append(('補正前',g,met(D.err,D.est,D.meas)))
    c=corr_insample(R)[D.index] if False else None
# in-sample
t=apply(R,R,'block')
for g,D in [('gen',t[t.grp=='gen']),('exp',t[t.grp!='gen'])]: rows.append(('補正後・同一データ（現行）',g,met(D.adj,D.esta,D.meas)))
# LOPO
parts=[apply(R[R.pref!=p],R[R.pref==p],'block') for p in R.pref.unique()]; t=pd.concat(parts)
for g,D in [('gen',t[t.grp=='gen']),('exp',t[t.grp!='gen'])]: rows.append(('都道府県押し出し',g,met(D.adj,D.esta,D.meas)))
lopo=t
# leave-one-block-out (gen: national mean of other blocks)
parts=[apply(R[R.blk!=b],R[R.blk==b],'national') for b in R.blk.unique()]; t=pd.concat(parts)
for g,D in [('gen',t[t.grp=='gen']),('exp',t[t.grp!='gen'])]: rows.append(('地域ブロック押し出し',g,met(D.adj,D.esta,D.meas)))
# temporal
tr=R[R.year<=2021]; te=R[R.year>=2022]; t=apply(tr,te,'block')
for g,D in [('gen',t[t.grp=='gen']),('exp',t[t.grp!='gen'])]: rows.append(('年度分割(2019-21→2022-23)',g,met(D.adj,D.esta,D.meas)))
# point-grouped 10-fold (no point in both)
rng=np.random.default_rng(1); pids=R.pid.unique(); fold=dict(zip(pids,rng.integers(0,10,len(pids)))); R['f']=R.pid.map(fold)
t=pd.concat([apply(R[R.f!=k],R[R.f==k],'block') for k in range(10)])
for g,D in [('gen',t[t.grp=='gen']),('exp',t[t.grp!='gen'])]: rows.append(('地点単位10分割',g,met(D.adj,D.esta,D.meas)))
O=pd.DataFrame([{'scheme':a,'grp':g,**m} for a,g,m in rows]).round(3)
pd.set_option('display.width',250); print(O.to_string())
O.to_csv(sys.argv[2]+'/cv_summary.csv',index=False,encoding='utf-8-sig')
# LOPO by block
G=lopo[lopo.grp=='gen']
pb=G.groupby('blk').apply(lambda x: pd.Series(met(x.adj,x.esta,x.meas))).round(2); print(pb.to_string())
pp=G.groupby('pref').apply(lambda x: pd.Series(met(x.adj,x.esta,x.meas))).round(2)
print('pref LOPO bias range',pp.bias.min(),pp.bias.max(),' |bias|>2 prefs:',(pp.bias.abs()>2).sum(), 'of',len(pp))
pp.to_csv(sys.argv[2]+'/cv_lopo_by_pref.csv',encoding='utf-8-sig')
# Highway by road type
H=R[R.grp!='gen']; print(H.groupby('rt').err.agg(['count','mean','std']).round(2).to_string())
print(R[R.grp=='gen'].groupby('rt').err.agg(['count','mean','std']).round(2).to_string())
# error vs gap, distance, meas level
G=R[R.grp=='gen']; print(G.groupby(pd.cut(G.meas,[0,55,60,65,70,90])).err.agg(['count','mean','std']).round(2).to_string())
