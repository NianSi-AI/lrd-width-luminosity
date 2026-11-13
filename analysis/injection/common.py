"""Task10 deterministic M3; log10 luminosity and width, covariance in dex^2.

M3 likelihood/priors and 801x601 midpoint cells are unchanged. Only posterior
summarization replaces random draws with marginal CDFs and direct convolution.
All writers refuse existing paths; imports do not run legacy writer scripts.
"""
from pathlib import Path
import json, csv, sys, hashlib
import numpy as np
from scipy.special import ndtr, logsumexp
from scipy.optimize import brentq

O=Path('work/injection').resolve()
B=Path('.')
ROOT=B
sys.path.insert(0,str(B))
Q=np.array([.025,.16,.5,.84,.975])
def native(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {str(k):native(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [native(x) for x in v]
    return v
def read(p):return json.loads(Path(p).read_text())
def save(p,x):
    p=Path(p);p=p if p.is_absolute() else O/p
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(native(x),f,ensure_ascii=False,indent=2,allow_nan=False)
def tsv(p,rows):
    p=Path(p);p=p if p.is_absolute() else O/p
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with p.open('x',newline='') as f:
        w=csv.DictWriter(f,keys,delimiter='\t');w.writeheader()
        w.writerows({k:json.dumps(native(v),ensure_ascii=False) if isinstance(v,(dict,list,np.ndarray)) else v for k,v in r.items()} for r in rows)
def quant(a,p,qs=Q):
    """Cell masses -> CDF at boundaries, including the lower support edge."""
    a=np.asarray(a); p=np.asarray(p,dtype=float); p=p/p.sum()
    d=np.diff(a); assert np.allclose(d,d[0]); step=d[0]
    edges=np.r_[a[0]-step/2,a+step/2]
    cdf=np.r_[0,np.cumsum(p)]; cdf[-1]=1
    return np.interp(qs,cdf,edges)

cache={}
def difference(g1,g2):
    s1,p1=g1['s'],g1['ps'];s2,p2=g2['s'],g2['ps']
    assert len(s1)==len(s2) and np.array_equal(s1,s2)
    p=np.convolve(p1,p2[::-1]);axis=np.arange(-(len(s1)-1),len(s1))*(s1[1]-s1[0])
    return dict(quantiles=Q,LRD_minus_nonLRD_dex=np.interp(Q,np.cumsum(p)/p.sum(),axis),P_lt0=float(p[axis<0].sum()/p.sum())),(axis,p)
def posterior(rows,nb=801,ns=601,intercept=False):
    assert len(rows)>=4
    x=np.array([r['x'] for r in rows]);y=np.array([r['y'] for r in rows]);C=np.array([r['cov'] for r in rows])
    assert np.isfinite(C).all() and np.isfinite(x+y).all()
    assert np.linalg.eigvalsh(C).min()>-1e-10
    b=-1+(np.arange(nb)+.5)*2/nb;s=(np.arange(ns)+.5)*.6/ns
    bb=b[:,None];ss=s[None,:];shape=(nb,ns)
    sw=np.zeros(shape);sr=sw.copy();sr2=sw.copy();slv=sw.copy()
    for xi,yi,ci in zip(x,y,C):
        V=ss**2+ci[1,1]+bb**2*ci[0,0]-2*bb*ci[0,1]
        assert V.min()>0
        r=yi-bb*xi;w=1/V;sw+=w;sr+=w*r;sr2+=w*r*r;slv+=np.log(V)
    ah=sr/sw;asd=1/np.sqrt(sw);lo=ndtr((2-ah)/asd);hi=ndtr((4.5-ah)/asd);norm=hi-lo
    lp=-.5*(slv+sr2-sr*sr/sw)-.5*np.log(sw)+np.log(np.maximum(norm,1e-300))
    p=np.exp(lp-logsumexp(lp));pb=p.sum(1);ps=p.sum(0)
    post={'slope':quant(b,pb),'intrinsic_scatter_dex':quant(s,ps)}
    if intercept:
        def cdf(a):return np.sum(p*(ndtr((a-ah)/asd)-lo)/np.maximum(norm,1e-300))
        post['intercept_logkms_at_xcenter']=[brentq(lambda a:cdf(a)-q,2,4.5,xtol=1e-9) for q in Q]
    slope=np.dot(x-x.mean(),y-y.mean())/np.sum((x-x.mean())**2)
    resid=y-y.mean()-slope*(x-x.mean())
    su=dict(N=len(rows),quantiles=Q,posterior=post,OLS_slope=slope,OLS_intercept=y.mean()-slope*x.mean(),observed_OLS_residual_std_dex=np.std(resid,ddof=1),grid=[nb,ns],edge_mass={'slope':pb[:3].sum()+pb[-3:].sum(),'scatter_upper':ps[-3:].sum()})
    return su,dict(b=b,s=s,pb=pb,ps=ps)
def bxrows(rr):
    return [dict(source_id=r['id'],x=r['logLbol']-45,y=r['y'],cov=[[r['sx']**2,0],[0,r['sy']**2]]) for r in rr]
def checks():
    # Homoscedastic fixture recovers the known slope and vertical intrinsic SD.
    rng=np.random.default_rng(702);x=np.linspace(-1.5,1,200);y=3.3+.19*x+rng.normal(0,.05,len(x))
    rows=[dict(x=a,y=b,cov=[[1e-12,0],[0,1e-8]]) for a,b in zip(x,y)]
    su,g=posterior(rows);assert abs(su['posterior']['slope'][2]-.19)<.02
    assert abs(su['posterior']['intrinsic_scatter_dex'][2]-.05)<.01
    d,(axis,p)=difference(g,g);assert np.allclose(p,p[::-1],rtol=1e-12,atol=1e-15)
    assert abs(d['LRD_minus_nonLRD_dex'][2])<=.6/601
    return {'synthetic':su,'self_difference':d,'convolution_symmetric':True}
