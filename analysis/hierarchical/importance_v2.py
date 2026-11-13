"""Importance integration of the exact hierarchical posterior, not a MAP proxy."""
import numpy as np
from scipy.special import gammaln,logsumexp
from model import Q,native
DF=5

def tdraw(rng,mean,scale,n):
    ch=np.linalg.cholesky(scale);z=rng.normal(size=(n,len(mean)))@ch.T
    return mean+z/np.sqrt(rng.chisquare(DF,size=n)/DF)[:,None]

def tlogpdf(z,mean,scale):
    p=len(mean);r=z-mean;chol=np.linalg.cholesky(scale);v=np.linalg.solve(chol,r.T).T;mahal=np.sum(v*v,axis=1)
    return gammaln((DF+p)/2)-gammaln(DF/2)-p/2*np.log(DF*np.pi)-np.log(np.diag(chol)).sum()-(DF+p)/2*np.log1p(mahal/DF)

def wquant(x,w,q=Q):
    order=np.argsort(x);xx=x[order];ww=w[order];ww=ww/ww.sum();pos=np.cumsum(ww)-ww/2
    return np.interp(q,pos,xx)

def integrate(model,seed=712,n=8192,start=None,pilot=2048):
    assert n >= 4 and n % 2 == 0, "Equal proposal weights require an even draw count"
    rng=np.random.default_rng(seed);opt=model.map(start);mean=opt.x;scale=model.covariance_at(mean);scale=(scale+scale.T)/2+np.eye(len(mean))*1e-10
    z=tdraw(rng,mean,scale,pilot);lw=model.logp(z)-tlogpdf(z,mean,scale);w=np.exp(lw-logsumexp(lw));mean2=np.sum(z*w[:,None],axis=0);r=z-mean2;scale2=np.einsum('ni,n,nj->ij',r,w,r)*.8+.2*scale;scale2=(scale2+scale2.T)/2+np.eye(len(mean))*1e-9
    z=np.vstack((tdraw(rng,mean,scale,n//2),tdraw(rng,mean2,scale2,n-n//2)));logq=np.logaddexp(tlogpdf(z,mean,scale),tlogpdf(z,mean2,scale2))-np.log(2)
    logw=model.logp(z)-logq;w=np.exp(logw-logsumexp(logw));ess=1/np.sum(w*w);maxw=float(w.max());result={'method':'exact_posterior_importance_integration','N':model.n,'draws':n,'pilot_draws':pilot,'seed':seed,'ESS':float(ess),'max_weight':maxw,'diagnostic_pass':bool(ess>=400 and maxw<=.02),'MAP_success':bool(opt.success),'posterior':{}}
    for i,name in enumerate(model.param_names[:2+model.k]):
        vec=np.exp(z[:,i]) if name=='log_sigma' else z[:,i];result['posterior']['sigma_dex' if name=='log_sigma' else name]=wquant(vec,w)
    return result,z,w,opt.x
