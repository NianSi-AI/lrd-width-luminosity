"""Rank-normalized/folded split Rhat and Geyer bulk/tail ESS.

Vehtari et al. style diagnostics; no external MCMC package required.
Tail ESS uses 2.5/97.5 percent indicators and is reported separately from
median-oriented bulk ESS. Chains are genuinely independent MH runs.
"""
import numpy as np
from scipy.special import ndtri
from scipy.stats import rankdata

def split(ch):
    n=ch.shape[1]//2
    return np.concatenate((ch[:,:n],ch[:,n:2*n]),axis=0)

def ranknorm(x):
    return ndtri((rankdata(x.ravel()).reshape(x.shape)-.375)/(x.size+.25))

def rhat(x):
    n=x.shape[1];W=np.var(x,axis=1,ddof=1).mean();vp=(n-1)/n*W+np.var(np.mean(x,axis=1),ddof=1)
    return float(np.sqrt(vp/W)) if W>0 else (1. if vp==0 else float('inf'))

def ess(x):
    m,n=x.shape;zz=x-x.mean(1,keepdims=True);size=1<<(2*n-1).bit_length();f=np.fft.rfft(zz,size,axis=1);ac=np.fft.irfft(f*np.conj(f),size,axis=1)[:,:n]/n
    W=np.var(x,axis=1,ddof=1).mean();vp=(n-1)/n*W+np.var(x.mean(1),ddof=1)
    if vp<=0:return float(m*n)
    rho=1-(W-ac.mean(0))/vp;rho[0]=1.;pair=rho[:2*(n//2):2]+rho[1:2*(n//2):2];bad=np.flatnonzero(pair<=0);stop=bad[0] if len(bad) else len(pair)
    tau=max(1.,-1+2*np.sum(np.minimum.accumulate(pair[:stop])))
    return float(m*n/tau)

def all_diagnostics(chain,names=None):
    ch=split(np.asarray(chain,float));p=ch.shape[-1];out={}
    for j in range(p):
        x=ch[:,:,j];zn=ranknorm(x);zf=ranknorm(np.abs(x-np.median(x)));lo,hi=np.quantile(x,[.025,.975]);lower=ess((x<=lo).astype(float));upper=ess((x>=hi).astype(float))
        out[names[j] if names is not None else str(j)]={'rank_folded_Rhat':max(rhat(zn),rhat(zf)),'bulk_ESS':ess(zn),'tail_lower_ESS':lower,'tail_upper_ESS':upper,'tail_ESS':min(lower,upper)}
    return {'parameters':out,'max_Rhat':max(v['rank_folded_Rhat'] for v in out.values()),'min_bulk_ESS':min(v['bulk_ESS'] for v in out.values()),'min_tail_ESS':min(v['tail_ESS'] for v in out.values())}
