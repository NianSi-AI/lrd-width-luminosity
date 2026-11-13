"""Compiled implementation of the unchanged two-predictor Model target/MH sampler.

Only a numerical acceleration, intended for F partial/BC simulations. It does
not implement a new likelihood, prior, estimator, or scientific selection.
Independent pointwise comparisons with model.Model are required before use.
"""
import time
import numpy as np
try:
 from numba import njit
 NUMBA_AVAILABLE=True
except ImportError:
 NUMBA_AVAILABLE=False
 # Identical arithmetic without JIT acceleration. Optional numba speeds full runs.
 def njit(*args,**kwargs):
  return lambda function:function
from model import Model,Q,diagnostics
from diagnostics_v2 import all_diagnostics

@njit(cache=False)
def logp2(z,X,y,C,bs,sdscale,muscale,flat):
    a=z[0];b0=z[1];b1=z[2];ls=z[3]
    if not(2<a<4.5 and -bs[0]<b0<bs[0] and -bs[1]<b1<bs[1] and -13<ls<np.log(.6)):return -np.inf
    sig2=np.exp(2*ls);lp=ls
    if flat:
        for i in range(len(y)):
            v=sig2+C[i,2,2]+b0*b0*C[i,0,0]+b1*b1*C[i,1,1]+2*b0*b1*C[i,0,1]-2*b0*C[i,0,2]-2*b1*C[i,1,2]
            if v<=0:return -np.inf
            r=y[i]-a-b0*X[i,0]-b1*X[i,1];lp-=.5*(np.log(2*np.pi*v)+r*r/v)
        return lp
    mu0=z[4];mu1=z[5];ld0=z[6];ld1=z[7];raw=z[8]
    if not(-10<mu0<10 and -10<mu1<10 and -10<ld0<3 and -10<ld1<3 and -5<raw<5):return -np.inf
    s0=np.exp(ld0);s1=np.exp(ld1);rho=np.tanh(raw)
    lp-=.5*((mu0/muscale[0])**2+(mu1/muscale[1])**2)
    lp+=ld0-2*np.log1p(s0*s0/(3*sdscale[0]*sdscale[0]))+ld1-2*np.log1p(s1*s1/(3*sdscale[1]*sdscale[1]))+np.log(1-rho*rho)
    o00=s0*s0;o11=s1*s1;o01=s0*s1*rho;o02=o00*b0+o01*b1;o12=o01*b0+o11*b1;o22=sig2+b0*o02+b1*o12;my=a+b0*mu0+b1*mu1
    for i in range(len(y)):
        v00=o00+C[i,0,0];v01=o01+C[i,0,1];v11=o11+C[i,1,1];v02=o02+C[i,0,2];v12=o12+C[i,1,2];v22=o22+C[i,2,2]
        l00=np.sqrt(v00);l10=v01/l00;l20=v02/l00;tmp=v11-l10*l10
        if tmp<=0:return -np.inf
        l11=np.sqrt(tmp);l21=(v12-l20*l10)/l11;tmp=v22-l20*l20-l21*l21
        if tmp<=0:return -np.inf
        l22=np.sqrt(tmp);u0=(X[i,0]-mu0)/l00;u1=(X[i,1]-mu1-l10*u0)/l11;u2=(y[i]-my-l20*u0-l21*u1)/l22
        lp-=np.log(l00)+np.log(l11)+np.log(l22)+.5*(u0*u0+u1*u1+u2*u2)+1.5*np.log(2*np.pi)
    return lp

@njit(cache=False)
def run_mh(z0,cov,X,y,C,bs,sdscale,muscale,flat,seed,steps,burn,thin):
    np.random.seed(seed);p=len(z0);m=4;chol=np.linalg.cholesky(cov+np.eye(p)*1e-10);z=np.empty((m,p));lp=np.empty(m)
    for i in range(m):
        z[i]=z0
        for j in range(100):
            v=z0+.4*chol@np.random.normal(0,1,p)
            if np.isfinite(logp2(v,X,y,C,bs,sdscale,muscale,flat)):z[i]=v;break
        lp[i]=logp2(z[i],X,y,C,bs,sdscale,muscale,flat)
    logscale=np.log(2.38/np.sqrt(p));acc=np.zeros(m);trial=0;hist=np.empty((burn,m,p));nkeep=(steps-burn+thin-1)//thin;kept=np.empty((m,nkeep,p));ki=0
    for it in range(steps):
        scl=.1 if np.random.random()<.05 else np.exp(logscale)
        for j in range(m):
            v=z[j]+scl*chol@np.random.normal(0,1,p);lpn=logp2(v,X,y,C,bs,sdscale,muscale,flat)
            if np.log(np.random.random())<lpn-lp[j]:z[j]=v;lp[j]=lpn;acc[j]+=1
        trial+=1
        if it<burn:
            hist[it]=z
            if (it+1)%250==0:
                fraction=np.mean(acc)/trial;logscale+=min(.15,max(-.15,(fraction-.234)*.5));acc[:]=0;trial=0
            if (it+1)%1000==0 and it>=1999:
                start=(it+1)//4;h=hist[start:it+1].reshape(((it+1-start)*m,p));avg=np.zeros(p)
                for row in range(len(h)):avg+=h[row]
                avg/=len(h);cen=h-avg;cc=(cen.T@cen)/(len(h)-1);cc=.9*cc+.1*cov+np.eye(p)*1e-9;chol=np.linalg.cholesky(cc)
        elif (it-burn)%thin==0:
            for j in range(m):kept[j,ki]=z[j]
            ki+=1
    return kept,acc/max(trial,1)

def sample_compiled(model,seed=102,steps=40000,burn=10000,thin=5,start=None):
    assert model.d==2 and model.k==2 and not model.fixed
    opt=model.map(start);z0=opt.x;cov=model.covariance_at(z0);t=time.monotonic()
    ch,acc=run_mh(z0,cov,np.ascontiguousarray(model.X),np.ascontiguousarray(model.y),np.ascontiguousarray(model.C),model.bs,model.sdscale,model.muscale,model.flat,int(seed%(2**32-1)),steps,burn,thin)
    physical=ch.copy();physical[:,:,3]=np.exp(ch[:,:,3]);diag=diagnostics(physical[:,:,:4]);alldiag=all_diagnostics(ch,model.param_names);flat=physical.reshape(-1,model.ndim)
    out={'N':model.n,'method':'flat_latent_predictor_compiled_MCMC' if model.flat else 'hierarchical_gaussian_population_compiled_MCMC','quantiles':Q,'posterior':{name:np.quantile(flat[:,j],Q) for j,name in enumerate(['a']+[model.names[k] for k in model.free]+['sigma_dex'])},'sampling':{'steps':steps,'burn':burn,'thin':thin,'chains':4,'seed':seed,'acceptance':acc,'diagnostics_regression':diag,'diagnostics_all_parameters_rank_folded':alldiag,'MAP_success':bool(opt.success),'MAP_message':str(opt.message),'compiled_loop_seconds':time.monotonic()-t},'converged':bool(alldiag['max_Rhat']<1.01 and alldiag['min_bulk_ESS']>=400)}
    return out,ch,z0

def selfcheck():
    rng=np.random.default_rng(8928);diff=[]
    for flat in [False,True]:
        for names in [['line','x_ph'],['x_ph','t']]:
            X=rng.normal(size=(11,2));y=3.2+.2*X[:,0]+rng.normal(0,.08,11);C=np.repeat(np.array([[.01,0,.001],[0,.003,0],[.001,0,.005]])[None],11,axis=0);m=Model(X,y,C,names,flat=flat);z=m.initial()
            for _ in range(100):
                v=z+rng.normal(0,.15,len(z));expected=m.logp(v);actual=logp2(v,X,y,C,m.bs,m.sdscale,m.muscale,flat)
                assert np.isfinite(expected)==np.isfinite(actual)
                if np.isfinite(expected):diff.append(abs(expected-actual))
    assert max(diff)<1e-8,max(diff)
    return {'passed':True,'pointwise_target_checks':len(diff),'max_absolute_logp_difference':max(diff),'exact_same_prior_and_likelihood':True,'numba_acceleration_available':NUMBA_AVAILABLE}

@njit(cache=False)
def logp2_batch(z,X,y,C,bs,sdscale,muscale,flat):
    out=np.empty(len(z))
    for i in range(len(z)):out[i]=logp2(z[i],X,y,C,bs,sdscale,muscale,flat)
    return out

class CompiledModel(Model):
    """Drop-in Model for exactly two free predictors; unchanged scientific target."""
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        assert self.d==2 and self.k==2 and not self.fixed
    def logp(self,z):
        zz=np.asarray(z,float)
        if zz.ndim==1:return logp2(zz,self.X,self.y,self.C,self.bs,self.sdscale,self.muscale,self.flat)
        return logp2_batch(zz,self.X,self.y,self.C,self.bs,self.sdscale,self.muscale,self.flat)
