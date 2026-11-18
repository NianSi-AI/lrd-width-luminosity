"""Gaussian-population errors-in-variables regression.

Kelly (2007)-type hierarchy with one Gaussian predictor population, not the
flat latent-X likelihood. Log quantities are base10; their errors are dex.
All dimensions include actual source covariance, not a fitted noise inflation.
Independent NumPy/SciPy implementation.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import gaussian_kde

Q=np.array([.025,.16,.5,.84,.975])
def native(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [native(v) for v in x]
    return x

def save(path,obj):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(native(obj),f,ensure_ascii=False,indent=2,allow_nan=False)

class Model:
    def __init__(self,X,y,C,names,fixed=None,flat=False):
        self.X=np.asarray(X,float);self.y=np.asarray(y,float);self.C=np.asarray(C,float);self.names=names
        self.n,self.d=self.X.shape;self.flat=flat;self.fixed={} if fixed is None else dict(fixed)
        self.free=[j for j in range(self.d) if j not in self.fixed]
        self.k=len(self.free);self.sdscale=np.array([.1 if z=='t' else .7 for z in names]);self.muscale=np.array([.5 if z=='t' else 3. for z in names])
        self.bs=np.array([3. if z=='t' else 1. for z in names]);self.D=np.column_stack((self.X,self.y))
        assert self.n>=4 and self.C.shape==(self.n,self.d+1,self.d+1)
        assert np.isfinite(self.D).all() and np.linalg.eigvalsh(self.C).min()>-1e-10
        self.pairs=[(i,j) for i in range(1,self.d) for j in range(i)]
        self.ndim=2+self.k+(0 if flat else 2*self.d+len(self.pairs))
        self.param_names=['a']+[names[j] for j in self.free]+['log_sigma']
        if not flat:self.param_names+=['mu_'+x for x in names]+['logsd_'+x for x in names]+[f'partial_rho_{i}{j}' for i,j in self.pairs]
        self.bounds=[(2,4.5)]+[(-self.bs[j],self.bs[j]) for j in self.free]+[(-13,np.log(.6))]
        if not flat:self.bounds += [(-10,10)]*self.d+[(-10,3)]*self.d+[(-5,5)]*len(self.pairs)

    def unpack(self,z):
        z=np.atleast_2d(z);m=len(z);a=z[:,0];beta=np.zeros((m,self.d))
        for j,v in self.fixed.items():beta[:,j]=v
        beta[:,self.free]=z[:,1:1+self.k];logs=z[:,1+self.k];s=np.exp(logs)
        if self.flat:return a,beta,s,None,None
        i=2+self.k;mu=z[:,i:i+self.d];sd=np.exp(z[:,i+self.d:i+2*self.d]);raw=z[:,i+2*self.d:]
        L=np.zeros((m,self.d,self.d));L[:,0,0]=1.;idx=0
        for row in range(1,self.d):
            prod=np.ones(m)
            for col in range(row):
                r=np.tanh(raw[:,idx]);idx+=1;L[:,row,col]=r*prod;prod*=np.sqrt(1-r*r)
            L[:,row,row]=prod
        L*=sd[:,:,None];Omega=L@np.swapaxes(L,-1,-2)
        return a,beta,s,mu,Omega

    def logp(self,z):
        one=np.ndim(z)==1;z=np.atleast_2d(z);m=len(z);out=np.full(m,-np.inf)
        bnd=np.array(self.bounds);ok=np.all((z>bnd[:,0])&(z<bnd[:,1]),axis=1)
        if not ok.any():return out[0] if one else out
        zz=z[ok];a,beta,s,mu,Om=self.unpack(zz);prior=zz[:,1+self.k].copy() # uniform sigma, log-coordinate Jacobian
        if self.flat:
            V=s[:,None]**2+self.C[None,:,-1,-1]+np.einsum('mi,nij,mj->mn',beta,self.C[:,:-1,:-1],beta)-2*beta@self.C[:,:-1,-1].T
            res=self.y[None,:]-a[:,None]-beta@self.X.T
            ll=-.5*np.sum(np.log(2*np.pi*V)+res*res/V,axis=1)
        else:
            i=2+self.k;logsd=zz[:,i+self.d:i+2*self.d];sd=np.exp(logsd);raw=zz[:,i+2*self.d:]
            prior-=.5*np.sum((mu/self.muscale)**2,axis=1)
            # Half Student t3 scale; normalization constants cancel in all fits.
            prior+=np.sum(logsd-2*np.log1p(sd**2/(3*self.sdscale**2)),axis=1)
            for p,(_,col) in enumerate(self.pairs):
                alpha=1+(self.d-col-2)/2
                prior+=alpha*np.log(np.maximum(1-np.tanh(raw[:,p])**2,1e-300))
            means=np.concatenate((mu,(a+np.sum(beta*mu,axis=1))[:,None]),axis=1)
            joint=np.zeros((len(zz),self.d+1,self.d+1));joint[:,:-1,:-1]=Om
            ob=np.einsum('mij,mj->mi',Om,beta);joint[:,:-1,-1]=ob;joint[:,-1,:-1]=ob
            joint[:,-1,-1]=s*s+np.sum(beta*ob,axis=1)
            V=joint[:,None,:,:]+self.C[None,:,:,:]
            sign,ld=np.linalg.slogdet(V);res=self.D[None,:,:]-means[:,None,:]
            sol=np.linalg.solve(V,res[...,None])[...,0]
            ll=-.5*np.sum(ld+np.sum(res*sol,axis=-1)+(self.d+1)*np.log(2*np.pi),axis=1)
            ll[np.any(sign<=0,axis=1)]=-np.inf
        out[ok]=ll+prior
        return out[0] if one else out

    def initial(self):
        X=np.column_stack((np.ones(self.n),self.X[:,self.free]));y=self.y.copy()
        for j,v in self.fixed.items():y-=v*self.X[:,j]
        b=np.linalg.lstsq(X,y,rcond=None)[0]
        b[0]=np.clip(b[0],2.1,4.4)
        for k,j in enumerate(self.free):b[k+1]=np.clip(b[k+1],-.7*self.bs[j],.7*self.bs[j])
        resid=y-X@b;sig=np.clip(np.std(resid),.03,.3);z=list(b)+[np.log(sig)]
        if not self.flat:
            cov=np.atleast_2d(np.cov(self.X,rowvar=False));cov-=np.mean(self.C[:,:-1,:-1],axis=0)
            eig,U=np.linalg.eigh(cov);cov=(U*np.maximum(eig,1e-4))@U.T
            sd=np.sqrt(np.diag(cov));sd=np.maximum(sd,self.sdscale*.1);corr=cov/np.outer(sd,sd)
            corr=.8*corr+.2*np.eye(self.d);L=np.linalg.cholesky(corr)
            raw=[]
            for i in range(1,self.d):
                prod=1.
                for j in range(i):
                    r=np.clip(L[i,j]/prod,-.95,.95);raw.append(np.arctanh(r));prod*=np.sqrt(1-r*r)
            z+=list(np.mean(self.X,axis=0))+list(np.log(sd))+raw
        return np.array(z)

    def map(self,start=None):
        z=self.initial() if start is None else np.array(start);bb=[(l+1e-7,h-1e-7) for l,h in self.bounds]
        opt=minimize(lambda p:-self.logp(p),z,method='L-BFGS-B',bounds=bb,options={'maxiter':1800,'ftol':1e-11,'gtol':1e-6})
        if not opt.success:
            opt2=minimize(lambda p:-self.logp(p),opt.x,method='L-BFGS-B',bounds=bb,options={'maxiter':2400,'ftol':1e-10})
            if opt2.fun<opt.fun+1e-5:opt=opt2
        return opt

    def covariance_at(self,z):
        p=len(z);steps=2e-4*np.maximum(1,np.abs(z));H=np.zeros((p,p));f0=-self.logp(z)
        for i in range(p):
            ei=np.eye(p)[i]*steps[i]
            H[i,i]=(-self.logp(z+ei)+2*self.logp(z)-self.logp(z-ei))/(steps[i]**2)
            for j in range(i):
                ej=np.eye(p)[j]*steps[j]
                H[i,j]=H[j,i]=(-self.logp(z+ei+ej)+self.logp(z+ei-ej)+self.logp(z-ei+ej)-self.logp(z-ei-ej))/(4*steps[i]*steps[j])
        if not np.isfinite(H).all():return np.diag(np.full(p,.02))
        e,U=np.linalg.eigh(H);C=(U*(1/np.maximum(e,.05)))@U.T
        return C

def diagnostics(chain):
    # chain shape independent chain, draw, reported parameter.
    m,n,p=chain.shape;z=np.concatenate((chain[:,:n//2],chain[:,n//2:2*(n//2)]),axis=0);nn=z.shape[1]
    W=np.var(z,axis=1,ddof=1).mean(0);B=nn*np.var(z.mean(1),axis=0,ddof=1);v=(nn-1)/nn*W+B/nn
    rh=np.sqrt(v/W)
    # Mean per-chain autocorrelation, Geyer initial-positive paired sum.
    zz=chain-chain.mean(1,keepdims=True);size=1<<(2*n-1).bit_length();f=np.fft.rfft(zz,size,axis=1)
    ac=np.fft.irfft(f*np.conj(f),size,axis=1)[:,:n];ac/=np.arange(n,0,-1)[None,:,None];ac/=ac[:,0:1,:]
    ac=ac.mean(0);ess=[]
    for j in range(p):
        pair=ac[1:(n//2)*2-1:2,j]+ac[2:(n//2)*2:2,j];nonpos=np.flatnonzero(pair<=0);stop=int(nonpos[0]) if len(nonpos) else len(pair)
        tau=max(1.,1+2*np.sum(np.minimum.accumulate(pair[:stop])))
        ess.append(m*n/tau)
    return dict(split_Rhat=rh,ESS=np.array(ess),max_Rhat=float(rh.max()),min_ESS=float(min(ess)))

def sample(model,seed=102,steps=40000,burn=10000,thin=5,start=None):
    rng=np.random.default_rng(seed);opt=model.map(start);z0=opt.x;cov=model.covariance_at(z0);p=len(z0);m=4
    chol=np.linalg.cholesky(cov+1e-10*np.eye(p));z=np.repeat(z0[None,:],m,axis=0)
    for i in range(m):
        for _ in range(100):
            zz=z0+.4*chol@rng.normal(size=p)
            if np.isfinite(model.logp(zz)):z[i]=zz;break
    lp=model.logp(z);scale=2.38/np.sqrt(p);kept=[];hist=[];accepted=np.zeros(m);trial=np.zeros(m);scale_log=np.log(scale)
    for it in range(steps):
        small=rng.uniform()<.05;scl=.1 if small else np.exp(scale_log)
        proposal=z+scl*(rng.normal(size=(m,p))@chol.T);lpn=model.logp(proposal);acc=np.log(rng.uniform(size=m))<lpn-lp
        z[acc]=proposal[acc];lp[acc]=lpn[acc];accepted+=acc;trial+=1
        if it<burn:
            hist.append(z.copy())
            if (it+1)%250==0:
                frac=accepted/trial;scale_log+=np.clip((frac.mean()-.234)*.5,-.15,.15);accepted*=0;trial*=0
            if (it+1)%1000==0 and it>=1999:
                h=np.array(hist)[int(len(hist)*.25):].reshape(-1,p);cc=np.cov(h,rowvar=False)
                # Shrink towards curvature covariance to stabilize weak directions.
                cc=.9*cc+.1*cov+np.eye(p)*1e-9;chol=np.linalg.cholesky(cc)
        elif (it-burn)%thin==0:kept.append(z.copy())
    ch=np.transpose(np.array(kept),(1,0,2));physical=ch.copy();physical[:,:,1+model.k]=np.exp(ch[:,:,1+model.k]);diag=diagnostics(physical[:,:,:2+model.k])
    flat=physical.reshape(-1,p);result={'N':model.n,'method':'hierarchical_gaussian_population_MCMC' if not model.flat else 'flat_latent_predictor_MCMC','quantiles':Q,'posterior':{},'sampling':{'steps':steps,'burn':burn,'thin':thin,'chains':m,'seed':seed,'acceptance':accepted/np.maximum(trial,1),'diagnostics_regression':diag,'MAP_success':bool(opt.success),'MAP_message':str(opt.message)}}
    for k,name in enumerate(model.param_names[:2+model.k]):result['posterior']['sigma_dex' if name=='log_sigma' else name]=np.quantile(flat[:,k],Q)
    result['sampling']['parameter_diagnostics_names']=['a']+[model.names[j] for j in model.free]+['sigma_dex']
    result['converged']=bool(diag['max_Rhat']<1.05 and diag['min_ESS']>=400)
    if not model.flat:
        i=2+model.k
        result['population_posterior']={}
        for j,nm in enumerate(model.names):
            result['population_posterior']['mu_'+nm]=np.quantile(flat[:,i+j],Q)
            result['population_posterior']['sd_'+nm]=np.quantile(np.exp(flat[:,i+model.d+j]),Q)
    return result,ch,opt.x

def hpd_rank(chain,model,point=(.25,1),names=('x_ph','t')):
    ids=[model.param_names.index(n) for n in names];xy=chain.reshape(-1,chain.shape[-1])[:,ids]
    rng=np.random.default_rng(765);samplexy=xy[rng.choice(len(xy),min(3000,len(xy)),replace=False)]
    out={}
    for fac in [.7,1.,1.4]:
        kde=gaussian_kde(samplexy.T,bw_method=lambda obj:obj.scotts_factor()*fac);level=kde(np.array(point)[:,None])[0]
        dens=kde(samplexy.T);out[str(fac)]=float(np.mean(dens>=level))
    return {'rank':out['1.0'],'bandwidth_sensitivity':out,'method':'posterior KDE density rank; 3000 fixed-seed draws max','point':point}

def selfcheck():
    rng=np.random.default_rng(100923);X=rng.multivariate_normal([-.3,.02],[[.3,.01],[.01,.006]],400);true=np.array([.25,1.]);y=3.3+X@true+rng.normal(0,.06,len(X));C=np.repeat(np.diag([.03**2,.006**2,.015**2])[None],len(X),axis=0);obs=np.column_stack((X,y))+rng.multivariate_normal(np.zeros(3),C[0],len(X));model=Model(obs[:,:2],obs[:,2],C,['x_ph','t']);z=model.initial();lp=model.logp(z);assert np.isfinite(lp)
    opt=model.map();aa,bb,ss,mu,Om=model.unpack(opt.x)
    assert abs(bb[0,0]-.25)<.04 and abs(bb[0,1]-1)<.15 and abs(ss[0]-.06)<.02,(bb,ss)
    # Exact Gaussian marginal likelihood construction cross-check with scipy.
    from scipy.stats import multivariate_normal
    A=np.vstack((np.eye(2),bb[0]));mean=np.r_[mu[0],aa[0]+bb[0]@mu[0]];cov=A@Om[0]@A.T;cov[-1,-1]+=ss[0]**2
    direct=sum(multivariate_normal.logpdf(o,mean=mean,cov=cov+c) for o,c in zip(obs,C))
    V=cov[None]+C;res=obs-mean;ind=-.5*np.sum(np.linalg.slogdet(V)[1]+np.einsum('ni,ni->n',res,np.linalg.solve(V,res[...,None])[...,0])+3*np.log(2*np.pi))
    assert abs(direct-ind)<1e-8
    return {'passed':True,'synthetic_N':400,'truth':{'b':.25,'c':1.,'sigma':.06},'MAP':{'b':bb[0,0],'c':bb[0,1],'sigma':ss[0]},'marginal_loglik_difference':direct-ind}
