"""Frozen-template injection recovery. Seeds and four scatter scenarios are retained from the original experiment.
Uses the original T2 functions without changing fit objective, bounds or seeds.
Widths km/s, flux CGS. Gaussian and two-sided exponential share true intrinsic FWHM.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,time,traceback,concurrent.futures,sys,copy
import numpy as np
from lrdwidth import linefit as c, fit_batch as fit
from templates import FROZEN as R, MU, template
from pathlib import Path
D=Path('work/injection/gaussian')
class reg:
 native=staticmethod(fit.native)

def one(args):
 pi,si,rep,i=args;dest=D/f'p{pi}s{si}r{rep:03d}i{i:02d}.json'
 if dest.exists():return dict(status='cached',name=dest.name)
 start=time.monotonic();r=R[i];base,p=template(i);d=dict(base);kind=['single_gaussian','exponential'][pi];rng=np.random.default_rng(np.random.SeedSequence([20261112,pi,si,rep,i]));scatter=[.1,.15,.2,.05][si]
 yy=MU+scatter*rng.normal() if si<3 else np.log10(2250)+.19*r['x_logL43']+scatter*rng.normal();width=10**yy;pars=p.copy();pars['W_broad1']=width;names=c.layout('Halpha',kind,False);pa=[pars[n] for n in names]
 noiseless=c.components(pa,d,kind,False)['total'];d['y']=noiseless+rng.normal(0,d['err'])
 try:
  nulls=[fit.null_model(d,a) for a in [False,True]];null=nulls[1] if nulls[0]['BIC']-nulls[1]['BIC']>6 else nulls[0]
  candidates=[c.fit_model(d,'single_gaussian',a,seeds=3) for a in [False,True]];no,ab=candidates;f=ab if no['BIC']-ab['BIC']>6 else no;f=fit.metrics(f,d,null)
  bounds=set(f['bound_parameters']);checks=dict(broadSNR=f['broad_SNR_conditional'] is not None and f['broad_SNR_conditional']>=5,deltaBIC=f['delta_BIC_null']>=6,widthLSF=f['FWHM_intrinsic_kms']>f['LSF_FWHM_kms'],widthBoundary='W_broad1' not in bounds,chi2=f['reduced_chi2']<2,converged=f['success'],absBoundary=not(f['absorption'] and bounds.intersection({'tau_abs','W_abs'})))
  names=f['names'];pp=np.array(f['params']);ii=[names.index('F_broad1'),names.index('W_broad1')];cov=np.array(f['covariance_conditional'])[np.ix_(ii,ii)];J=np.diag(1/(np.log(10)*pp[ii]));logcov=J@cov@J.T
  cov_ok=np.isfinite(logcov).all() and np.linalg.eigvalsh(logcov).min()>0
  out=dict(status='fitted',profile=kind,scenario=si,replicate=rep,template=i,source_id=r['source_id'],cls=r['lrd_class_stage2'],true_x=r['x_logL43'],true_y=yy,true_FWHM_kms=width,true_broad_flux_cgs=p['F_broad1']*c.FUNIT,true_narrow_flux_cgs=p['F_narrow']*c.FUNIT,selected=all(checks.values()),quality_checks=checks,recovered_x=r['x_logL43']+float(np.log10(f['broad_flux_erg_s_cm2']/(p['F_broad1']*c.FUNIT))),recovered_y=float(np.log10(f['FWHM_intrinsic_kms'])),cov_xy_log10=logcov.tolist(),covariance_PD=bool(cov_ok),main=f,candidates=candidates,null_candidates=nulls,noise_seed=[20261112,pi,si,rep,i],elapsed_s=time.monotonic()-start)
 except Exception as e:out=dict(status='execution_failed',profile=kind,scenario=si,replicate=rep,template=i,source_id=r['source_id'],error=str(e),traceback=traceback.format_exc(),elapsed_s=time.monotonic()-start)
 tmp=dest.with_suffix('.writing');tmp.write_text(json.dumps(reg.native(out),ensure_ascii=False,allow_nan=False));tmp.rename(dest);return dict(status=out['status'],name=dest.name,seconds=round(out['elapsed_s'],2))
def main():
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path('work/injection/gaussian'));p.add_argument('--nrep',type=int,default=100);p.add_argument('--limit',type=int);a=p.parse_args()
 c.self_checks()
 global D
 D=a.out;D.mkdir(parents=True,exist_ok=True)
 tasks=[(pi,si,rep,i) for rep in range(a.nrep) for pi in range(2) for si in range(4) for i in range(len(R))]
 if a.limit:tasks=tasks[:a.limit]
 for task in tasks:print(one(task),flush=True)
if __name__=='__main__':main()
