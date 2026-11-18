"""M4 exact seeded spectra, five independently measured width definitions.

Gaussian recoveries are reused from the existing M4 cache after checking seed,
injected width, flux and source identity. Other fits use original T2 optimization.
Primary selection follows M4 Gaussian QC, identical membership for all widths.
An all-29-template diagnostic is also stored; none is labeled a population model.
"""
from measure import *
from lrdwidth import fit_batch as fit
import concurrent.futures,os,time,traceback,argparse
from templates import FROZEN, MU, template
M=Path('work/injection')
KINDS=['single_gaussian','lorentzian','exponential','nonparametric_W50','nonparametric_FWHM']

def one(args):
    pi,si,rep,i=args;name=f'p{pi}s{si}r{rep:03d}i{i:02d}.json';dest=O/'W2_measurements'/name
    if dest.exists():return read(dest)
    t0=time.monotonic();r=FROZEN[i];base,p=template(i);d=dict(base);kind=['single_gaussian','exponential'][pi]
    seed=[20261112,pi,si,rep,i];rng=np.random.default_rng(np.random.SeedSequence(seed));scatter=.05 if si==3 else .15
    ytrue=np.log10(2250)+.19*r['x_logL43']+scatter*rng.normal() if si==3 else MU+scatter*rng.normal()
    pars=p.copy();pars['W_broad1']=10**ytrue;pa=[pars[n] for n in core.layout('Halpha',kind,False)]
    noiseless=core.components(pa,d,kind,False)['total'];d['y']=noiseless+rng.normal(0,d['err'])
    cached=read(M/'gaussian'/name)
    assert cached['noise_seed']==seed and cached['source_id']==r['source_id'] and cached['status']=='fitted'
    assert cached['true_y']==ytrue and cached['true_broad_flux_cgs']==p['F_broad1']*core.FUNIT
    results={};chosen={};npstatus={}
    try:
        g=cached['main'];chosen['single_gaussian']=g
        null=min(cached['null_candidates'],key=lambda a:a['BIC'])
        # null choice changes only detection metadata; broad fit likelihood unaffected.
        for model in ['lorentzian','exponential']:
            no=core.fit_model(d,model,False,seeds=3);ab=core.fit_model(d,model,True,seeds=3)
            f=ab if no['BIC']-ab['BIC']>6 else no;chosen[model]=fit.metrics(f,d,null)
        for model,f in chosen.items():
            row=make_row(r,log_measurement(f));row.update(converged=f['success'],bound_parameters=f['bound_parameters'],absorption=f['absorption'],BIC=f['BIC']);results[model]=row
        npfit=core.nonparametric(d,g,nboot=0,seed=0)
        bs,_=np_bootstrap(d,npfit,200,np.random.SeedSequence(seed+[10,2]))
        for model,metric in [('nonparametric_W50','W50_kms'),('nonparametric_FWHM','FWHM_kms')]:
            if npfit['metrics'][metric] is None:raise ValueError('Nonparametric '+metric+' missing')
            mm=dict(width_kms=npfit['metrics'][metric],flux_cgs=npfit['metrics']['flux_window_erg_s_cm2'],cov_log_flux_width=bs[model]['cov_log_flux_width'])
            results[model]=make_row(r,mm);npstatus[model]=bs[model]
        out=dict(status='complete',source_id=r['source_id'],profile=kind,scenario=si,rep=rep,template=i,seed=seed,true_y=ytrue,M4_Gaussian_selected=cached['selected'],rows=results,np_covariance_scope='fixed-alpha fixed-subtraction 200 draws',np_bootstrap=npstatus,old_Gaussian_cache=str(M/'gaussian'/name),selected_fits={k:{a:v[a] for a in ['kind','absorption','names','params','covariance_conditional','BIC','success','bound_parameters']} for k,v in chosen.items()})
    except Exception as e:
        out=dict(status='failed',source_id=r['source_id'],profile=kind,scenario=si,rep=rep,template=i,seed=seed,rows=results,error=str(e),traceback=traceback.format_exc(),M4_Gaussian_selected=cached['selected'])
    out['elapsed_s']=time.monotonic()-t0;save(dest,out);return out
def main():
 p=argparse.ArgumentParser();p.add_argument('--work',type=Path,default=Path('work/injection'));p.add_argument('--scenario',type=int,choices=[1,3],required=True);p.add_argument('--nrep',type=int,default=100);p.add_argument('--limit',type=int);a=p.parse_args()
 core.self_checks()
 global M,O
 M=O=a.work.resolve();(O/'W2_measurements').mkdir(parents=True,exist_ok=True)
 ids=[i for i,r in enumerate(FROZEN) if r['lrd_class_stage2']=='LRD'];assert len(ids)==29
 tasks=[(pi,a.scenario,rep,i) for rep in range(a.nrep) for pi in range(2) for i in ids]
 if a.limit:tasks=tasks[:a.limit]
 for job in tasks:print(one(job)['status'],job,flush=True)
if __name__=='__main__':main()
