"""Preserve independent successful definitions when another measurement fails.

No criterion is changed: original >=160/200 finite bootstrap realizations and
closed central half-maximum crossings are required. Failed attempts remain in
W2_measurements. This repair reconstructs NP draws to retain W50 diagnostics;
it never substitutes censored FWHM or drops a source to force a fit.
"""
from width_recovery import *
from scipy.optimize import nnls
from lrdwidth.paths import DATA

def repair(path):
    old=read(path);dest=O/'W2_np_diagnostics'/path.name
    if dest.exists():return read(dest)
    pi=0 if old['profile']=='single_gaussian' else 1;si,rep,i=old['scenario'],old['rep'],old['template'];r=FROZEN[i]
    base,p=template(i);d=dict(base);seed=old['seed'];rng=np.random.default_rng(np.random.SeedSequence(seed))
    ytrue=np.log10(2250)+.19*r['x_logL43']+.05*rng.normal() if si==3 else MU+.15*rng.normal()
    pa=p.copy();pa['W_broad1']=10**ytrue;pars=[pa[n] for n in core.layout('Halpha',old['profile'],False)]
    d['y']=core.components(pars,d,old['profile'],False)['total']+rng.normal(0,d['err'])
    cached=read(M/'gaussian'/path.name);g=cached['main'];out={'rows':{},'valid':{},'counts':{},'failure':{}}
    try:
        n=core.nonparametric(d,g,nboot=0,seed=0);knots=np.array(n['knots_kms']);sol=np.array(n['intrinsic_knots']);fine=np.array(n['v_fine'])
        basis=np.maximum(1-abs(d['grid'][:,None]-knots)/125,0);A=np.asarray(d['K']@basis);D2=np.diff(np.eye(len(knots)),n=2,axis=0);alpha=n['cv_alpha_grid'][n['selected_alpha_index']]
        mat=np.vstack([A/d['err'][:,None],np.sqrt(alpha)*D2]);pred=A@sol;rng=np.random.default_rng(np.random.SeedSequence(seed+[10,2]));samples=[]
        for _ in range(200):
            yy=pred+rng.normal(0,d['err']);xx=nnls(mat,np.r_[yy/d['err'],np.zeros(len(D2))],maxiter=3000)[0];m=core.widths(fine,np.interp(fine,knots,xx))
            samples.append([m['flux_window_erg_s_cm2'],m['W50_kms'],m['FWHM_kms'] or np.nan])
        a=np.array(samples)
        for j,kind,metric in [(1,'nonparametric_W50','W50_kms'),(2,'nonparametric_FWHM','FWHM_kms')]:
            mask=np.isfinite(a[:,[0,j]]).all(1)&(a[:,[0,j]]>0).all(1);N=int(mask.sum());central=n['metrics'][metric]
            valid=N>=160 and central is not None;out['valid'][kind]=valid;out['counts'][kind]=N
            if central is not None and N>=20:
                C=np.cov(np.log10(a[mask][:,[0,j]]),rowvar=False,ddof=1)
                out['rows'][kind]=make_row(r,dict(width_kms=central,flux_cgs=n['metrics']['flux_window_erg_s_cm2'],cov_log_flux_width=C))
            if not valid:out['failure'][kind]={'central_width_missing':central is None,'N_valid_bootstrap':N,'minimum_required':160}
    except Exception as e:out['error']=str(e)
    save(dest,out);return out
def set_work(path):
    global O,M
    O=M=Path(path).resolve()

def finalize(si):
    paths=sorted((O/'W2_measurements').glob(f'p*s{si}r*i*.json'));assert len(paths)==5800
    failed=[p for p in paths if read(p)['status']!='complete'];(O/'W2_np_diagnostics').mkdir(exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=6,initializer=set_work,initargs=(str(O),)) as ex:list(ex.map(repair,failed))
    bykey={}
    for p in paths:
        r=read(p)
        if r['status']!='complete':
            d=read(O/'W2_np_diagnostics'/p.name)
            r['rows'].update({k:v for k,v in d['rows'].items() if d['valid'].get(k,False)})
            r['np_invalid']=d.get('failure',{})
        bykey[(r['profile'],r['rep'],r['template'])]=r
    stat=[];ids=[i for i,r in enumerate(FROZEN) if r['lrd_class_stage2']=='LRD']
    for profile in ['single_gaussian','exponential']:
        for rep in range(100):
            batch=[bykey[profile,rep,i] for i in ids]
            for selection in ['M4_Gaussian_QC','all29']:
                subset=[r for r in batch if selection=='all29' or r['M4_Gaussian_selected']]
                for kind in KINDS:
                    missing=[r['source_id'] for r in subset if kind not in r['rows'] or kind in r.get('np_invalid',{})]
                    row=dict(scenario=si,profile=profile,rep=rep,selection=selection,definition=kind,N=len(subset),missing_sources=missing,status='complete')
                    if missing or len(subset)<4:row['status']='incomplete';stat.append(row);continue
                    su,g=posterior([r['rows'][kind] for r in subset]);row.update(sigma_median=su['posterior']['intrinsic_scatter_dex'][2],sigma_quantiles=su['posterior']['intrinsic_scatter_dex'],OLS_SD=su['observed_OLS_residual_std_dex'],slope_median=su['posterior']['slope'][2]);stat.append(row)
            if rep%20==0:print('W2 finalizing',si,profile,rep,flush=True)
    save(f'W2_replicates_by_definition_s{si}.json',stat)
    w1=read(DATA/'widths/injection_observed.json');summary=[]
    for profile in ['single_gaussian','exponential']:
        for sel in ['M4_Gaussian_QC','all29']:
            for kind in KINDS:
                allr=[r for r in stat if r['profile']==profile and r['selection']==sel and r['definition']==kind];rs=[r for r in allr if r['status']=='complete'];obs=w1[kind]['own_luminosity']['LRD'];inj=.05 if si==3 else .15
                row=dict(scenario=si,profile=profile,selection=sel,definition=kind,Nrequested=100,Nvalid=len(rs),Nincomplete=100-len(rs),incomplete_replicates=[r['rep'] for r in allr if r['status']!='complete'])
                if rs:
                    for metric,val in [('sigma_median',obs['posterior']['intrinsic_scatter_dex'][2]),('OLS_SD',obs['observed_OLS_residual_std_dex'])]:
                        arr=np.array([r[metric] for r in rs]);row[metric]={'p16_p50_p84':np.quantile(arr,[.16,.5,.84]),'median_minus_injected':np.median(arr)-inj,'observed_W1':val,'count_le_observed':int((arr<=val).sum()),'denominator':len(rs)}
                    row['N_selected_p16_p50_p84']=np.quantile([r['N'] for r in rs],[.16,.5,.84])
                summary.append(row)
    save(f'W2_final_summary_s{si}.json',summary)
    save(f'W2_finalization_QA_s{si}.json',{'N_spectra':len(paths),'N_initial_partial':len(failed),'no_source_dropped_to_complete_definition':True,'minimum_usable_bootstraps_retained':160,'supersedes_for_reporting':f'W2_summary_s{si}.json','scope':'Definition-specific availability; missing spectra retain their fixed selection membership and invalidate only that definition/replicate.'})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path,default=Path('work/injection'));p.add_argument('--scenario',type=int,choices=[1,3],required=True);args=p.parse_args();O=M=args.work.resolve();finalize(args.scenario)
