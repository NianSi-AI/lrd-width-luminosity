"""Uniform fits for the 204-spectrum parent inventory. No rows are dropped for low quality.
Main width is ALWAYS the intrinsic single-Gaussian broad FWHM at default LSF.
Absorption is phenomenological positive transmission, accepted only by deltaBIC>6.
All stored fluxes are CGS; line-width/velocity units km/s. No population regression.
"""
from pathlib import Path
import json,hashlib,copy,time,traceback,sys,os,concurrent.futures,datetime
import numpy as np
from scipy.optimize import least_squares
from lrdwidth import linefit as c
O=Path('.')
W=Path('work/jwst')
LINES=['Halpha','Hbeta','Pabeta','Pagamma','HeI10830','OI8446']
PILOT={}

def native(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,dict):return {k:native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 return x

def null_model(d,absorption):
 # Same narrow/continuum/nuisance parameterization, with broad amplitude fixed zero.
 names,xx,lo,hi=c.setup(d,'single_gaussian',absorption,0);xx[5]=0;xx[6]=0;xx[7]=1000
 free=[i for i,n in enumerate(names) if n not in ['F_broad1','v_broad1','W_broad1']];best=None
 for seed in range(3):
  _,start,_,_=c.setup(d,'single_gaussian',absorption,seed);start=np.clip(start,lo+1e-8,hi-1e-8)
  def residual(p):
   x=xx.copy();x[free]=p;return (c.components(x,d,'single_gaussian',absorption)['total']-d['y'])/d['err']
  f=least_squares(residual,start[free],bounds=(lo[free],hi[free]),x_scale='jac',max_nfev=1200,ftol=2e-7,xtol=2e-7,gtol=2e-7)
  if best is None or np.sum(f.fun**2)<np.sum(best.fun**2):best=f
 x=xx.copy();x[free]=best.x;chi=float(np.sum(best.fun**2));return dict(absorption=absorption,chi2=chi,BIC=chi+len(free)*np.log(len(d['y'])),k=len(free),N=len(d['y']),names=names,params=x.tolist(),success=bool(best.success))

def metrics(f,d,null):
 co=c.components(f['params'],d,f['kind'],f['absorption']);p=dict(zip(f['names'],f['params']));m=c.widths(d['grid'],co['intrinsic_broad'])
 f['FWHM_intrinsic_kms']=float(p['W_broad1']) if f['kind']!='double_gaussian' else m['FWHM_kms'] if m else None
 f['ve_model_kms']=p['W_broad1']/(2*np.log(2)) if f['kind']=='exponential' else None
 f['delta_BIC_null']=null['BIC']-f['BIC']
 # Retain both the emitted analytic flux and actually transmitted finite-window flux.
 v=d['grid'];b=co['intrinsic_broad'];n=1000*p['F_narrow']*c.gauss(v,p['v_narrow'],p['W_narrow'])
 transmission=np.exp(-p['tau_abs']*np.exp(-.5*((v-p['v_abs'])/(p['W_abs']/c.SQ))**2)) if f['absorption'] else np.ones_like(v)
 mask=abs(v)<=5000
 f['broad_transmitted_flux_pm5000_erg_s_cm2']=float(np.trapezoid((b*transmission)[mask],v[mask])/1000*c.FUNIT)
 f['total_target_transmitted_flux_pm5000_erg_s_cm2']=float(np.trapezoid(((n+b)*transmission)[mask],v[mask])/1000*c.FUNIT)
 lsf=c.C/d['R_center']*d['lsf_scale'];conditions=dict(broad_SNR_ge5=(f['broad_SNR_conditional'] or 0)>=5,delta_BIC_null_ge6=f['delta_BIC_null']>=6,broad_FWHM_gt_LSF=f['FWHM_intrinsic_kms'] is not None and f['FWHM_intrinsic_kms']>lsf,all_parameters_away_from_bounds=not f['bound_parameters'],reduced_chi2_lt2=f['reduced_chi2']<2)
 f['quality_conditions']=conditions;f['quality_pass_five_conditions']=all(conditions.values());f['usable_numerical_fit']=f['success'] and np.isfinite(f['chi2']);f['quality_pass']=f['quality_pass_five_conditions'] and f['usable_numerical_fit'];f['LSF_FWHM_kms']=lsf
 return f

def fit_line(target,line):
 base=c.load_line(O/target['input_path'],target['z_fit'],line,1.)
 out=dict(file_id=target['file_id'],file=target['file'],source_id=target['source_id'],line=line,z_fit=target['z_fit'],status=base['status'],cases=[])
 if base['status']!='ready':out['coverage']=base;return out
 out['status']='fitted';out['data']={k:base[k] for k in ['v','lo','hi','y','err','quoted_err','R_center','grating','valid_fraction','center_um','window','Nvalid','wavelength_um']}
 cache=PILOT.get((target['file'],line));cache_ok=cache is not None and np.isclose(cache['literature']['z'],target['z_fit'],rtol=0,atol=1e-12) and np.array_equal(np.array(cache['data']['y']),base['y']) and np.array_equal(np.array(cache['data']['err']),base['err'])
 out['pilot_cache_candidates_reused']=bool(cache_ok)
 for scale in [1.,.7,1.3]:
  d=base if scale==1 else c.load_line(O/target['input_path'],target['z_fit'],line,scale)
  nulls=[null_model(d,False),null_model(d,True)];null=nulls[1] if nulls[0]['BIC']-nulls[1]['BIC']>6 else nulls[0]
  case=dict(lsf_FWHM_scale=scale,null_candidates=nulls,null=null,all_fits=[],selected_by_model=[],main=None,nonparametric_main=None)
  cached_case=next((q for q in cache['cases'] if q['lsf_FWHM_scale']==scale),None) if cache_ok else None
  for kind in c.MODELS:
   candidates=[]
   for absorption in [False,True]:
    cached=next((q for q in cached_case['all_fits'] if q['kind']==kind and q['absorption']==absorption),None) if cached_case else None
    f=copy.deepcopy(cached) if cached is not None else c.fit_model(d,kind,absorption,seeds=3)
    # New null and quality definitions are always recalculated, including for cached candidates.
    f=metrics(f,d,null);f['candidate_origin']='pilot_v3_identical_data_and_model' if cached is not None else 'full_v1_fit';candidates.append(f);case['all_fits'].append(f)
   no,ab=candidates;selected=ab if no['BIC']-ab['BIC']>6 else no;selected['absorption_BIC_improvement']=no['BIC']-ab['BIC'];case['selected_by_model'].append(selected)
  case['main']=case['selected_by_model'][0]
  # NP numbers remain ancillary; they never replace the prescribed Gaussian main width.
  main=case['main']
  if main['delta_BIC_null']>=6 and (main['broad_SNR_conditional'] or 0)>0:
   try:case['nonparametric_main']=c.nonparametric(d,main,nboot=80 if scale==1 and main['quality_pass'] else 0,seed=20261111)
   except Exception as e:case['nonparametric_failure']=dict(error=str(e),traceback=traceback.format_exc())
  out['cases'].append(case)
 return out

def one_file(target):
 start=time.monotonic();dest=W/(target['file_id']+'.json')
 if dest.exists():return dict(file_id=target['file_id'],status='already_saved',file=target['file'])
 output=dict(target=target,status='completed',lines=[],input_sha256=hashlib.sha256((O/target['input_path']).read_bytes()).hexdigest())
 for line in LINES:
  try:output['lines'].append(fit_line(target,line))
  except Exception as e:
   output['lines'].append(dict(file_id=target['file_id'],file=target['file'],source_id=target['source_id'],line=line,status='execution_failed',error=str(e),traceback=traceback.format_exc()));output['status']='completed_with_failures'
 output['elapsed_s']=time.monotonic()-start
 temp=dest.with_suffix('.writing');temp.write_text(json.dumps(native(output),ensure_ascii=False,indent=2));temp.rename(dest)
 return dict(file_id=target['file_id'],file=target['file'],status=output['status'],elapsed_s=round(output['elapsed_s'],1),fitted_lines=sum(x['status']=='fitted' for x in output['lines']),Ha_quality=next((x['cases'][0]['main']['quality_pass'] for x in output['lines'] if x['line']=='Halpha' and x['status']=='fitted'),None))
