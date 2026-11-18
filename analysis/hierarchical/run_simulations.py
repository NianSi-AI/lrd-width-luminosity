"""Task11a F6/F7: fixed-design measurement-error calibration.

All log variables are dex: line=log L_Halpha/1e43, ph=log Lph/1e45,
t=log T/5000, y=log width/(km/s). The exact same posterior Model is
integrated for observation and every simulated dataset. No MAP statistic.
"""
from pathlib import Path
import os,sys,json,time,hashlib,argparse,traceback
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import numpy as np
HERE=Path(__file__).resolve().parent;HIER=HERE;sys.path.insert(0,str(HIER))
from lrdwidth.paths import DATA
from model import native
from numba_sampler import CompiledModel as Model
from robust_fit_v3 import fit
SAMPLES=['main_E','Curtis_A_own','Curtis_B_own']
SCATTERS=['measured','0.05','0.075','0.10'];HYPS=['ph','line'];METHODS=['hierarchical','flat']
GLOBAL={}

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save_new(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(native(x),f,ensure_ascii=False,separators=(',',':'),allow_nan=False)

def prepare(data):
    out={}
    for key in SAMPLES:
        rows=data[key];X=np.array([[r['line_x'],r['x'],r['t']] for r in rows]);y=np.array([r['y'] for r in rows]);C=np.zeros((len(rows),4,4))
        for i,r in enumerate(rows):
            C[i]=np.diag([r['line_sx']**2,r['sx']**2,r['st']**2,r['sy']**2]);C[i,0,3]=C[i,3,0]=r['line_width_cov']
        ev,U=np.linalg.eigh(C);assert ev.min()>-1e-10;L=U*np.sqrt(np.maximum(ev,0))[:,None,:]
        assert np.allclose(L@np.swapaxes(L,-1,-2),C,rtol=1e-9,atol=1e-12)
        out[key]={'X':X,'y':y,'C':C,'noise_factor':L,'ids':[r.get('source_id',r['temperature_id']) for r in rows]}
    return out

def init(data_path,base_path,outpath):
    GLOBAL.update(data=prepare(read(data_path)),base=read(base_path),out=Path(outpath),data_sha=sha(data_path),base_sha=sha(base_path))

def fixture(job):
    si,sci,hi,rep=job;key=SAMPLES[si];hyp=HYPS[hi];scenario=SCATTERS[sci];d=GLOBAL['data'][key];base=GLOBAL['base'][key][hyp]
    ss=np.random.SeedSequence([2026111211,si,sci,hi,rep]);children=ss.spawn(5);rng=np.random.default_rng(children[0]);infseeds=[int(s.generate_state(1)[0]) for s in children[1:]]
    sigma=base['sigma'] if scenario=='measured' else float(scenario);X=d['X']
    law=base['a']+.25*X[:,1]+X[:,2] if hyp=='ph' else base['a']+base['b_line']*X[:,0]
    intrinsic=rng.normal(0,sigma,len(X));error=np.einsum('nij,nj->ni',d['noise_factor'],rng.normal(size=(len(X),4)))
    obs=np.column_stack((X,law+intrinsic))+error
    return d,obs,law,sigma,infseeds

def run_one(job):
    si,sci,hi,rep=job;key=SAMPLES[si];hyp=HYPS[hi];scenario=SCATTERS[sci];dest=GLOBAL['out']/'replicates'/key/scenario/hyp/f'{rep:04d}.json'
    if dest.exists():
        r=read(dest);assert r['input_sha256']==GLOBAL['data_sha'] and r['baseline_sha256']==GLOBAL['base_sha'],'Existing replicate has different frozen input'
        return (job,True,r.get('all_fits_pass',False),0.,'cached')
    t=time.monotonic();res={'sample':key,'N':len(GLOBAL['data'][key]['X']),'scenario':scenario,'hypothesis':hyp,'replicate':rep,'seed_sequence':[2026111211,si,sci,hi,rep],'input_sha256':GLOBAL['data_sha'],'baseline_sha256':GLOBAL['base_sha'],'algorithm_version':'v3_compiled_same_target','algorithm_hashes':{'compiled':sha(HIER/'numba_sampler.py'),'robust':sha(HIER/'robust_fit_v3.py')},'fits':{}}
    try:
        d,obs,law,sigma,seeds=fixture(job);res.update(sigma_truth=sigma,simulated_observations=obs,observed_columns=['line','x_ph','t','y'],generating_mean_y=law,inference_seeds=seeds)
        for mi,method in enumerate(METHODS):
            idx=[0,1,3];C=d['C'][:,idx][:,:,idx];model=Model(obs[:,:2],obs[:,3],C,['line','x_ph'],flat=method=='flat')
            try:res['fits'][method+'_partial']=fit(model,seed=seeds[mi],draws=4096)
            except Exception as e:res['fits'][method+'_partial']={'diagnostic_pass':False,'method':'exception','error':type(e).__name__+': '+str(e),'traceback':traceback.format_exc()}
        if hyp=='line':
            for mi,method in enumerate(METHODS):
                idx=[1,2,3];C=d['C'][:,idx][:,:,idx];model=Model(obs[:,1:3],obs[:,3],C,['x_ph','t'],flat=method=='flat')
                try:res['fits'][method+'_bc']=fit(model,seed=seeds[2+mi],draws=4096)
                except Exception as e:res['fits'][method+'_bc']={'diagnostic_pass':False,'method':'exception','error':type(e).__name__+': '+str(e),'traceback':traceback.format_exc()}
        res['all_fits_pass']=all(v.get('diagnostic_pass',False) for v in res['fits'].values())
    except Exception as e:res.update(all_fits_pass=False,error=type(e).__name__+': '+str(e),traceback=traceback.format_exc())
    res['seconds']=time.monotonic()-t;save_new(dest,res)
    return(job,False,res['all_fits_pass'],res['seconds'],'ok')

def selfcheck():
    # Test measurement-noise construction independently of astrophysical result.
    rng=np.random.default_rng(219);C=np.array([[.04,0,0,.002],[0,.01,0,0],[0,0,.001,0],[.002,0,0,.003]])
    ev,U=np.linalg.eigh(C);L=U*np.sqrt(ev)[None,:];noise=rng.normal(size=(100000,4))@L.T;est=np.cov(noise,rowvar=False)
    assert np.max(np.abs(np.diag(est)/np.diag(C)-1))<.02 and abs(est[0,3]-C[0,3])<.0002
    assert np.allclose(L@L.T,C)
    return {'passed':True,'known_covariance':C,'simulated_covariance':est,'draws':100000,'note':'measurement noise includes independent errors on ph and t, plus same-fit line-y covariance'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--selfcheck-only',action='store_true');ap.add_argument('--data',default=str(DATA/'photosphere/hierarchical/datasets.json'));ap.add_argument('--baseline',default=str(DATA/'photosphere/hierarchical/generating_baselines.json'));ap.add_argument('--out',default='work/hierarchical');ap.add_argument('--workers',type=int,default=6);ap.add_argument('--n',type=int,default=2000);ap.add_argument('--limit',type=int,default=None);args=ap.parse_args();
    if args.selfcheck_only:print(native(selfcheck()));return
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    data=read(args.data);base=read(args.baseline);assert all(base[k]['precision_target_pass'] and base[k]['same_members'] for k in SAMPLES)
    plan={'sample_order':SAMPLES,'scenario_order':SCATTERS,'hypothesis_order':HYPS,'n_per_cell':args.n,'data_path':str(Path(args.data).resolve()),'data_sha256':sha(args.data),'baseline_path':str(Path(args.baseline).resolve()),'baseline_sha256':sha(args.baseline),'methods':METHODS,'workers':args.workers,'model_sha256':sha(HIER/'model.py'),'importance_sha256':sha(HIER/'importance_v2.py'),'robust_fit_sha256':sha(HIER/'robust_fit_v3.py'),'compiled_target_sha256':sha(HIER/'numba_sampler.py'),'selfcheck':selfcheck()}
    if (out/'execution_manifest_v3.json').exists():
        old=read(out/'execution_manifest_v3.json');assert all(old[k]==plan[k] for k in ['data_sha256','baseline_sha256','model_sha256','importance_sha256','robust_fit_sha256','n_per_cell'])
    else:save_new(out/'execution_manifest_v3.json',plan)
    # Interleave samples/hypotheses: early diagnostics do not only cover one cell.
    jobs=[(si,sci,hi,r) for r in range(args.n) for sci in range(4) for hi in range(2) for si in range(3)]
    if args.limit is not None:jobs=jobs[:args.limit]
    t=time.monotonic();done=0;failed=0;cached=0;iterator=iter(jobs)
    with ProcessPoolExecutor(max_workers=args.workers,initializer=init,initargs=(args.data,args.baseline,args.out)) as pool:
        pending={}
        for _ in range(args.workers*2):
            try:j=next(iterator)
            except StopIteration:break
            pending[pool.submit(run_one,j)]=j
        while pending:
            ready,_=wait(pending,return_when=FIRST_COMPLETED)
            for f in ready:
                job=pending.pop(f)
                try:
                    _,iscached,passed,seconds,msg=f.result();cached+=iscached;failed+=not passed
                except Exception as e:
                    failed+=1;print('WORKER_EXCEPTION',job,type(e).__name__,str(e),flush=True)
                done+=1
                if done%100==0 or done==len(jobs):print('progress',done,'/',len(jobs),'allfit_failures',failed,'cached',cached,'elapsed_min',round((time.monotonic()-t)/60,2),flush=True)
                try:j=next(iterator)
                except StopIteration:continue
                pending[pool.submit(run_one,j)]=j
    final={'jobs_finished':done,'jobs_planned':len(jobs),'allfit_failure_count':failed,'cached':cached,'wall_seconds':time.monotonic()-t}
    save_new(out/('execution_completed_v3.json' if args.limit is None else f'pilot_{args.limit}_completed.json'),final);print(native(final),flush=True)
if __name__=='__main__':main()
