"""Summarize all fixed F replicate denominators without deleting failures."""
from pathlib import Path
import json,sys,collections,hashlib
import numpy as np
from scipy.stats import beta
O=Path(sys.argv[1]) if len(sys.argv)>1 else Path('work/hierarchical')
from lrdwidth.paths import DATA
from model import native,save
SAMPLES=['main_E','Curtis_A_own','Curtis_B_own'];SCENARIOS=['measured','0.05','0.075','0.10'];METHODS=['hierarchical','flat'];Q=[.025,.16,.5,.84,.975]

def read(p):return json.loads(Path(p).read_text())
def interval(k,n):
    return [0. if k==0 else float(beta.ppf(.025,k,n-k+1)),1. if k==n else float(beta.ppf(.975,k+1,n-k))]
def tail(x,threshold,side):
    k=int(np.sum(x<=threshold if side=='le' else x>=threshold));n=len(x)
    return {'count':k,'denominator':n,'fraction':k/n,'binomial_95_interval':interval(k,n),'comparison':side,'threshold':float(threshold)}
def summary(values):return {'quantiles':Q,'values':np.quantile(values,Q),'N':len(values)}

def main():
    obs=read(DATA/'photosphere/hierarchical/observed_thresholds.json');base=read(DATA/'photosphere/hierarchical/generating_baselines.json');data=read(DATA/'photosphere/hierarchical/datasets.json');manifest=read(O/'execution_manifest_v3.json');planned=manifest['n_per_cell'];assert planned>=2000
    assert hashlib.sha256((DATA/'photosphere/hierarchical/datasets.json').read_bytes()).hexdigest()==manifest['data_sha256']
    assert hashlib.sha256((DATA/'photosphere/hierarchical/generating_baselines.json').read_bytes()).hexdigest()==manifest['baseline_sha256']
    allout={};failures=[];methods=collections.Counter();attempts=collections.Counter();algorithms=collections.Counter();ess=[];maxw=[];durations=[];source_counts={}
    for name in SAMPLES:
        allout[name]={'N_sources':len(data[name]),'source_ids':[r.get('source_id',r['temperature_id']) for r in data[name]],'scenarios':{}}
        for scenario in SCENARIOS:
            records={};scenario_result={'sigma_ph':base[name]['ph']['sigma'] if scenario=='measured' else float(scenario),'sigma_line':base[name]['line']['sigma'] if scenario=='measured' else float(scenario),'methods':{}}
            for hyp in ['ph','line']:
                records[hyp]=[]
                for ri in range(planned):
                    path=O/'replicates'/name/scenario/hyp/f'{ri:04d}.json'
                    if not path.exists():failures.append({'sample':name,'scenario':scenario,'hypothesis':hyp,'replicate':ri,'failure':'missing file'});records[hyp].append(None);continue
                    r=read(path);assert r['input_sha256']==manifest['data_sha256'] and r['baseline_sha256']==manifest['baseline_sha256'];records[hyp].append(r);algorithms[r.get('algorithm_version','v2_original_target')]+=1;durations.append(r['seconds'])
                    for fitname,v in r['fits'].items():
                        methods[v['method']]+=1;attempts[str(len(v.get('attempts',[])))]+=1
                        if 'ESS' in v:ess.append(v['ESS']);maxw.append(v['max_weight'])
                        if not v.get('diagnostic_pass',False):failures.append({'sample':name,'scenario':scenario,'hypothesis':hyp,'replicate':ri,'fit':fitname,'failure':'posterior numerical diagnostics failed','path':str(path.relative_to(O)),'attempts':v.get('attempts',[])})
            for method in METHODS:
                out={};arr={};validids={};badids={}
                for hyp in ['ph','line']:
                    good=[];bad=[];vals=[]
                    for ri,r in enumerate(records[hyp]):
                        fit={} if r is None else r.get('fits',{}).get(method+'_partial',{})
                        if fit.get('diagnostic_pass',False) and np.isfinite(fit['posterior']['x_ph'][2]):good.append(ri);vals.append(fit['posterior']['x_ph'][2])
                        else:bad.append(ri)
                    arr[hyp]=np.array(vals);validids[hyp]=good;badids[hyp]=bad
                    out[hyp]={'N_planned':planned,'N_valid':len(good),'failed_replicate_ids':bad,'b2_posterior_median_distribution':summary(vals) if vals else None}
                out['observed_b2']=obs[name][method]['b2'];complete=all(len(arr[hyp])==planned for hyp in arr);out['complete']=complete
                if complete:
                    out['p_ph']=tail(arr['ph'],out['observed_b2'],'le');out['p_line']=tail(arr['line'],out['observed_b2'],'ge');cut=float(np.quantile(arr['line'],.95));k=int(np.sum(arr['ph']>cut));out['power']={'critical_b2_Hline_p95':cut,'correct_Hph_count':k,'denominator':planned,'fraction':k/planned,'binomial_95_interval_conditional_on_estimated_cut':interval(k,planned),'note':'interval conditions on the Monte Carlo-estimated Hline cutoff; does not include cutoff estimation uncertainty'}
                else:
                    out['status']='未完成：固定2000个重复中存在未解决失败，不删除失败形成科学p值'
                    for hyp,side in [('ph','le'),('line','ge')]:
                        known=int(np.sum(arr[hyp]<=out['observed_b2'] if side=='le' else arr[hyp]>=out['observed_b2']));missing=len(badids[hyp]);out['p_'+hyp+'_failure_bounds']={'known_tail_count':known,'unknown_count':missing,'planned_denominator':planned,'possible_fraction_bounds':[known/planned,(known+missing)/planned]}
                bc=[];bc_bad=[]
                for ri,r in enumerate(records['line']):
                    fit={} if r is None else r.get('fits',{}).get(method+'_bc',{})
                    if fit.get('diagnostic_pass',False):bc.append([fit['posterior']['x_ph'][2],fit['posterior']['t'][2]])
                    else:bc_bad.append(ri)
                bco={'N_planned':planned,'N_valid':len(bc),'failed_replicate_ids':bc_bad,'complete':len(bc)==planned}
                if bc:
                    bc=np.array(bc);bco['b_distribution']=summary(bc[:,0]);bco['c_distribution']=summary(bc[:,1]);bco['observed_c']=obs[name][method]['c'];bco['observed_c_99_interval']=obs[name][method]['c_99_interval'];bco['observed_c_99_excludes_1']=bool(not(bco['observed_c_99_interval'][0]<=1<=bco['observed_c_99_interval'][1]))
                    if len(bc)==planned:
                        bco['c_lower_tail']=tail(bc[:,1],bco['observed_c'],'le');bco['c_upper_tail']=tail(bc[:,1],bco['observed_c'],'ge');bco['c_equal_tail_two_sided_fraction']=min(1.,2*min(bco['c_lower_tail']['fraction'],bco['c_upper_tail']['fraction']))
                out['Hline_bc']=bco;scenario_result['methods'][method]=out
            allout[name]['scenarios'][scenario]=scenario_result
    result={'scope':'F6-F7 fixed-design calibration; simulated statistics are actual posterior medians','N_per_hypothesis_scenario_sample':planned,'sample_order':SAMPLES,'scenarios':SCENARIOS,'methods':METHODS,'generating_parameter_source':'../F/generating_baselines.json','observed_statistic_source':'../F/observed_thresholds.json','observed_thresholds_sha256':hashlib.sha256((DATA/'photosphere/hierarchical/observed_thresholds.json').read_bytes()).hexdigest(),'data_sha256':manifest['data_sha256'],'baseline_sha256':manifest['baseline_sha256'],'results':allout,'numerical_audit':{'all_posteriors_valid':len(failures)==0,'failures':failures,'fit_method_counts':dict(methods),'number_of_attempts_distribution':dict(attempts),'replicate_algorithm_counts':dict(algorithms),'importance_ESS_min_p16_median_p84':np.quantile(ess,[0,.16,.5,.84]) if ess else None,'importance_max_weight_max':max(maxw) if maxw else None,'replicate_seconds_p16_p50_p84':np.quantile(durations,[.16,.5,.84])},'scientific_limits':['Fixed observed predictor medians used as simulation truths; no targeting/classification/selection simulator.','Gaussianized published percentile errors; missing Lph-T and cross-study covariances remain zero assumptions.','Main E has N=5, below the preregistered N>=20 gate regardless of simulated tails.','Curtis A own-Halpha is13-source matched subset; the full18-source BX-Lbol diagnostic is not used as a proxy for own Halpha.','A and B overlap and are not independent evidence.']}
    save(O/'F_simulation_results.json',result)
    print('complete',len(failures)==0,'fit count',sum(methods.values()),'failures',len(failures),flush=True)
if __name__=='__main__':main()
