"""Same robust posterior integration as v2; compiled exact-target MCMC fallback.

No scientific prior or likelihood changes. Existing installed numba compiles
the verified two-predictor log density and MH loop. Every algorithm attempt
and failed diagnostic is returned for the simulation audit trail.
"""
import time
from importance_v2 import integrate
from numba_sampler import sample_compiled
from diagnostics_v2 import all_diagnostics

def fit(model,seed=20261112,draws=8192,min_ess=400,max_weight=.02):
    assert draws>=4 and draws%2==0
    attempts=[];start=None;t=time.monotonic()
    for j,n in enumerate(sorted(set([draws,max(draws,32768),max(draws,131072)]))):
        try:
            result,z,w,start=integrate(model,seed=seed+7919*j,n=n,start=start)
            result['method']='exact_posterior_importance_compiled_target'
            passed=result['ESS']>=min_ess and result['max_weight']<=max_weight
            attempts.append({k:result[k] for k in ['method','draws','ESS','max_weight','MAP_success']})
            if passed:
                result.update(diagnostic_pass=True,attempts=attempts,seconds=time.monotonic()-t,algorithm_version='v3_compiled_same_target')
                return result
        except Exception as err:
            attempts.append({'method':'importance','draws':n,'error':type(err).__name__+': '+str(err)});start=None
    result=None
    for j,(steps,burn,thin) in enumerate([(40000,10000,5),(120000,30000,10),(240000,60000,10)]):
        try:
            result,ch,start=sample_compiled(model,seed=seed+100000+7919*j,steps=steps,burn=burn,thin=thin,start=start)
            all_diag=result['sampling']['diagnostics_all_parameters_rank_folded']
            passed=all_diag['max_Rhat']<1.01 and all_diag['min_bulk_ESS']>=400 and result['converged']
            attempts.append({'method':result['method'],'steps':steps,'converged_all':bool(passed),'all_diagnostics':all_diag})
            result['converged']=bool(passed)
            if passed:break
        except Exception as err:
            attempts.append({'method':'MCMC','steps':steps,'error':type(err).__name__+': '+str(err)});start=None
    if result is None:result={'N':model.n,'method':'failed','posterior':{},'converged':False}
    result.update(diagnostic_pass=bool(result['converged']),attempts=attempts,seconds=time.monotonic()-t,algorithm_version='v3_compiled_same_target')
    return result
