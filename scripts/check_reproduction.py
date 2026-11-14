"""Short checks of the released scientific paths; no downloads or long production runs."""
import sys,json,hashlib,subprocess
from pathlib import Path
import numpy as np
from lrdwidth.paths import ROOT,DATA
sys.path.insert(0,str(ROOT/'analysis/injection'))
from common import posterior
fixture=json.loads((DATA/'measurements/recovery_check.json').read_text());s,g=posterior(fixture['rows'])
assert np.allclose(s['posterior']['intrinsic_scatter_dex'],fixture['sigma_quantiles'],rtol=0,atol=1e-10)
assert abs(s['posterior']['slope'][2]-fixture['slope_median'])<1e-10
print('PASS corrected recovery posterior from frozen per-source measurements')
for scenario in [1,3]:
 assert json.loads((DATA/f'widths/injection_recovery_s{scenario}.json').read_text())==json.loads((DATA/f'measurements/W2_final_summary_s{scenario}.json').read_text())
sys.path.insert(0,str(ROOT/'analysis/hierarchical'))
from numba_sampler import selfcheck
print('PASS hierarchical target:',selfcheck())
subprocess.run([sys.executable,str(ROOT/'analysis/hierarchical/run_simulations.py'),'--selfcheck-only'],check=True)
for script in ['make_profile_counts.py','make_measurement_summary.py','make_tables.py','make_numbers.py']:
 subprocess.run([sys.executable,str(ROOT/'scripts'/script)],cwd=ROOT,check=True)
for p,expected in json.loads((DATA/'measurements/paper_reference_hashes.json').read_text()).items():
 assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==expected,p
print('PASS 637 macro values and 366 manuscript/table definitions unchanged')
