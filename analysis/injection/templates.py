"""Frozen Halpha templates loaded from external spectra.
Velocities km/s; y/err in 1e-18 CGS per 1000 km/s. Only author fit parameters
are redistributed. Published DJA pixels remain link-only.
"""
import json,os,hashlib
from pathlib import Path
import numpy as np
from lrdwidth.paths import DATA
from lrdwidth import linefit as core
TEMPLATES=json.loads((DATA/'jwst/injection_parameters.json').read_text())
FROZEN=[t['source'] for t in TEMPLATES]
assert len(FROZEN)==44 and sum(r['lrd_class_stage2']=='LRD' for r in FROZEN)==29
MU=float(np.mean([r['y_logFWHM_kms'] for r in FROZEN if r['lrd_class_stage2']=='LRD']))
CACHE={}
def template(i):
 if i not in CACHE:
  t=TEMPLATES[i];p=Path(os.environ.get('LRDWIDTH_SPECTRA','work/jwst_spectra'))/t['input_path']
  if not p.is_file():raise FileNotFoundError(f'{p}: set LRDWIDTH_SPECTRA; see third_party/external_data.tsv')
  assert hashlib.sha256(p.read_bytes()).hexdigest()==t['input_sha256'],'Spectrum version differs from frozen input'
  d=core.load_line(p,float(t['source']['z_fit']),'Halpha')
  assert hashlib.sha256(np.asarray(d['err'],dtype='<f8').tobytes()).hexdigest()==t['error_array_sha256'],'Noise preprocessing changed'
  CACHE[i]=(d,t['parameters'])
 d,parameters=CACHE[i]
 return dict(d),parameters.copy()
