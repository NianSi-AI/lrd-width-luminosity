"""List or download external spectra at their frozen byte versions. No pixels in git."""
import argparse,csv,json,hashlib,urllib.request
from pathlib import Path
from lrdwidth.paths import DATA,ROOT
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path('work/jwst_spectra'));p.add_argument('--download',action='store_true');p.add_argument('--frozen-only',action='store_true');p.add_argument('--limit',type=int)
p.add_argument('--qso1',action='store_true',help='also the A2744-QSO1 spectrum used by analysis/jwst/ (not one of the 204 targets)');a=p.parse_args()
rows=json.loads((DATA/'jwst/targets.json').read_text());assert len(rows)==204
if a.frozen_only:
 ids={r['source']['file_id'] for r in json.loads((DATA/'jwst/injection_parameters.json').read_text())};rows=[r for r in rows if r['file_id'] in ids];assert len(rows)==44
if a.limit:rows=rows[:a.limit]
if a.qso1:
 q=json.loads((DATA/'qso1/qso1_fit.json').read_text());name=Path(q['spectrum']).name
 ext=next(r for r in csv.DictReader((ROOT/'third_party/external_data.tsv').open(),delimiter='\t') if r['file']==name);assert ext['sha256']==q['sha256']
 rows=rows+[dict(file=name,input_path='spectra/'+name,download_url=ext['download_url'],sha256=q['sha256'])]
for r in rows:
 dest=a.out/r['input_path'];assert dest.resolve().is_relative_to(a.out.resolve())
 if not a.download:print(r['download_url'],dest);continue
 if dest.exists():data=dest.read_bytes()
 else:
  with urllib.request.urlopen(r['download_url'],timeout=120) as response:data=response.read()
 assert hashlib.sha256(data).hexdigest()==r['sha256'],f"Spectrum version changed: {r['file']}"
 if not dest.exists():dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
 print('verified',r['file'])
