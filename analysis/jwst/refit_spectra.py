"""Uniform JWST batch fits: flux erg/s/cm^2, intrinsic widths km/s.
Uses the original forward Gaussian-LSF pixel integration, model comparison and QC.
Spectra remain external: see third_party/external_data.tsv. No automatic download.
"""
import argparse,json
from pathlib import Path
from lrdwidth import fit_batch as fit
from lrdwidth.paths import DATA

def main():
 p=argparse.ArgumentParser();p.add_argument('--spectra',type=Path,required=True,help='directory containing the manifest input_path tree');p.add_argument('--out',type=Path,default=Path('work/jwst'));p.add_argument('--limit',type=int);p.add_argument('--line',choices=fit.LINES);a=p.parse_args()
 fit.c.self_checks()
 targets=json.loads((DATA/'jwst/targets.json').read_text());assert len(targets)==204
 targets=targets[:a.limit] if a.limit else targets
 missing=[a.spectra/t['input_path'] for t in targets if not (a.spectra/t['input_path']).is_file()]
 if missing:raise SystemExit(f'{len(missing)} missing spectra; first: {missing[0]}. See third_party/external_data.tsv')
 fit.O=a.spectra;fit.W=a.out;a.out.mkdir(parents=True,exist_ok=True)
 if a.line:fit.LINES=[a.line]
 for t in targets:print(fit.one_file(t),flush=True)
if __name__=='__main__':main()
