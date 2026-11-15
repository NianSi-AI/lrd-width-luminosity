"""Uniform JWST batch fits: flux erg/s/cm^2, intrinsic widths km/s.
Uses the original forward Gaussian-LSF pixel integration, model comparison and QC.
Spectra remain external: see third_party/external_data.tsv. No automatic download.

Default settings reproduce the paper's fits (LSF x1, x0.7, x1.3; +-5000 km/s; 200 nonparametric bootstrap draws).
Sensitivity variants (write each to its own --out directory; every record stores its fit_settings):
  --frozen-only --line Halpha --lsf-scales 0.5                 halved LSF width (data/widths/source_fits_lsf_half/)
  --frozen-only --line Halpha --window-kms 6800                +-6800 km/s window, [S II] outside (win6800_noSII)
  --frozen-only --line Halpha --window-kms 8000 --mask-sii     +-8000 km/s window, [S II] masked (win8000_maskSII)
"""
import argparse,json
from pathlib import Path
from lrdwidth import fit_batch as fit
from lrdwidth.paths import DATA

def main():
 p=argparse.ArgumentParser();p.add_argument('--spectra',type=Path,required=True,help='directory containing the manifest input_path tree');p.add_argument('--out',type=Path,default=Path('work/jwst'));p.add_argument('--limit',type=int);p.add_argument('--line',choices=fit.LINES)
 p.add_argument('--frozen-only',action='store_true',help='only the 44 frozen-sample spectra (data/jwst/injection_parameters.json)')
 p.add_argument('--lsf-scales',type=float,nargs='+',default=fit.LSF_SCALES,help='LSF FWHM multipliers (default: %(default)s)')
 p.add_argument('--window-kms',type=float,help='fit half-width in km/s (default 5000; 4500 for Hbeta)')
 p.add_argument('--mask-sii',action='store_true',help='exclude [S II] 6716,6731 from the Halpha window')
 p.add_argument('--mask-halfwidth-kms',type=float,default=500.,help='half-width of the [S II] masks (default: %(default)s km/s; not recorded for the stored +-8000 km/s run)')
 p.add_argument('--nboot-np',type=int,default=fit.NBOOT_NP,help='nonparametric bootstrap draws (default: %(default)s)')
 a=p.parse_args()
 if a.mask_sii and a.line!='Halpha':p.error('--mask-sii applies to --line Halpha only')
 fit.c.self_checks()
 targets=json.loads((DATA/'jwst/targets.json').read_text());assert len(targets)==204
 if a.frozen_only:
  ids={r['source']['file_id'] for r in json.loads((DATA/'jwst/injection_parameters.json').read_text())};targets=[t for t in targets if t['file_id'] in ids];assert len(targets)==44
 targets=targets[:a.limit] if a.limit else targets
 missing=[a.spectra/t['input_path'] for t in targets if not (a.spectra/t['input_path']).is_file()]
 if missing:raise SystemExit(f'{len(missing)} missing spectra; first: {missing[0]}. See third_party/external_data.tsv')
 fit.O=a.spectra;fit.W=a.out;a.out.mkdir(parents=True,exist_ok=True)
 if a.line:fit.LINES=[a.line]
 fit.LSF_SCALES=list(a.lsf_scales);fit.WINDOW_KMS=a.window_kms;fit.NBOOT_NP=a.nboot_np
 if a.mask_sii:fit.MASK_REST_UM=fit.SII_REST_UM;fit.MASK_HALFWIDTH_KMS=a.mask_halfwidth_kms
 print('fit settings:',fit.settings(),flush=True)
 for t in targets:print(fit.one_file(t),flush=True)
if __name__=='__main__':main()
