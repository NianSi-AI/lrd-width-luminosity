"""Rebuild BIC profile counts from the 44 author-measured per-source records.
SNR bins and strict BIC comparisons follow the original width analysis.
"""
import argparse,json
from pathlib import Path
from lrdwidth.paths import DATA
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=DATA/'measurements/W1_BIC_counts.json');a=p.parse_args()
kinds=['single_gaussian','lorentzian','exponential','double_gaussian'];table=[]
for path in sorted((DATA/'widths/source_fits').glob('*.json')):
 r=json.loads(path.read_text());bic={k:r[k]['BIC'] for k in kinds};best=min(bic,key=bic.get);sn=r['single_gaussian']['broad_SNR']
 table.append(dict(cls=r['single_gaussian']['cls'],SNR_bin='<10' if sn<10 else '10-20' if sn<=20 else '>20',best=best,delta=bic['single_gaussian']-bic[best],single_range=max(bic[k] for k in kinds[:3])-min(bic[k] for k in kinds[:3])))
assert len(table)==44
counts=[]
for cl in ['LRD','non-LRD']:
 for sbin in ['all','<10','10-20','>20']:
  rows=[r for r in table if r['cls']==cl and (sbin=='all' or r['SNR_bin']==sbin)]
  counts.append(dict(cls=cl,SNR_bin=sbin,N=len(rows),best_counts={k:sum(r['best']==k for r in rows) for k in kinds},Gaussian_deltaBIC_gt6=sum(r['delta']>6 for r in rows),Gaussian_deltaBIC_gt10=sum(r['delta']>10 for r in rows),three_single_within6=sum(r['single_range']<6 for r in rows)))
assert sum(r['N'] for r in counts if r['SNR_bin']=='all')==44
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(counts,indent=2)+'\n');print('profile counts:44 sources')
