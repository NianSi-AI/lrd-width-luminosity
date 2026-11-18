"""Build manuscript measurement summaries from BIC counts and recovery aggregates.
No posterior or empirical quantile is silently recomputed with a new convention.
The source map retains exact JSON pointers and physical units (dex, km/s, counts).
"""
import argparse,json
from pathlib import Path
from lrdwidth.paths import DATA,RESULTS

def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,default=DATA/'measurements');p.add_argument('--out',type=Path,default=RESULTS/'measurement_summary.json');a=p.parse_args()
 mapping=json.loads((DATA/'measurements/summary_map.json').read_text());numbers={};cache={}
 for k,m in mapping.items():
  if m['input'] not in cache:cache[m['input']]=json.loads((a.inputs/m['input']).read_text())
  value=cache[m['input']]
  for part in m['pointer'].strip('/').split('/'):
   part=part.replace('~1','/').replace('~0','~');value=value[int(part)] if isinstance(value,list) else value[part]
  assert isinstance(value,(int,float)),k
  numbers[k]={'value':value,'unit':m['unit']}
 assert len(numbers)==887,len(numbers)
 result={'quantile_convention':CONVENTION,'numbers':numbers};a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=1)+'\n')
 print('measurement summary fields:',len(numbers))
CONVENTION='Uniform midpoint cell mass; CDF at true cell boundaries with leading zero; empirical sample/replicate quantiles remain NumPy linear. Difference distributions retain original discrete convolution.'
if __name__=='__main__':main()
