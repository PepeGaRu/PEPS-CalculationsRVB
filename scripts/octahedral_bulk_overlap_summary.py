"""Recompute diagnostic, ratio-aware blocked statistics from saved runs."""
import json,numpy as np
from pathlib import Path
src=Path(__file__).resolve().parents[1]/'results/octahedral_bulk_overlap';summary=[]
for p in sorted(src.glob('mcmc_*.json')):
 r=json.load(open(p));data=r['records'];data=data[len(data)//4:];n=len(data);w=np.array([d.get('weight',1.) for d in data]);sign=np.array([d['sign'] for d in data]);den=w*sign;bs=np.array_split(np.arange(n),min(20,n//5));item={k:r[k] for k in ['L','steps','seed','initial','seconds','moves','negative_proposals']};item.update(file=p.name,bias=r.get('bias',1),burn_in_steps=r['steps']//4,measurements=n,mean_sign=float(np.sum(w*sign)/np.sum(w)),importance_weight_effective_samples=float(np.sum(w)**2/np.sum(w*w)),diagnostic_trimer_difference_fraction=float(np.mean([d['trimer_difference']>0 for d in data])),diagnostic_trimer_difference_weighted_fraction=float(np.sum(w*np.array([d['trimer_difference']>0 for d in data]))/np.sum(w)),correlations={})
 for key in sorted(k for k in data[0] if k.startswith('spin_')):
  a=w*np.array([d[key] for d in data]);est=float(np.sum(a)/np.sum(den));blockratios=np.array([(sum(a)-sum(a[b]))/(sum(den)-sum(den[b])) for b in bs]);err=float(np.sqrt((len(bs)-1)/len(bs)*sum((blockratios-blockratios.mean())**2)));nz=sum(d[key]!=0 for d in data);item['correlations'][key]=dict(value=est,block_standard_error=err,nonzero_measurements=nz,nonzero_seen=bool(nz))
 summary.append(item)
result=dict(status='Exploratory finite periodic boxes; no established bulk correlation length. Delete-one-block jackknife errors are statistical diagnostics only and do not include possible mixing bias.',observable='C(r)=<S_i dot S_j>, same physical edge axis, displacement parallel or perpendicular to that axis; randomized origins and axes.',sampling='Target |<bra|ket>| times bounded bias; all observables and signs reweighted by inverse bias. Each individual contraction exact.',runs=summary)
json.dump(result,open(src/'summary.json','w'),indent=2)
for r in summary:
 print(r['file'],'negative',r['negative_proposals'],'mismatch sampled/weighted',r['diagnostic_trimer_difference_fraction'],r['diagnostic_trimer_difference_weighted_fraction'],'weightESS',r['importance_weight_effective_samples'])
 for key,v in r['correlations'].items():print(' ',key,v)
