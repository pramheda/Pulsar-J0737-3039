"""Compare the same saved draws at two time spacings and coordinate choices."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'results'
out={}
for branch in ['published','corrected_coordinates']:
    a=np.load(p/f'mc_{branch}_n32768_step0.005.npz')
    b=np.load(p/f'mc_{branch}_n1024_step0.0025.npz')
    interior=(a['minimum_time_myr'][:1024]>0)&(a['minimum_time_myr'][:1024]<21.74)
    out[branch]={'n':1024,'event_disagreements':int(np.sum(a['events'][:1024]!=b['events'])),
                 'max_minimum_difference_pc':np.max(abs(a['minimum_pc'][:1024]-b['minimum_pc']),axis=0).tolist(),
                 'max_interior_minimum_difference_pc':float(np.max(abs(a['minimum_pc'][:1024]-b['minimum_pc'])[interior])),
                 'max_time_difference_myr':np.max(abs(a['minimum_time_myr'][:1024]-b['minimum_time_myr']),axis=0).tolist(),
                 'first1024_counts_step005':a['events'][:1024].sum(axis=0).tolist(),
                 'first1024_counts_step0025':b['events'].sum(axis=0).tolist()}
    assert out[branch]['event_disagreements']==0
a=np.load(p/'mc_published_n32768_step0.005.npz')
b=np.load(p/'mc_corrected_coordinates_n32768_step0.005.npz')
out['paired_inputs']={k:bool(np.array_equal(a[k],b[k])) for k in ['binary_distance_kpc','radial_heliocentric_kms','weight_uniform_heliocentric']}
assert all(out['paired_inputs'].values())
out['corrected_minima_global_pc']=b['minimum_pc'].min(axis=0).tolist()
out['corrected_minima_quantiles_pc']=np.quantile(b['minimum_pc'],[0,.01,.16,.5,.84,1],axis=0).tolist()
(p/'timestep_validation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
