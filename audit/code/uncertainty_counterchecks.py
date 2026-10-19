"""Explicit adversarial uncertainty alternatives; never inferred as author inputs."""
from pathlib import Path
import sys,json,time
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];CAMPAIGN=ROOT.parent
sys.path.insert(0,str(CAMPAIGN/'trajectory/code'))
from galactic_trace import sky,integrate,coordinates,TIME_MYR,BIN,timing_binary,CLUSTERS
import numpy as np
from scipy.stats import beta

def run(name,n=4096,position='correct'):
    start=time.monotonic();z=np.random.default_rng(2026101909).standard_normal((n,12))
    components=np.searchsorted(np.cumsum([.463,.477,.060]),np.random.default_rng(2026101910).random(n))
    rv=np.array([43.1,107.6,406.1])[components]*z[:,3]
    b=dict(BIN) if position=='printed' else dict(BIN,ra_deg=timing_binary()['ra_deg'],dec_deg=timing_binary()['dec_deg'])
    bd=.74+.06*z[:,0];pra=-2.57+.03*z[:,1];pde=2.08+.04*z[:,2]
    clusters=[];ages=[];input_errors=[]
    for j,c in enumerate(CLUSTERS):
        raw=c['raw_record'];ps=[c['pm_ra_standard_error_masyr'],c['pm_dec_standard_error_masyr']]
        rs=c['rv_standard_error_kms']
        if name in ['rv_member_scatter','rv_pm_member_scatter']:rs=c['rv_member_scatter_kms']
        if name in ['pm_member_scatter','rv_pm_member_scatter']:ps=[float(raw['s_pmRA']),float(raw['s_pmDE'])]
        age=c['age_q84_Myr'] if name=='age84' else c['age_Myr'];ages.append(age)
        ds=(c['distance_q84_pc']-c['distance_q16_pc'])/2000;o=4+4*j
        clusters.append((c['distance_pc']/1000+ds*z[:,o],c['pm_ra_cosdec_masyr']+ps[0]*z[:,o+1],
                         c['pm_dec_masyr']+ps[1]*z[:,o+2],c['heliocentric_radial_velocity_kms']+rs*z[:,o+3]))
        input_errors.append(dict(name=c['name'],distance_sigma_kpc=ds,pmra_sigma=ps[0],pmdec_sigma=ps[1],rv_sigma=rs,age_horizon_myr=age))
    step=.005;times=np.arange(0,max(ages)+step,step);events=np.zeros((n,2),bool);mins=np.zeros((n,2));min_times=np.zeros((n,2))
    for first in range(0,n,64):
        last=min(n,first+64);m=last-first;sl=slice(first,last)
        coord=sky(np.repeat([b['ra_deg']]+[c['ra_deg'] for c in CLUSTERS],m),
                  np.repeat([b['dec_deg']]+[c['dec_deg'] for c in CLUSTERS],m),
                  np.concatenate([bd[sl]]+[c[0][sl] for c in clusters]),np.concatenate([pra[sl]]+[c[1][sl] for c in clusters]),
                  np.concatenate([pde[sl]]+[c[2][sl] for c in clusters]),np.concatenate([rv[sl]]+[c[3][sl] for c in clusters]))
        xyz=coordinates(integrate(coord,-times),-times/TIME_MYR).reshape(3,m,len(times),3)
        for j,c in enumerate(CLUSTERS):
            valid=times<=ages[j];d2=np.sum((xyz[0]-xyz[j+1])**2,axis=-1)[:,valid]*1e6
            k=d2.argmin(axis=1);ii=np.arange(m);d=np.sqrt(d2[ii,k])
            events[sl,j]=d<c['encounter_threshold_3r50_pc'];mins[sl,j]=d;min_times[sl,j]=times[k]
    out=dict(name=name,n=n,position=position,seed=2026101909,
             radial_prior='Zero-centered heliocentric 3-Gaussian mixture',cluster_error_inputs=input_errors,
             counts=events.sum(axis=0).tolist(),probabilities=events.mean(axis=0).tolist(),union_count=int(events.any(axis=1).sum()),
             minimum_pc=mins.min(axis=0).tolist(),binomial_upper95=[float(beta.ppf(.95,int(k)+1,n-int(k))) for k in events.sum(axis=0)],
             runtime_seconds=time.monotonic()-start,
             limitations=['Member dispersion is not the uncertainty of the cluster mean; these branches test possible implementation ambiguity, not preferred measurement errors.',
                          'This is a bounded diagnostic, not a full posterior under cluster evolution or source covariance.',
                          'Grid minimum used; baseline convergence tests show grid effect negligible at5kyr for these thresholds.'])
    (ROOT/'results'/f'uncertainty_{name}.json').write_text(json.dumps(out,indent=2)+'\n')
    np.savez_compressed(ROOT/'results'/f'uncertainty_{name}.npz',events=events,minimum_pc=mins.astype('f4'),minimum_time_myr=min_times.astype('f4'))
    print(json.dumps(out),flush=True)
    return out

if __name__=='__main__':
    out=[run(name,4096 if name!='age84' else 1024) for name in ['rv_member_scatter','pm_member_scatter','rv_pm_member_scatter','age84']]
    (ROOT/'results/uncertainty_counterchecks.json').write_text(json.dumps(out,indent=2)+'\n')
