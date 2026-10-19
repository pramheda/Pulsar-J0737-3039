"""Shared-trajectory importance-sampled test, not an exact publication replay."""
from galactic_trace import ROOT,RO,VO,ZO,FRAME,TIME_MYR,CLUSTERS,BIN,timing_binary,sky,integrate,coordinates
import argparse,json,time
import numpy as np
from scipy.stats import norm,t as student_t,beta as beta_distribution
from galpy.potential import MWPotential2014,vcirc
from astropy import units as u

ALPHA=np.array([.463,.477,.060]);SIGMA=np.array([43.1,107.6,406.1])

def peculiar_projection(ra,dec,distance,pmra,pmdec):
    zero=sky(ra,dec,distance,pmra,pmdec,np.zeros(len(distance))).transform_to(FRAME)
    one=sky(ra,dec,distance,pmra,pmdec,np.ones(len(distance))).transform_to(FRAME)
    xyz=zero.cartesian.xyz.to_value(u.kpc).T*np.array([-1,1,1])
    vel=zero.velocity.d_xyz.to_value(u.km/u.s).T*np.array([-1,1,1])
    sight=(one.velocity.d_xyz-zero.velocity.d_xyz).to_value(u.km/u.s).T*np.array([-1,1,1])
    assert np.max(abs(np.linalg.norm(sight,axis=1)-1))<1e-12
    R=np.hypot(xyz[:,0],xyz[:,1]);circular=vcirc(MWPotential2014,R/RO,use_physical=False)*VO
    mean=np.stack([-circular*xyz[:,1]/R,circular*xyz[:,0]/R,np.zeros(len(R))],axis=-1)
    peculiar=vel-mean;radial=np.einsum('ij,ij->i',peculiar,sight)
    tangent=peculiar-radial[:,None]*sight
    return -radial,np.linalg.norm(tangent,axis=1)

def importance_priors(rv,center,transverse):
    proposal=norm.pdf(rv[:,None],scale=SIGMA)@ALPHA
    marginal=norm.pdf(rv[:,None]-center[:,None],scale=SIGMA)@ALPHA
    beta=ALPHA/SIGMA**2*np.exp(-.5*(transverse[:,None]/SIGMA)**2);beta/=beta.sum(axis=1)[:,None]
    conditional=np.sum(beta*norm.pdf(rv[:,None]-center[:,None],scale=SIGMA),axis=1)
    a=np.sqrt(51.2**2+transverse**2)
    cauchy=student_t.pdf(rv-center,df=3,scale=a/np.sqrt(3))
    uniform=(abs(rv)<=200)/400
    return dict(gaussian_heliocentric=np.ones(len(rv)),uniform_heliocentric=uniform/proposal,
                gaussian_peculiar_marginal=marginal/proposal,
                gaussian_peculiar_conditional=conditional/proposal,
                cauchy_peculiar_conditional=cauchy/proposal)

def summarize(events,weights,direct_prior='gaussian_heliocentric'):
    names=list(weights);rows=[]
    f=np.column_stack([events[:,0],events[:,1],events.all(axis=1),events.any(axis=1)]).astype(float)
    for name,w in weights.items():
        mass=w.mean();p=np.sum(w[:,None]*f,axis=0)/w.sum()
        se=np.sqrt(np.var(w[:,None]*(f-p),axis=0,ddof=1)/len(w))/mass
        baseline_limits=None
        if name==direct_prior:
            baseline_limits={}
            for column,label in enumerate(['Theia_58','OC_0450','intersection','true_union']):
                count=int(f[:,column].sum());number=len(w)
                baseline_limits[label]=dict(count=count,one_sided_95_lower=0. if count==0 else float(beta_distribution.ppf(.05,count,number-count+1)),
                                             one_sided_95_upper=1. if count==number else float(beta_distribution.ppf(.95,count+1,number-count)))
        rows.append(dict(prior=name,weight_mean=float(mass),weight_se=float(w.std(ddof=1)/np.sqrt(len(w))),
                         effective_sample_size=float(w.sum()**2/(w@w)),
                         probabilities=dict(zip(['Theia_58','OC_0450','intersection','true_union'],p.tolist())),
                         monte_carlo_se=dict(zip(['Theia_58','OC_0450','intersection','true_union'],se.tolist())),
                         product_union=float(1-(1-p[0])*(1-p[1])),event_covariance=float(p[2]-p[0]*p[1]),
                         baseline_binomial_one_sided_95=baseline_limits))
    return rows

def run(n=4096,position='timing',step=.005,chunk=64,seed=2026101909,radial_prior='gaussian_heliocentric'):
    start=time.monotonic();b=BIN if position=='published' else timing_binary()
    if position=='corrected_coordinates':
        corrected=timing_binary()
        b=dict(BIN,ra_deg=corrected['ra_deg'],dec_deg=corrected['dec_deg'],
               position_source='Kramer2021 position propagated to2016; ALL other binary inputs and draws IDENTICAL to published-coordinate branch')
    rng=np.random.default_rng(seed);z=rng.standard_normal((n,12))
    rngc=np.random.default_rng(seed+1);component_uniform=rngc.random(n)
    components=np.searchsorted(np.cumsum(ALPHA),component_uniform)
    bd=b['distance_kpc']+b['distance_sigma_kpc']*z[:,0]
    bpra=b['pmra']+b['pmra_sigma']*z[:,1];bpdec=b['pmdec']+b['pmdec_sigma']*z[:,2]
    assert np.all(bd>0)
    center,vt=peculiar_projection(b['ra_deg'],b['dec_deg'],bd,bpra,bpdec)
    if radial_prior=='gaussian_heliocentric':
        rv=SIGMA[components]*z[:,3]
        weights=importance_priors(rv,center,vt)
    else:
        if radial_prior=='uniform_heliocentric':
            rv=400*component_uniform-200
        elif radial_prior=='gaussian_peculiar_marginal':
            rv=center+SIGMA[components]*z[:,3]
        elif radial_prior=='gaussian_peculiar_conditional':
            probabilities=ALPHA/SIGMA**2*np.exp(-.5*(vt[:,None]/SIGMA)**2)
            probabilities/=probabilities.sum(axis=1)[:,None]
            components=np.sum(component_uniform[:,None]>np.cumsum(probabilities,axis=1),axis=1)
            rv=center+SIGMA[components]*z[:,3]
        elif radial_prior=='cauchy_peculiar_conditional':
            rv=center+np.sqrt(51.2**2+vt**2)/np.sqrt(3)*np.random.default_rng(seed+2).standard_t(3,n)
        else:
            raise ValueError(radial_prior)
        weights={radial_prior:np.ones(n)}
    cluster_data=[]
    for j,c in enumerate(CLUSTERS):
        offset=4+4*j
        dsigma=(c['distance_q84_pc']-c['distance_q16_pc'])/2/1000
        cluster_data.append((c['distance_pc']/1000+dsigma*z[:,offset],
                             c['pm_ra_cosdec_masyr']+c['pm_ra_standard_error_masyr']*z[:,offset+1],
                             c['pm_dec_masyr']+c['pm_dec_standard_error_masyr']*z[:,offset+2],
                             c['heliocentric_radial_velocity_kms']+c['rv_standard_error_kms']*z[:,offset+3]))
    maxage=max(c['age_Myr'] for c in CLUSTERS)
    lookback=np.arange(0,maxage+step,step);tt=-lookback;natural=tt/TIME_MYR
    events=np.zeros((n,2),bool);minimum=np.zeros((n,2));minimum_time=np.zeros((n,2))
    grid_flags=np.zeros((n,2),bool)
    for first in range(0,n,chunk):
        last=min(n,first+chunk);m=last-first;sl=slice(first,last)
        coord=sky(np.repeat([b['ra_deg']]+[c['ra_deg'] for c in CLUSTERS],m),
                  np.repeat([b['dec_deg']]+[c['dec_deg'] for c in CLUSTERS],m),
                  np.concatenate([bd[sl]]+[d[0][sl] for d in cluster_data]),
                  np.concatenate([bpra[sl]]+[d[1][sl] for d in cluster_data]),
                  np.concatenate([bpdec[sl]]+[d[2][sl] for d in cluster_data]),
                  np.concatenate([rv[sl]]+[d[3][sl] for d in cluster_data]))
        o=integrate(coord,tt);xyz=coordinates(o,natural).reshape(3,m,len(tt),3)
        for j,c in enumerate(CLUSTERS):
            allowed=lookback<=c['age_Myr'];d2=np.sum((xyz[0,:,allowed,:]-xyz[j+1,:,allowed,:])**2,axis=-1)*1e6
            # Boolean advanced indexing puts time first; restore row=realization.
            if d2.shape!=(m,int(allowed.sum())):d2=d2.T
            ix=np.arange(m);k=np.argmin(d2,axis=1);y=d2[ix,k]
            l=d2[ix,np.maximum(k-1,0)];r=d2[ix,np.minimum(k+1,d2.shape[1]-1)]
            den=l+r-2*y;interior=(k>0)&(k<d2.shape[1]-1)&(den>0)
            shift=np.zeros(m);reduction=np.zeros(m)
            shift[interior]=np.clip((l[interior]-r[interior])/(2*den[interior]),-1,1)
            reduction[interior]=(r[interior]-l[interior])**2/(8*den[interior])
            d=np.sqrt(np.maximum(0,y-reduction));mint=(k+shift)*step
            minimum[sl,j]=d;minimum_time[sl,j]=mint
            grid_flags[sl,j]=y<c['encounter_threshold_3r50_pc']**2
            events[sl,j]=d<c['encounter_threshold_3r50_pc']
        if first==0 or last==n or last%1024==0:
            print('progress',position,last,n,'seconds',round(time.monotonic()-start,1),flush=True)
    summary=summarize(events,weights,radial_prior)
    tag=f'mc_{position}_n{n}_step{step:g}'
    if radial_prior!='gaussian_heliocentric':tag+='_'+radial_prior
    np.savez_compressed(ROOT/'results'/f'{tag}.npz',events=events,minimum_pc=minimum.astype(np.float32),
                        minimum_time_myr=minimum_time.astype(np.float32),binary_distance_kpc=bd,
                        radial_heliocentric_kms=rv,peculiar_radial_center_kms=center,transverse_peculiar_kms=vt,
                        **{'weight_'+key:value for key,value in weights.items()})
    out=dict(scope='Independent conditional two-cluster calculation; no original execution files recovered.',
             position_choice=position,binary=b,n=n,step_myr=step,seed=seed,chunk=chunk,
             direct_radial_prior=radial_prior,
             proposals='Unconditional zero-centered heliocentric three-Gaussian mixture; alternatives importance-reweighted.' if radial_prior=='gaussian_heliocentric' else 'Direct independent samples from named radial prior; all non-radial normal draws shared with baseline.',
             covariance_model='Fixed RA/Dec; independent Gaussian distance/pmRA/pmDec and cluster-mean RV errors; symmetric cluster distance error; fixed median ages.',
             peculiar_frame='Astropy Galactocentric with given solar velocity; local circular velocity from MWPotential2014 at each sampled position, no vertical mean.',
             effective_observation_epoch=2016.,summary=summary,
             half_sample_summary=summarize(events[:n//2],{k:v[:n//2] for k,v in weights.items()},radial_prior),
             grid_to_quadratic_event_changes=int(np.count_nonzero(events!=grid_flags)),
             raw_event_counts=events.sum(axis=0).tolist(),raw_intersection_count=int(events.all(axis=1).sum()),
             prior_projection=dict(center_quantiles=np.quantile(center,[.16,.5,.84]).tolist(),vt_quantiles=np.quantile(vt,[.16,.5,.84]).tolist()),
             runtime_seconds=time.monotonic()-start,
             limitations=['Two clusters only; neither all-catalog union nor posterior probability of origin.',
                          'Original coordinate, radial-frame and covariance implementation unrecovered.',
                          'No Gaia cluster-center systematic covariance or birth/ejection population selection.',
                          'No conclusion that changing a prior identifies the physical truth.',
                          'A zero-event weighted estimate is not a zero true probability; exact binomial limits only apply to the unweighted proposal prior.',
                          'Cauchy importance weights have finite-sample ESS only; Gaussian proposal does not certify far-tail convergence of polynomial Cauchy tails.'])
    (ROOT/'results'/f'{tag}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(tag=tag,summary=summary,runtime=out['runtime_seconds']),indent=2),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=4096)
    p.add_argument('--position',choices=['published','timing','corrected_coordinates'],default='timing')
    p.add_argument('--step',type=float,default=.005);p.add_argument('--chunk',type=int,default=64)
    p.add_argument('--prior',choices=['gaussian_heliocentric','uniform_heliocentric','gaussian_peculiar_marginal','gaussian_peculiar_conditional','cauchy_peculiar_conditional'],default='gaussian_heliocentric')
    a=p.parse_args();run(a.n,a.position,a.step,a.chunk,radial_prior=a.prior)
