"""Independent Cartesian MWPotential2014 integration of retained MC draws.

Analytic normalized forces are implemented here, without calling galpy inside
the ODE. Shared inputs, Astropy coordinate transforms, and potential definition
mean this validates numerical implementation, not independent physical evidence.
"""
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.special import gammainc
from astropy import units as u
from astropy.coordinates import SkyCoord,Galactocentric,CartesianDifferential
from astropy.time import Time

ROOT=Path(__file__).resolve().parents[1]
CAMPAIGN=ROOT.parent
RO=8.122;VO=220.
TIME_MYR=(RO*u.kpc/(VO*u.km/u.s)).to_value(u.Myr)
FRAME=Galactocentric(galcen_distance=RO*u.kpc,z_sun=.0208*u.kpc,
    galcen_v_sun=CartesianDifferential([12.9,245.6,7.78]*u.km/u.s))
CLUSTERS=json.loads((CAMPAIGN/'trajectory/data/key_clusters.json').read_text())['clusters']
from galpy.potential import MWPotential2014,evaluateRforces,evaluatezforces
from galpy.orbit import Orbit

def force(xyz):
    """Dimensionless Cartesian acceleration; circular-force fractions .05,.6,.35."""
    x,y,z=np.moveaxis(xyz,-1,0);r=np.linalg.norm(xyz,axis=-1)
    bulge=.05*gammainc(.6,(r/(1.9/8))**2)/gammainc(.6,(1/(1.9/8))**2)
    mass=lambda v:np.log1p(v/2)-v/(2+v)
    halo=.35*mass(r)/mass(1.)
    acc=-(bulge+halo)[...,None]*xyz/r[...,None]**3
    a=3/8;b=.28/8;B=np.sqrt(z*z+b*b)
    gm=.6*(1+(a+b)**2)**1.5
    den=(x*x+y*y+(a+B)**2)**1.5
    acc-=gm*np.stack([x,y,(a+B)*z/B],axis=-1)/den[...,None]
    return acc

def state(ra,dec,d,pmra,pmdec,rv):
    c=SkyCoord(ra=np.asarray(ra)*u.deg,dec=np.asarray(dec)*u.deg,distance=np.asarray(d)*u.kpc,
        pm_ra_cosdec=np.asarray(pmra)*u.mas/u.yr,pm_dec=np.asarray(pmdec)*u.mas/u.yr,
        radial_velocity=np.asarray(rv)*u.km/u.s).transform_to(FRAME)
    # Keep Astropy's x handedness in independent ODE. Axisymmetry permits this.
    return np.concatenate([c.cartesian.xyz.to_value(u.kpc).T/RO,
                          c.velocity.d_xyz.to_value(u.km/u.s).T/VO],axis=-1)

def select(data,number=32):
    """Deliberately include events, non-events, threshold-nearest and extreme RV."""
    evt=data['events'];p=data['minimum_pc'];rv=data['radial_heliocentric_kms'];out=[]
    def add(values):
        for v in values:
            if int(v) not in out:out.append(int(v))
    for j,c in enumerate(CLUSTERS):
        add(np.argsort(abs(p[:,j]-c['encounter_threshold_3r50_pc']))[:4])
        add(np.flatnonzero(evt[:,j])[:4]);add(np.argsort(p[:,j])[:2])
    add(np.argsort(rv)[:2]);add(np.argsort(rv)[-2:])
    add(np.flatnonzero(evt.all(axis=1))[:3]);add(np.flatnonzero(~evt.any(axis=1))[:3])
    add(np.linspace(0,len(evt)-1,number,dtype=int))
    return np.array(out[:number])

def draw_states(position,n,indices):
    z=np.random.default_rng(2026101909).standard_normal((n,12))
    components=np.searchsorted(np.cumsum([.463,.477,.060]),np.random.default_rng(2026101910).random(n))
    rv=np.array([43.1,107.6,406.1])[components]*z[:,3]
    ra,dec,pra,pde=144.448,-30.6278,-2.57,2.08
    if position!='published':
        c=SkyCoord(ra='07h37m51.248115s',dec='-30d39m40.70485s',distance=.74*u.kpc,
            pm_ra_cosdec=-2.567*u.mas/u.yr,pm_dec=2.082*u.mas/u.yr,
            radial_velocity=66.25*u.km/u.s,obstime=Time(55045.,format='mjd'))
        c=c.apply_space_motion(new_obstime=Time(2016.,format='jyear'))
        ra,dec=float(c.ra.deg),float(c.dec.deg)
        if position=='timing':pra,pde=-2.567,2.082
    bd=.74+.06*z[:,0];pmra=pra+.03*z[:,1];pmdec=pde+.04*z[:,2]
    rows=[state(ra,dec,bd[indices],pmra[indices],pmdec[indices],rv[indices])]
    for j,c in enumerate(CLUSTERS):
        k=4+4*j
        ds=(c['distance_q84_pc']-c['distance_q16_pc'])/2000
        rows.append(state(c['ra_deg'],c['dec_deg'],c['distance_pc']/1000+ds*z[indices,k],
            c['pm_ra_cosdec_masyr']+c['pm_ra_standard_error_masyr']*z[indices,k+1],
            c['pm_dec_masyr']+c['pm_dec_standard_error_masyr']*z[indices,k+2],
            c['heliocentric_radial_velocity_kms']+c['rv_standard_error_kms']*z[indices,k+3]))
    return np.stack(rows),bd,rv

def check(position,n=32768):
    path=CAMPAIGN/'trajectory/results'/f'mc_{position}_n{n}_step0.005.npz'
    data=np.load(path);indices=select(data);initial,bd,rv=draw_states(position,n,indices)
    assert np.array_equal(bd,data['binary_distance_kpc'])
    assert np.array_equal(rv,data['radial_heliocentric_kms'])
    shape=initial.shape
    def rhs(t,y):
        s=y.reshape(-1,6)
        return np.concatenate([s[:,3:],force(s[:,:3])],axis=1).ravel()
    tmax=max(c['age_Myr'] for c in CLUSTERS)
    sol=solve_ivp(rhs,[0,-tmax/TIME_MYR],initial.ravel(),method='DOP853',
                  rtol=2e-12,atol=2e-14,dense_output=True)
    assert sol.success
    def xyz(t):return sol.sol(-np.asarray(t)/TIME_MYR).reshape(shape+np.shape(t))[:,:,:3,...]*RO
    rows=[]
    for j,c in enumerate(CLUSTERS):
        horizon=np.floor(c['age_Myr']/.005)*.005;grid=np.linspace(0,horizon,2401)
        pos=xyz(grid);dist=np.linalg.norm(pos[0]-pos[j+1],axis=1)*1000
        for k,idx in enumerate(indices):
            p=int(np.argmin(dist[k]));bounds=(grid[max(p-1,0)],grid[min(p+1,len(grid)-1)])
            def objective(t):
                q=xyz(t);return float(np.linalg.norm(q[0,k]-q[j+1,k])*1000)
            opt=minimize_scalar(objective,bounds=bounds,method='bounded',options={'xatol':1e-10})
            options=[(0.,objective(0.)),(horizon,objective(horizon)),(float(opt.x),float(opt.fun))]
            mint,mind=min(options,key=lambda x:x[1]);flag=mind<c['encounter_threshold_3r50_pc']
            rows.append(dict(draw_index=int(idx),cluster=c['name'],ode_minimum_pc=mind,
                stored_galpy_minimum_pc=float(data['minimum_pc'][idx,j]),
                difference_pc=mind-float(data['minimum_pc'][idx,j]),
                ode_minimum_time_myr=mint,stored_galpy_minimum_time_myr=float(data['minimum_time_myr'][idx,j]),
                time_difference_myr=mint-float(data['minimum_time_myr'][idx,j]),
                ode_event=bool(flag),stored_galpy_event=bool(data['events'][idx,j]),
                threshold_pc=c['encounter_threshold_3r50_pc'],horizon_myr=horizon,cluster_median_age_myr=c['age_Myr']))
    # Direct path check for every selected binary and cluster using fresh galpy.
    s=initial.reshape(-1,6).copy();s[:,0]*=-1;s[:,3]*=-1
    x,y,z,vx,vy,vz=s.T;R=np.hypot(x,y)
    cyl=np.column_stack([R,(x*vx+y*vy)/R,(x*vy-y*vx)/R,z,vz,np.arctan2(y,x)])
    ts=-np.linspace(0,tmax,65)/TIME_MYR
    orbit=Orbit(cyl,ro=RO,vo=VO)
    orbit.integrate(ts,MWPotential2014,method='dop853_c',rtol=1e-11,atol=1e-13,progressbar=False,numcores=1)
    gx=np.stack([orbit.x(ts),orbit.y(ts),orbit.z(ts)],axis=1).reshape((3,len(indices),3,len(ts)))
    gx[:,:,0,:]*=-1
    path_diff=np.linalg.norm(gx-xyz(-ts*TIME_MYR),axis=2)*1000
    disagreements=sum(x['ode_event']!=x['stored_galpy_event'] for x in rows)
    out=dict(position=position,input_path=str(path.relative_to(CAMPAIGN)),
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),selection=indices.tolist(),
        n_selected=int(len(indices)),draw_reconstruction_exact=True,
        max_direct_path_difference_pc=float(path_diff.max()),
        max_minimum_difference_pc=max(abs(x['difference_pc']) for x in rows),
        max_minimum_time_difference_myr=max(abs(x['time_difference_myr']) for x in rows),
        event_disagreements=disagreements,rows=rows,ode_function_evaluations=sol.nfev)
    assert out['max_direct_path_difference_pc']<1e-4
    assert out['max_minimum_difference_pc']<.01
    assert disagreements==0
    return out

if __name__=='__main__':
    start=time.monotonic()
    probes=np.random.default_rng(314159).normal(size=(100,3));probes[:,0]+=1.3
    a=force(probes);rr=np.hypot(probes[:,0],probes[:,1]);zz=probes[:,2]
    radial=evaluateRforces(MWPotential2014,rr,zz,use_physical=False)
    vertical=evaluatezforces(MWPotential2014,rr,zz,use_physical=False)
    expected=np.column_stack([radial*probes[:,0]/rr,radial*probes[:,1]/rr,vertical])
    ferr=float(abs(a-expected).max());assert ferr<1e-12
    result=dict(scope='Numerical implementation validation of selected controlled Monte Carlo trajectories; not original author-code execution or independent observations.',
        analytic_force_max_absolute_difference=ferr,
        potential_definition='Bovy MWPotential2014: bulge alpha1.8 cutoff1.9/8, MiyamotoNagai a3/8 b0.28/8, NFW scale2; circular-force fractions .05 .6 .35 at R1 z0.',
        time_unit_myr=TIME_MYR,branches=[check('published'),check('corrected_coordinates')],runtime_seconds=time.monotonic()-start,
        limitations=['Selection deliberately targets event and threshold boundaries; it is not an independently sampled event-rate estimate.',
        'Astropy frame transformation and physical potential definition are shared; independent force algebra and Cartesian integrator test implementation.',
        'Stored MC minima use 0.005Myr grid with quadratic interpolation and float32 serialization; independent minima use dense SciPy output and bounded minimization.',
        'Both minimizations stop at the last 0.005Myr sample below the cluster median age, matching the source MC allowed-time mask exactly.'])
    (ROOT/'results/ode_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='branches'},indent=2))
    print(json.dumps([{k:v for k,v in b.items() if k!='rows'} for b in result['branches']],indent=2))
