"""Invert vector-figure curves under a declared potential and endpoint-time assumption.

The fit is a forensic numerical consistency test, not an observation or recovered code.
"""
from pathlib import Path
import sys,json
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];CAMPAIGN=ROOT.parent
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gammainc
from scipy.optimize import least_squares,minimize_scalar
from astropy import units as u
from astropy.coordinates import SkyCoord,Galactocentric,CartesianRepresentation,CartesianDifferential
RO=8.122;VO=220.;TIME=(RO*u.kpc/(VO*u.km/u.s)).to_value(u.Myr)
K=(1*u.kpc/u.Myr).to_value(u.km/u.s)

def force(xyz):
    x,y,z=xyz/RO;r=np.linalg.norm(xyz)/RO
    bulge=.05*gammainc(.6,(r/(1.9/8))**2)/gammainc(.6,(1/(1.9/8))**2)
    mass=lambda v:np.log1p(v/2)-v/(2+v)
    halo=.35*mass(r)/mass(1.)
    acc=-(bulge+halo)*np.array([x,y,z])/r**3
    a=3/8;b=.28/8;B=np.sqrt(z*z+b*b);gm=.6*(1+(a+b)**2)**1.5
    acc-=gm*np.array([x,y,(a+B)*z/B])/(x*x+y*y+(a+B)**2)**1.5
    return acc*RO/TIME**2
def rhs(t,y):return np.concatenate([y[3:],force(y[:3])])
def integrate(state,times):
    s=solve_ivp(rhs,(0,times[-1]),state,t_eval=times,method='DOP853',rtol=3e-12,atol=3e-14)
    assert s.success;return s.y.T
def icrs(state,solar,zo=.0208):
    frame=Galactocentric(galcen_distance=8.122*u.kpc,z_sun=zo*u.kpc,galcen_v_sun=CartesianDifferential(solar*u.km/u.s))
    position=state[:3]*[-1,1,1];vel=state[3:]*[-1,1,1]*K
    c=SkyCoord(CartesianRepresentation(position*u.kpc,differentials=CartesianDifferential(vel*u.km/u.s)),frame=frame).icrs
    return dict(ra_deg=float(c.ra.deg),dec_deg=float(c.dec.deg),distance_kpc=float(c.distance.to_value(u.kpc)),
                pmra=float(c.pm_ra_cosdec.to_value(u.mas/u.yr)),pmdec=float(c.pm_dec.to_value(u.mas/u.yr)),rv=float(c.radial_velocity.to_value(u.km/u.s)))

def run():
    d=json.loads((ROOT/'results/vector_figure.json').read_text());rows=[]
    for name,path in d['paths'].items():
        xy=np.array(path['xy']);xz=np.array(path['xz']);points=np.column_stack([(xy[:,0]+xz[:,0])/2,xy[:,1],xz[:,1]])
        initial=np.r_[points[0],-(points[1]-points[0])/.05,.05]
        def fun(params):return ((integrate(params[:6],-np.arange(len(points))*params[6])[:,:3]-points)*1000).ravel()
        fit=least_squares(fun,initial,bounds=(np.r_[[-np.inf]*6,.02],np.r_[[np.inf]*6,.08]),x_scale=np.r_[[1]*3,[.1]*3,.05],xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=200)
        nominal_last=7.43 if name=='Theia_58' else 8.04
        lookback_first=nominal_last-(len(points)-1)*fit.x[6]
        state=integrate(fit.x[:6],np.array([0.,lookback_first]))[-1]
        frames={label:icrs(state,np.array(solar),zo) for label,solar,zo in [('stated',[12.9,245.6,7.78],.0208),('schoenrich',[11.1,232.24,7.25],.0208),('schoenrich_z25',[11.1,232.24,7.25],.025)]}
        sensitivity=[];phase_free=None
        if name=='binary':
            for first in [4.,5.,6.,6.59,7.,8.,9.,10.]:
                now=integrate(fit.x[:6],np.array([0.,first]))[-1]
                sensitivity.append(dict(assumed_first_lookback_myr=first,**icrs(now,np.array([12.9,245.6,7.78]))))
            # Remove the chosen endpoint epoch: scan every nonnegative first-vertex
            # lookback within the oldest of the two median cluster ages, rounded up.
            target=SkyCoord(ra=114.46352849386929*u.deg,dec=-30.661303192105727*u.deg)
            phases=np.linspace(0,24,241)
            states=integrate(fit.x[:6],phases)
            def separation(state):
                sky=icrs(state,np.array([12.9,245.6,7.78]))
                return float(SkyCoord(ra=sky['ra_deg']*u.deg,dec=sky['dec_deg']*u.deg).separation(target).deg)
            angles=np.array([separation(s) for s in states]);k=int(angles.argmin())
            def fun_phase(t):return separation(integrate(fit.x[:6],np.array([0.,t]))[-1]) if t>0 else separation(fit.x[:6])
            opt=minimize_scalar(fun_phase,bounds=(phases[max(0,k-1)],phases[min(len(phases)-1,k+1)]),method='bounded')
            phase_free=dict(first_lookback_range_myr=[0,24],minimum_angular_separation_from_primary_deg=float(opt.fun),
                            minimizing_first_lookback_myr=float(opt.x),minimum_galpy_z_kpc=float(states[:,2].min()),
                            scope='Fixed potential and stated spatial frame; frees only temporal placement of the reconstructed red curve. Not a global optimization over plotting conventions.')
        rows.append(dict(name=name,n_points=len(points),inferred_time_spacing_myr=float(fit.x[6]),
                         rms_curve_residual_pc=float(np.sqrt(np.mean(np.sum(fit.fun.reshape(-1,3)**2,axis=1)))),
                         fitted_initial_cartesian_kpc_and_kpc_per_myr=fit.x[:6].tolist(),
                         assumed_last_lookback_myr=nominal_last,inferred_first_lookback_myr=float(lookback_first),
                         inferred_present_galpy_cartesian_kpc_and_kpc_per_myr=state.tolist(),inferred_present_icrs_by_frame=frames,
                         first_epoch_sensitivity=sensitivity,phase_free_skycheck=phase_free))
    out=dict(scope='Inversion of actual published vector centerline under MWPotential2014; post hoc diagnostic, not recovered execution inputs.',
             assumptions=['Path vertices are samples at a uniform interval; fitted interval checks this assumption.',
                          'The terminal Theia and OC vertices correspond to captioned7.43 and8.04Myr respectively; red terminal vertex shares OC8.04Myr. This is tested as an interpretation, not stated in source code.',
                          'Physical potential is MWPotential2014 withro8.122,vo220; source text does not identify all plotting settings.'],results=rows)
    (ROOT/'results/vector_inversion.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':run()
