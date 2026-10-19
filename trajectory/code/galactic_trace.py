"""Conditional two-cluster traceback with explicit frames and shared binary states."""
from pathlib import Path
import os,sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT)
from astropy import units as u
from astropy.coordinates import SkyCoord,Galactocentric,CartesianDifferential
from astropy.time import Time
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,evaluateRforces,evaluatezforces
from galpy.util import conversion
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

RO=8.122;VO=220.;ZO=.0208
FRAME=Galactocentric(galcen_distance=RO*u.kpc,z_sun=ZO*u.kpc,
                    galcen_v_sun=CartesianDifferential([12.9,245.6,7.78]*u.km/u.s))
TIME_MYR=conversion.time_in_Gyr(VO,RO)*1000
CLUSTER_INPUT=ROOT/'data/key_clusters.json'
CLUSTERS=json.loads(CLUSTER_INPUT.read_text())['clusters']
BIN=dict(ra_deg=144.448,dec_deg=-30.6278,distance_kpc=.74,distance_sigma_kpc=.06,
         pmra=-2.57,pmra_sigma=.03,pmdec=2.08,pmdec_sigma=.04)

def sky(ra,dec,distance,pmra,pmdec,rv):
    return SkyCoord(ra=np.asarray(ra)*u.deg,dec=np.asarray(dec)*u.deg,distance=np.asarray(distance)*u.kpc,
                    pm_ra_cosdec=np.asarray(pmra)*u.mas/u.yr,pm_dec=np.asarray(pmdec)*u.mas/u.yr,
                    radial_velocity=np.asarray(rv)*u.km/u.s,frame='icrs')

def state(coord):
    gc=coord.transform_to(FRAME)
    # Astropy's +x points Sun->GC; galpy's +x points GC->Sun.
    xyz=gc.cartesian.xyz.to_value(u.kpc).T*np.array([-1,1,1])
    vel=gc.velocity.d_xyz.to_value(u.km/u.s).T*np.array([-1,1,1])
    R=np.hypot(xyz[...,0],xyz[...,1]);phi=np.arctan2(xyz[...,1],xyz[...,0])
    vR=(xyz[...,0]*vel[...,0]+xyz[...,1]*vel[...,1])/R
    vT=(xyz[...,0]*vel[...,1]-xyz[...,1]*vel[...,0])/R
    return np.stack([R/RO,vR/VO,vT/VO,xyz[...,2]/RO,vel[...,2]/VO,phi],axis=-1)

def make_orbit(coord):return Orbit(state(coord),ro=RO,vo=VO,zo=ZO)

def coordinates(orbit,t):
    return np.stack([orbit.x(t),orbit.y(t),orbit.z(t)],axis=-1)

def integrate(coord,t_myr,rtol=1e-10,atol=1e-12):
    orbit=make_orbit(coord);orbit.integrate(np.asarray(t_myr)/TIME_MYR,MWPotential2014,
          method='dop853_c',rtol=rtol,atol=atol,progressbar=False,numcores=1)
    return orbit

def timing_binary():
    c=SkyCoord(ra='07h37m51.248115s',dec='-30d39m40.70485s',distance=.74*u.kpc,
               pm_ra_cosdec=-2.567*u.mas/u.yr,pm_dec=2.082*u.mas/u.yr,
               radial_velocity=66.25*u.km/u.s,obstime=Time(55045.,format='mjd'))
    c=c.apply_space_motion(new_obstime=Time(2016.,format='jyear'))
    return dict(BIN,ra_deg=float(c.ra.deg),dec_deg=float(c.dec.deg),pmra=-2.567,pmdec=2.082,
                position_source='Kramer2021TableII, MJD55045 propagated to2016.0; rounded distance0.74kpc retained')

def benchmark(position='published'):
    start=time.monotonic();rows=[]
    b=BIN if position=='published' else timing_binary()
    binary=sky(b['ra_deg'],b['dec_deg'],b['distance_kpc'],b['pmra'],b['pmdec'],66.25)
    allcoord=sky([b['ra_deg']]+[c['ra_deg'] for c in CLUSTERS],
                 [b['dec_deg']]+[c['dec_deg'] for c in CLUSTERS],
                 [b['distance_kpc']]+[c['distance_pc']/1000 for c in CLUSTERS],
                 [b['pmra']]+[c['pm_ra_cosdec_masyr'] for c in CLUSTERS],
                 [b['pmdec']]+[c['pm_dec_masyr'] for c in CLUSTERS],
                 [66.25]+[c['heliocentric_radial_velocity_kms'] for c in CLUSTERS])
    for step in [.01,.005,.0025]:
        tt=-np.arange(0,24+step/2,step);o=integrate(allcoord,tt);t=tt/TIME_MYR
        xyz=coordinates(o,t);E=o.E(t);L=o.Lz(t)
        for j,c in enumerate(CLUSTERS):
            dist=np.linalg.norm(xyz[0]-xyz[j+1],axis=-1)*1000;k=int(np.argmin(dist))
            objective=lambda x:float(np.linalg.norm(coordinates(o[0],x/TIME_MYR)-coordinates(o[j+1],x/TIME_MYR))*1000)
            opt=minimize_scalar(objective,bounds=sorted([tt[max(0,k-2)],tt[min(k+2,len(tt)-1)]]),method='bounded',options={'xatol':1e-9})
            rows.append(dict(step_myr=step,cluster=c['name'],grid_min_separation_pc=float(dist[k]),grid_lookback_myr=float(-tt[k]),
                             interpolated_min_separation_pc=float(opt.fun),interpolated_lookback_myr=float(-opt.x),
                             threshold_pc=c['encounter_threshold_3r50_pc'],is_encounter=bool(opt.fun<c['encounter_threshold_3r50_pc'])))
        energy=float(np.max(abs((E-E[:,0,None])/E[:,0,None])))
        angular=float(np.max(abs((L-L[:,0,None])/L[:,0,None])))
        assert energy<1e-8 and angular<1e-8
    # Cartesian solve_ivp is an independent integration algorithm using same potential.
    s=state(binary);R,vR,vT,z,vz,phi=s
    initial=[R*np.cos(phi),R*np.sin(phi),z,vR*np.cos(phi)-vT*np.sin(phi),vR*np.sin(phi)+vT*np.cos(phi),vz]
    def rhs(t,s):
        x,y,z,vx,vy,vz=s;radius=np.hypot(x,y)
        fR=evaluateRforces(MWPotential2014,radius,z,use_physical=False)
        fz=evaluatezforces(MWPotential2014,radius,z,use_physical=False)
        return [vx,vy,vz,fR*x/radius,fR*y/radius,fz]
    sol=solve_ivp(rhs,(0,-24/TIME_MYR),initial,t_eval=t,method='DOP853',rtol=2e-12,atol=2e-14)
    assert sol.success
    discrepancy=float(np.max(np.linalg.norm(sol.y[:3].T*RO-xyz[0],axis=-1))*1000)
    assert discrepancy<1e-4
    # Coordinate handedness check against galpy's separate SkyCoord input path.
    native=Orbit(binary.transform_to(FRAME),ro=np.sqrt(RO**2-ZO**2),vo=VO,zo=ZO)
    manual=make_orbit(binary)
    coordinate_discrepancy=float(np.max(abs(coordinates(native,0)-coordinates(manual,0))))
    assert coordinate_discrepancy<1e-10
    out=dict(scope='Fixed published Table1 binary state, catalog cluster centers, heliocentric vr66.25. Figure2 precise draw parameters unavailable.',
             binary=b,position_choice=position,radial_heliocentric_kms=66.25,published_figure2_lookback_myr={'Theia_58':7.43,'OC_0450':8.04},
             frame=dict(galcen_distance_kpc=RO,z_sun_kpc=ZO,galcen_v_sun_kms=[12.9,245.6,7.78],galpy_ro_kpc=RO,galpy_vo_kms=VO),
             result=rows,max_relative_energy_drift=energy,max_relative_Lz_drift=angular,
             independent_integrator_max_discrepancy_pc=discrepancy,coordinate_conversion_max_discrepancy_kpc=coordinate_discrepancy,
             runtime_seconds=time.monotonic()-start)
    (ROOT/'results'/f'benchmark_{position}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    benchmark('published')
    benchmark('timing')
