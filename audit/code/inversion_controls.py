"""Known-input recovery, multistart, and figure quantization sensitivity."""
from pathlib import Path
import sys,json
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
from invert_vector_tracks import integrate,icrs,K
import numpy as np
from scipy.optimize import least_squares
from astropy.coordinates import SkyCoord,Galactocentric,CartesianDifferential
from astropy import units as u

def fit_track(points,dt0=.05):
    initial=np.r_[points[0],-(points[1]-points[0])/dt0,dt0]
    fun=lambda p:((integrate(p[:6],-np.arange(len(points))*p[6])[:,:3]-points)*1000).ravel()
    fit=least_squares(fun,initial,bounds=(np.r_[[-np.inf]*6,.02],np.r_[[np.inf]*6,.08]),x_scale=np.r_[[1]*3,[.1]*3,.05],xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=200)
    # SVD only describes local identifiability under the fixed physical potential.
    scaled_jac=fit.jac*np.r_[[1]*3,[.1]*3,.05]
    singular=np.linalg.svd(scaled_jac,compute_uv=False)
    now=integrate(fit.x[:6],np.array([0.,8.04-(len(points)-1)*fit.x[6]]))[-1]
    return dict(start_dt=dt0,fitted_dt=float(fit.x[6]),rms_pc=float(np.sqrt(np.sum(fit.fun**2)/len(points))),
                scaled_jac_singular_values=singular.tolist(),parameters=fit.x.tolist(),reconstructed_icrs=icrs(now,np.array([12.9,245.6,7.78])))

def run():
    frame=Galactocentric(galcen_distance=8.122*u.kpc,z_sun=.0208*u.kpc,galcen_v_sun=CartesianDifferential([12.9,245.6,7.78]*u.km/u.s))
    quantum=np.array([.0025,.010,.0015])/1000
    output=[]
    for name,ra,dec in [('printed',144.448,-30.6278),('correct',114.46352849386929,-30.661303192105727)]:
        c=SkyCoord(ra=ra*u.deg,dec=dec*u.deg,distance=.74*u.kpc,pm_ra_cosdec=-2.57*u.mas/u.yr,pm_dec=2.08*u.mas/u.yr,radial_velocity=66.25*u.km/u.s).transform_to(frame)
        state=np.r_[c.cartesian.xyz.to_value(u.kpc)*[-1,1,1],c.velocity.d_xyz.to_value(u.km/u.s)*[-1,1,1]/K]
        lookback=np.r_[0.,6.59+np.arange(30)*.05]
        exact=integrate(state,-lookback)[1:,:3]
        quantized=np.round(exact/quantum)*quantum
        control=fit_track(quantized)
        output.append(dict(name=name,input=dict(ra=ra,dec=dec,d=.74,pmra=-2.57,pmdec=2.08,rv=66.25),true_dt=.05,
                           coordinate_quantum_pc=(quantum*1000).tolist(),max_quantization_error_pc=(np.max(abs(quantized-exact),axis=0)*1000).tolist(),result=control))
    d=json.loads((ROOT/'results/vector_figure.json').read_text())['paths']['binary']
    xy=np.array(d['xy']);xz=np.array(d['xz']);points=np.column_stack([(xy[:,0]+xz[:,0])/2,xy[:,1],xz[:,1]])
    starts=[fit_track(points,dt) for dt in [.025,.04,.05,.065,.079]]
    rng=np.random.default_rng(2026101931)
    perturb=[fit_track(points+rng.uniform(-.5,.5,points.shape)*quantum) for i in range(8)]
    out=dict(scope='Synthetic known-input recovery validates an inversion method; it is not additional astronomical evidence.',
             synthetic_controls=output,figure_multistart=starts,figure_bounded_quantization_perturbations=perturb,
             limitations=['Small local multistart study does not prove global uniqueness across different potentials or time samplings.',
                          'Sky-coordinate recovery assumes the final point has lookback8.04Myr.',
                          'Quantization perturbations test numerical sensitivity, not measurement errors or a calibrated confidence interval.',
                          'Potential normalization and uniform time sampling remain assumptions.'])
    (ROOT/'results/inversion_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'synthetic':output,'multistart_dt':[x['fitted_dt'] for x in starts],
                      'figure_ra_range':[min(x['reconstructed_icrs']['ra_deg'] for x in perturb),max(x['reconstructed_icrs']['ra_deg'] for x in perturb)],
                      'figure_dt_range':[min(x['fitted_dt'] for x in perturb),max(x['fitted_dt'] for x in perturb)]},indent=2))
if __name__=='__main__':run()
