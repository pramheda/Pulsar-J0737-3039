"""Compare explicit frame choices against the actual published vector curves."""
from pathlib import Path
import sys,json
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];CAMPAIGN=ROOT.parent
sys.path.insert(0,str(CAMPAIGN/'trajectory/code'))
import galactic_trace as g
import numpy as np
from scipy.optimize import minimize_scalar
from astropy.coordinates import Galactocentric,CartesianDifferential
from astropy import units as u
from galpy.util import conversion

def run():
    data=json.loads((ROOT/'results/vector_figure.json').read_text());out=[]
    for label,radial_scale,frame_distance,zo,solar in [
        ('stated',8.122,8.122,.0208,[12.9,245.6,7.78]),
        ('schoenrich_solar',8.122,8.122,.0208,[11.1,232.24,7.25]),
        ('schoenrich_solar_ro8',8.,8.122,.0208,[11.1,232.24,7.25]),
        ('schoenrich_solar_R8_ro8',8.,8.,.0208,[11.1,232.24,7.25]),
        ('schoenrich_solar_z25',8.122,8.122,.025,[11.1,232.24,7.25]),
        ('schoenrich_solar_ro8_z25',8.,8.122,.025,[11.1,232.24,7.25]),
        ('stated_solar_U_signflip',8.122,8.122,.0208,[-12.9,245.6,7.78]),
    ]:
        g.RO=radial_scale;g.TIME_MYR=conversion.time_in_Gyr(220.,radial_scale)*1000
        g.FRAME=Galactocentric(galcen_distance=frame_distance*u.kpc,z_sun=zo*u.kpc,galcen_v_sun=CartesianDifferential(solar*u.km/u.s))
        times=-np.linspace(0,16,1601);rows=[]
        for name,path in data['paths'].items():
            xy=np.array(path['xy']);xz=np.array(path['xz'])
            target=np.stack([(xy[:,0]+xz[:,0])/2,xy[:,1],xz[:,1]],axis=1)
            def compare(coord):
                orbit=g.integrate(coord,times);xyz=g.coordinates(orbit,times/g.TIME_MYR)
                assert np.all(np.diff(xyz[:,0])<0)
                fit=np.stack([np.interp(target[:,0],xyz[::-1,0],xyz[::-1,j]) for j in range(3)],axis=1)
                return dict(rms_pc=float(np.sqrt(np.mean(np.sum((fit-target)**2,axis=1)))*1000),
                            xyz_residual_mean_pc=np.mean((fit-target)*1000,axis=0).tolist(),
                            xyz_at_labelled_times=g.coordinates(orbit,-np.array([7.43,8.04])/g.TIME_MYR).tolist())
            if name!='binary':
                c=next(c for c in g.CLUSTERS if c['name']==name)
                row=compare(g.sky(c['ra_deg'],c['dec_deg'],c['distance_pc']/1000,c['pm_ra_cosdec_masyr'],c['pm_dec_masyr'],c['heliocentric_radial_velocity_kms']))
                rows.append(dict(name=name,**row))
            else:
                for branch,b in [('printed',g.BIN),('correct',dict(g.BIN,ra_deg=g.timing_binary()['ra_deg'],dec_deg=g.timing_binary()['dec_deg']))]:
                    fun=lambda d:compare(g.sky(b['ra_deg'],b['dec_deg'],d,b['pmra'],b['pmdec'],66.25))
                    opt=minimize_scalar(lambda d:fun(d)['rms_pc'],bounds=(.38,1.1),method='bounded',options={'xatol':1e-8})
                    rows.append(dict(name=name,position=branch,nominal=fun(.74),best_distance_kpc=float(opt.x),best_distance_rms_pc=float(opt.fun)))
        out.append(dict(name=label,potential_ro_kpc=radial_scale,frame_distance_kpc=frame_distance,z_sun_kpc=zo,solar_velocity_kms=solar,comparisons=rows))
    (ROOT/'results/frame_counterchecks.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in out:print(row['name'],[{k:v for k,v in r.items() if k in ['name','position','rms_pc','best_distance_kpc','best_distance_rms_pc']} for r in row['comparisons']],flush=True)

if __name__=='__main__':run()
