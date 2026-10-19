"""Small paired frame-convention control; underlying source draws unchanged."""
from pathlib import Path
import sys,json
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];CAMPAIGN=ROOT.parent
sys.path.insert(0,str(CAMPAIGN/'trajectory/code'))
import galactic_trace as g
from uncertainty_counterchecks import run
from astropy.coordinates import Galactocentric,CartesianDifferential
from astropy import units as u
import numpy as np

g.FRAME=Galactocentric(galcen_distance=8.122*u.kpc,z_sun=.0208*u.kpc,galcen_v_sun=CartesianDifferential([11.1,232.24,7.25]*u.km/u.s))
out=[]
for position in ['printed','correct']:
    r=run('schoenrich_frame_'+position,n=4096,position=position)
    old=np.load(CAMPAIGN/'trajectory/results'/('mc_published_n32768_step0.005.npz' if position=='printed' else 'mc_corrected_coordinates_n32768_step0.005.npz'))
    new=np.load(ROOT/'results'/f'uncertainty_schoenrich_frame_{position}.npz')
    r['stated_frame_first4096_counts']=old['events'][:4096].sum(axis=0).tolist()
    r['paired_flag_changes']=np.count_nonzero(old['events'][:4096]!=new['events'],axis=0).tolist()
    r['solar_velocity_kms']=[11.1,232.24,7.25]
    r['max_minimum_change_pc']=np.max(abs(old['minimum_pc'][:4096]-new['minimum_pc']),axis=0).tolist()
    out.append(r)
(ROOT/'results/frame_mc.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
