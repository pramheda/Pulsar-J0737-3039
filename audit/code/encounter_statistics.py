"""Audit closest-point versus all-inside-time definitions against Table4."""
from pathlib import Path
import sys,json,time
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];CAMPAIGN=ROOT.parent
sys.path.insert(0,str(CAMPAIGN/'trajectory/code'))
from galactic_trace import sky,integrate,coordinates,TIME_MYR,CLUSTERS,BIN
import numpy as np

def run():
    start=time.monotonic();saved=np.load(CAMPAIGN/'trajectory/results/mc_published_n32768_step0.005.npz')
    z=np.random.default_rng(2026101909).standard_normal((32768,12));selected=np.flatnonzero(saved['events'].any(axis=1))
    clusters=[]
    for j,c in enumerate(CLUSTERS):
        o=4+4*j;ds=(c['distance_q84_pc']-c['distance_q16_pc'])/2000
        clusters.append((c['distance_pc']/1000+ds*z[:,o],c['pm_ra_cosdec_masyr']+c['pm_ra_standard_error_masyr']*z[:,o+1],
                         c['pm_dec_masyr']+c['pm_dec_standard_error_masyr']*z[:,o+2],c['heliocentric_radial_velocity_kms']+c['rv_standard_error_kms']*z[:,o+3]))
    tt=-np.arange(0,23.605,.005);natural=tt/TIME_MYR
    collected=[{k:[] for k in ['time','separation','relative_speed','rv','distance','cluster_speed_relative_present_Sun','binary_speed_relative_present_Sun']} for j in range(2)]
    for first in range(0,len(selected),64):
        ids=selected[first:first+64];m=len(ids)
        coord=sky(np.repeat([BIN['ra_deg']]+[c['ra_deg'] for c in CLUSTERS],m),np.repeat([BIN['dec_deg']]+[c['dec_deg'] for c in CLUSTERS],m),
                  np.concatenate([saved['binary_distance_kpc'][ids]]+[c[0][ids] for c in clusters]),
                  np.concatenate([-2.57+.03*z[ids,1]]+[c[1][ids] for c in clusters]),
                  np.concatenate([2.08+.04*z[ids,2]]+[c[2][ids] for c in clusters]),
                  np.concatenate([saved['radial_heliocentric_kms'][ids]]+[c[3][ids] for c in clusters]))
        o=integrate(coord,tt);xyz=coordinates(o,natural).reshape(3,m,len(tt),3)
        vel=np.stack([o.vx(natural),o.vy(natural),o.vz(natural)],axis=-1).reshape(3,m,len(tt),3)
        for j,c in enumerate(CLUSTERS):
            delta=np.linalg.norm(xyz[0]-xyz[j+1],axis=-1)*1000
            mask=(delta<c['encounter_threshold_3r50_pc'])&(-tt[None,:]<=c['age_Myr'])
            rows,times=np.nonzero(mask);d=collected[j]
            d['time'].append(-tt[times]);d['separation'].append(delta[rows,times])
            d['relative_speed'].append(np.linalg.norm(vel[0,rows,times]-vel[j+1,rows,times],axis=-1))
            d['cluster_speed_relative_present_Sun'].append(np.linalg.norm(vel[j+1,rows,times]-np.array([-12.9,245.6,7.78]),axis=-1))
            d['binary_speed_relative_present_Sun'].append(np.linalg.norm(vel[0,rows,times]-np.array([-12.9,245.6,7.78]),axis=-1))
            d['rv'].append(saved['radial_heliocentric_kms'][ids[rows]]);d['distance'].append(saved['binary_distance_kpc'][ids[rows]])
    results=[]
    for j,c in enumerate(CLUSTERS):
        kept=saved['events'][:,j]
        close={key:np.quantile(saved[key][kept,j] if saved[key].ndim==2 else saved[key][kept],[.16,.5,.84]).tolist() for key in ['minimum_pc','minimum_time_myr','radial_heliocentric_kms','binary_distance_kpc']}
        concat={k:np.concatenate(v) for k,v in collected[j].items()}
        allpoints={k:np.quantile(v,[.16,.5,.84]).tolist() for k,v in concat.items()}
        results.append(dict(cluster=c['name'],event_histories=int(kept.sum()),inside_timestep_count=len(concat['time']),closest_history_quantiles=close,all_inside_timestep_quantiles=allpoints))
    out=dict(scope='Reinterpretation of same retained printed-coordinate histories; no new observations or original execution.',
             definitions=['Closest point: one value per encounter history.','All inside time steps: every5kyr sample inside3r50, thereby duration-weighting histories.'],
             limitation='Paper does not unambiguously state which encounter-parameter summary statistic is used; agreement cannot prove implementation.',results=results,runtime_seconds=time.monotonic()-start)
    (ROOT/'results/encounter_statistics.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':run()
