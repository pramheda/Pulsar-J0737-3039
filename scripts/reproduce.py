"""Run the standalone audit in an isolated copy, preserving archived evidence."""
from pathlib import Path
import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def run(full=False):
    start = time.monotonic()
    name = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    dest = ROOT/'reproduced'/name
    dest.mkdir(parents=True)
    for folder in ['trajectory','audit','validation']:
        shutil.copytree(ROOT/folder, dest/folder, ignore=shutil.ignore_patterns('__pycache__'))
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               PYTHONDONTWRITEBYTECODE='1', MPLCONFIGDIR=str(dest/'matplotlib-cache'))
    commands = []
    def execute(path, *args):
        command = [sys.executable, '-B', str(dest/path), *map(str,args)]
        print('Running', path, *args, flush=True)
        log = dest/('run_'+str(len(commands)).zfill(2)+'.log')
        with log.open('w') as f:
            result = subprocess.run(command,cwd=dest/'trajectory',env=env,stdout=f,stderr=subprocess.STDOUT)
        commands.append(dict(script=path,arguments=list(map(str,args)),returncode=result.returncode,log=log.name))
        if result.returncode:
            raise RuntimeError(f'{path} failed; inspect {log}')

    n = 32768 if full else 512
    for branch in ['published','corrected_coordinates']:
        execute('trajectory/code/monte_carlo.py','--n',n,'--position',branch)
        if full:
            execute('trajectory/code/monte_carlo.py','--n',n,'--position',branch,'--prior','uniform_heliocentric')
            execute('trajectory/code/monte_carlo.py','--n',1024,'--position',branch,'--step',.0025)
    comparisons = []
    for branch in ['published','corrected_coordinates']:
        fresh = np.load(dest/'trajectory/results'/f'mc_{branch}_n{n}_step0.005.npz',allow_pickle=False)
        old = np.load(ROOT/'trajectory/results'/f'mc_{branch}_n32768_step0.005.npz',allow_pickle=False)
        events_equal = np.array_equal(fresh['events'],old['events'][:n])
        max_error = float(np.max(abs(fresh['minimum_pc']-old['minimum_pc'][:n])))
        assert events_equal and max_error < .01, (branch, events_equal, max_error)
        comparisons.append(dict(branch=branch,n=n,event_flags_equal=events_equal,max_minimum_difference_pc=max_error))
    execute('audit/code/figure_vector.py')
    execute('audit/code/invert_vector_tracks.py')
    execute('audit/code/inversion_controls.py')
    old = json.loads((ROOT/'audit/results/vector_inversion.json').read_text())['results'][0]
    fresh = json.loads((dest/'audit/results/vector_inversion.json').read_text())['results'][0]
    ra_difference = abs(old['inferred_present_icrs_by_frame']['stated']['ra_deg']-fresh['inferred_present_icrs_by_frame']['stated']['ra_deg'])
    sky_difference = abs(old['phase_free_skycheck']['minimum_angular_separation_from_primary_deg']-fresh['phase_free_skycheck']['minimum_angular_separation_from_primary_deg'])
    assert ra_difference < .001 and sky_difference < .001, (ra_difference,sky_difference)
    execute('validation/code/ode_check.py')
    if full:
        for prior in ['gaussian_peculiar_conditional','cauchy_peculiar_conditional']:
            execute('trajectory/code/monte_carlo.py','--n',8192,'--position','corrected_coordinates','--prior',prior)
        execute('trajectory/code/galactic_trace.py')
        execute('trajectory/code/validate_saved.py')
        execute('trajectory/code/plot_results.py')
        for script in ['frame_counterchecks','encounter_statistics','uncertainty_counterchecks','frame_mc']:
            execute('audit/code/'+script+'.py')
    versions = {p:importlib.metadata.version(p) for p in ['numpy','scipy','astropy','matplotlib','galpy']}
    report = dict(status='PASS',mode='full' if full else 'quick',python=sys.version,
                  packages=versions,commands=commands,paired_reproduction=comparisons,
                  figure_inversion_ra_difference_deg=ra_difference,
                  figure_minimum_sky_separation_difference_deg=sky_difference,
                  seconds=time.monotonic()-start,
                  scope='Fresh rerun of this independent audit; does not recover the paper authors\' original execution.')
    (dest/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    print('Saved:',dest/'REPRODUCTION.json')


if __name__ == '__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--quick',action='store_true');g.add_argument('--full',action='store_true')
    a=p.parse_args();run(full=a.full)
