"""Verify hashes and numerical consistency of the exported J0737 evidence."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
from scipy.stats import beta
from check_coordinates import check as coordinate_check

ROOT = Path(__file__).resolve().parents[1]


def run(skip_hashes=False):
    checks = 0
    if not skip_hashes:
        for line in (ROOT/'SHA256SUMS').read_text().splitlines():
            sha, name = line.split('  ', 1)
            p = ROOT/name
            assert p.is_file() and not p.is_symlink(), name
            p.resolve().relative_to(ROOT)
            assert hashlib.sha256(p.read_bytes()).hexdigest() == sha, name
            checks += 1
    coords = coordinate_check()
    checks += 5
    base = ROOT/'trajectory/results'
    counts = {}
    thresholds = np.array([31.06958874, 53.34613749])
    for suffix in ['', '_uniform_heliocentric']:
        pair = []
        for branch in ['published', 'corrected_coordinates']:
            stem = f'mc_{branch}_n32768_step0.005{suffix}'
            with np.load(base/(stem+'.npz'), allow_pickle=False) as z:
                j = json.loads((base/(stem+'.json')).read_text())
                assert z['events'].shape == (32768, 2)
                assert z['events'].sum(axis=0).tolist() == j['raw_event_counts']
                assert np.array_equal(z['events'], z['minimum_pc'] < thresholds)
                pair.append({k:z[k].copy() for k in ['events','binary_distance_kpc','radial_heliocentric_kms']})
                counts[stem] = j['raw_event_counts']
                checks += 3
        for key in ['binary_distance_kpc','radial_heliocentric_kms']:
            assert np.array_equal(pair[0][key], pair[1][key]); checks += 1
        assert not pair[1]['events'].any(); checks += 1
    assert counts['mc_published_n32768_step0.005'] == [1725, 4514]
    assert counts['mc_published_n32768_step0.005_uniform_heliocentric'] == [1276, 4074]
    upper = float(beta.ppf(.95, 1, 32768))
    assert abs(upper - (1-.05**(1/32768))) < 1e-14
    checks += 3
    a = np.load(base/'mc_published_n32768_step0.005.npz', allow_pickle=False)
    assert int(a['events'].all(axis=1).sum()) == 1724; checks += 1
    for branch in ['published','corrected_coordinates']:
        a = np.load(base/f'mc_{branch}_n32768_step0.005.npz', allow_pickle=False)
        b = np.load(base/f'mc_{branch}_n1024_step0.0025.npz', allow_pickle=False)
        assert np.array_equal(a['events'][:1024], b['events']); checks += 1
    ode = json.loads((ROOT/'validation/results/ode_check.json').read_text())
    for branch in ode['branches']:
        assert branch['event_disagreements'] == 0
        assert branch['max_minimum_difference_pc'] < .0001
        checks += 2
    inv = json.loads((ROOT/'audit/results/vector_inversion.json').read_text())
    binary = next(x for x in inv['results'] if x['name'] == 'binary')
    assert 144.4 < binary['inferred_present_icrs_by_frame']['stated']['ra_deg'] < 144.7
    assert binary['phase_free_skycheck']['minimum_angular_separation_from_primary_deg'] > 24
    checks += 2
    extreme = np.load(ROOT/'audit/results/uncertainty_rv_pm_member_scatter.npz', allow_pickle=False)
    assert extreme['events'].sum(axis=0).tolist() == [0,1]; checks += 1
    vector = json.loads((ROOT/'audit/results/vector_figure.json').read_text())
    assert hashlib.sha256((ROOT/'audit'/vector['source']).read_bytes()).hexdigest() == vector['sha256']
    checks += 1
    return dict(status='PASS', checks=checks, hashes_checked=not skip_hashes,
                coordinate_check=coords, archived_event_counts=counts,
                zero_event_one_sided_95_upper_probability=upper,
                scope='Export integrity and retained-result consistency; fresh calculations are tested separately by reproduce.py.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--skip-hashes', action='store_true')
    args = p.parse_args()
    print(json.dumps(run(args.skip_hashes), indent=2))
