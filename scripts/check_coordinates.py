"""Check the coordinate discrepancy using only the Python standard library."""
from pathlib import Path
import json
import math
import re

ROOT = Path(__file__).resolve().parents[1]


def separation(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    a = math.sin((d2-d1)/2)**2 + math.cos(d1)*math.cos(d2)*math.sin((r2-r1)/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(1., max(0., a)))))


def check():
    p = json.loads((ROOT/'trajectory/data/J0737_primary_position.json').read_text())
    h, m, s = map(float, p['ra_hms'].split(':'))
    ra = 15*(h + m/60 + s/3600)
    d, m, s = map(float, p['dec_dms'].split(':'))
    dec = math.copysign(abs(d) + m/60 + s/3600, d)
    assert abs(ra-p['ra_deg']) < 1e-10
    assert abs(dec-p['dec_deg']) < 1e-10
    row = (ROOT/'sources/yu2026_table1_J0737.tsv').read_text().strip().split('\t')
    assert row[0] == 'J0737-3039AmidB'
    printed_ra, printed_dec = map(float, row[1:3])
    tex = (ROOT/'sources/yu2026_table1_J0737.tex').read_text()
    assert re.search(r'&\s*144\.448\s*&\s*\$-30\.6278\$', tex)
    angle = separation(printed_ra, printed_dec, ra, dec)
    dt = (51544.5 + 16*365.25 - p['position_epoch_MJD'])/365.25
    motion = math.hypot(p['mu_ra_cosdec_masyr'], p['mu_dec_masyr'])*dt/1000
    after_30 = separation(printed_ra-30, printed_dec, ra, dec)*3600
    assert 25 < angle < 26 and motion < .03 and after_30 > 120
    return dict(primary_ra_deg=ra, primary_dec_deg=dec,
                printed_ra_deg=printed_ra, printed_dec_deg=printed_dec,
                sky_separation_deg=angle, ra_coordinate_difference_deg=printed_ra-ra,
                approximate_motion_to_2016_arcsec=motion,
                residual_sky_error_after_subtracting_30_RA_degrees_arcsec=after_30,
                conclusion='Published position is incompatible with the primary timing position; an epoch correction cannot account for it.',
                limitation='Coordinates alone do not recover the inputs used by the original numerical executable.')


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
