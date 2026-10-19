# Claim-to-evidence map

| Claim | Inspect | Reproduce |
|---|---|---|
| The printed J0737 position conflicts with primary timing astrometry | `sources/yu2026_table1_J0737.tex`, `sources/yu2026_table1_J0737.tsv`, `trajectory/data/J0737_primary_position.json` | `python3 scripts/check_coordinates.py` |
| Correcting only RA/Dec removes the baseline sampled encounters | Paired `trajectory/results/mc_*_n32768_step0.005.{json,npz}` | `scripts/reproduce.py --full` |
| Uniform and conditional-peculiar priors were also tested | Named-prior JSON and NPZ files in `trajectory/results/` | `scripts/reproduce.py --full` |
| The numerical integrator and time resolution were checked | `validation/results/ode_check.json`, `trajectory/results/timestep_validation.json` | Full reproduction |
| The actual published curve is incompatible with the primary sky position under the declared model | `audit/sources/j0737_oc_orbits.pdf`, `audit/results/vector_figure.json`, `audit/results/vector_inversion.json` | Quick or full reproduction |
| Inversion can recover both correct and incorrect known-input positions | `audit/results/inversion_controls.json` | Quick or full reproduction |
| Plotting-frame and encounter-summary mismatches were not concealed | `audit/results/frame_counterchecks.json`, `audit/results/frame_mc.json`, `audit/results/encounter_statistics.json` | Full reproduction |
| Broader uncertainty assumptions do not justify universal zero probability | `audit/results/uncertainty_counterchecks.json` and individual arrays | Full reproduction |

The table, its TeX source, and its machine-readable rendering are versions of the same published information. They are not three independent measurements. Likewise, multiple numerical checks do not become multiple independent astronomical observations.

## Primary sources

1. Yu et al. (2026), *Kinematic Evidence for Open Cluster Origins of Galactic Binary Neutron Stars*, [ApJL 1000, L25](https://doi.org/10.3847/2041-8213/ae4bd1), [arXiv v1](https://arxiv.org/abs/2602.14145v1). Table 1: position; Tables 3–4: clusters and encounter summaries; Figure 2: trajectories. [Author source archive](https://arxiv.org/src/2602.14145v1).
2. [Publisher machine-readable Table 1](https://iopscience.iop.org/2041-8205/1000/1/L25/suppdata/apjlae4bd1t1_ascii.txt?doi=10.3847/2041-8213/ae4bd1). Only the J0737 row is redistributed here.
3. Kramer et al. (2021), [*Strong-Field Gravity Tests with the Double Pulsar*, arXiv:2112.06795](https://arxiv.org/abs/2112.06795), Tables II/IV. Position epoch MJD 55045; equinox J2000.
4. Hunt & Reffert catalogs: [2023, J/A+A/673/A114](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/673/A114/clusters) and [2024, J/A+A/686/A42](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/686/A42/clusters). Theia 58 ID 5495; OC 0450 ID 5015. The selected raw rows are in `sources/catalog/`.
5. [Catalog column definitions](https://cdsarc.cds.unistra.fr/ftp/J/A+A/673/A114/ReadMe): `e_RV`/`e_pm` are mean errors; `s_RV`/`s_pm` are member dispersions.

The exact source-document hashes and J0737 row locations are recorded in `provenance/EXPORT.json`. The original vector PDF has SHA-256 `6ff997a3e20be530dbb09b8c66082a6f8384b9db944e9c858b2a9127a564254f`.
