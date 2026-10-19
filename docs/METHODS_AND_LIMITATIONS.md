# Methods and limits of the inference

The controlled comparison uses the paper's distance, 0.74 ± 0.06 kpc, and proper motions, (−2.57 ± 0.03, 2.08 ± 0.04) mas/yr, in both branches. The corrected branch replaces only RA and Dec with the primary timing position propagated from MJD 55045 to 2016.0. A separate nominal `timing` benchmark also uses more precise proper motions; it is not the coordinate-only headline experiment.

The two cluster rows come from the Hunt 2023 and 2024 VizieR catalogs. The selected astrometric quantities agree between those releases and match the rounded cluster table in Yu et al. The encounter thresholds are three times the catalog half-number radius: **31.06958874 pc** for Theia 58 and **53.34613749 pc** for OC 0450. Median age horizons are **23.59673145 Myr** and **21.74984819 Myr**. Cluster radial-velocity errors are errors on the means, not member dispersions.

The physical model uses `galpy.MWPotential2014`, distance scale 8.122 kpc and velocity scale 220 km/s. The Astropy Galactocentric frame uses solar height 20.8 pc and solar velocity [12.9, 245.6, 7.78] km/s. Handedness is checked against galpy's separate coordinate conversion. Backward integration uses `dop853_c`, relative tolerance 1e−10 and absolute tolerance 1e−12, with 0.005 Myr output spacing and quadratic interpolation of interior closest approaches.

Each Monte Carlo row contains one binary state shared by both cluster comparisons. Independent Gaussian measurement errors, symmetrized cluster distance errors, fixed catalog center positions, and fixed median cluster ages define the baseline uncertainty model. PCG64 seeds 2026101909, 2026101910, and 2026101911 generate the recorded draws. The code preserves the joint encounter flags rather than assuming independent cluster events.

The baseline radial-velocity prior is a zero-centered heliocentric Gaussian mixture with weights [0.463, 0.477, 0.060] and standard deviations [43.1, 107.6, 406.1] km/s. The uniform alternative is [−200, 200] km/s. These reproduce declared-prior interpretations; the original implementation's radial-frame convention remains unresolved. The direct conditional-peculiar Gaussian and Cauchy-family tests use solar-reflex and local-circular-motion corrections. They are separate model assumptions, not measurements of the unknown radial velocity.

## Numerical checks

- Halving the output spacing for the same first 1,024 histories changes no encounter decisions. The largest boundary-minimum difference comes from the final allowed output time, not an interior crossing.
- An independently written analytic-force Cartesian SciPy DOP853 integrator checks 32 deliberately selected histories in each coordinate branch. All **128** pairwise encounter decisions agree; the maximum minimum-separation difference is about **1.89×10⁻⁵ pc**. This validates an implementation of the same dynamics and shared frame, not an independent observation.
- A different plausible solar-velocity convention changes zero encounter decisions in 4,096 paired histories per branch. It improves the absolute published-cluster-curve comparison but does not exactly reproduce it.
- The figure inversion is checked with synthetic correct- and printed-coordinate trajectories, five starting intervals, and eight bounded quantization perturbations. These are local algorithmic checks, not a calibrated uncertainty posterior.

## Results that limit the claim

The printed-coordinate replay is close to, but not identical to, the paper. With a uniform radial prior, Theia 58 gives **3.8940%** versus the reported **3.69%** (about 1.88 combined Monte Carlo standard errors under a simple binomial comparison). OC 0450 gives **12.4329%** versus **12.42%**. Neither the frame nor the radial prior was tuned to remove that mismatch.

One closest point per encountered history does not match the paper's separation medians. A post hoc duration-weighted statistic, collecting all in-threshold time samples, matches them better. Radial-velocity medians match a different summary convention better. This unresolved bookkeeping is retained in `audit/results/encounter_statistics.json`; correlated post hoc matches must not be multiplied as independent evidence.

Using both cluster member proper-motion and radial-velocity dispersions, rather than errors on their means, yields **one OC 0450 encounter in 4,096 corrected-coordinate histories**. Member dispersion is not the preferred error on a cluster center, but this is a counterexample to an unconditional claim of zero probability. Other tested uncertainty alternatives and extended age horizons are also retained.

The actual phase-space arrays, numerical integration code, and complete plotting configuration for the original production run were not recovered. The inspected author source archive contains manuscript and figure assets, not the numerical execution bundle. This is a bounded inspection, not proof that no code exists elsewhere.

A corrected all-cluster search, full cluster-center/source covariance treatment, chance-alignment null model, and binary-supernova/ejection model have not been performed here. The earlier literature's use of different astrometry is not automatically invalidated. None of these calculations measures the true birthplace or establishes a revised formation-channel fraction.

## Failed controls retained

The original nominal trajectory did not reproduce the published Figure 2 encounter times. The cluster curves also disagreed, exposing a plotting-frame ambiguity; that mismatch alone was not attributed to the pulsar coordinate. The independent ODE's first minimum comparison failed because it integrated to the exact cluster age while the baseline stopped at the last allowed output sample; matching the time domains resolved it. Both limitations are reflected in the retained numerical outputs and the reproduction code.
