# Pulsar J0737−3039: coordinate-error audit of Yu et al. (2026)

**The position printed for J0737−3039 in Yu et al. (2026) is wrong by about 25.719° on the sky. Changing only that position removes the sampled encounters with both highlighted birth-cluster candidates in the controlled calculations archived here. The numerical geometry of the published trajectory figure also points to the wrong sky region.**

This repository contains the J0737−3039 evidence, analysis code, saved Monte Carlo histories, and numerical controls. The paper is [*Kinematic Evidence for Open Cluster Origins of Galactic Binary Neutron Stars*, ApJL **1000**, L25](https://doi.org/10.3847/2041-8213/ae4bd1), by Wen-Jie Yu and colleagues; the inspected manuscript is [arXiv:2602.14145v1](https://arxiv.org/abs/2602.14145v1).

**What is established:** the published coordinate error, its causal effect in an explicit two-cluster reanalysis, and a strong conditional inconsistency between the published figure and the primary pulsar position. **What is not established:** the exact inputs used by every original calculation. The authors' numerical execution bundle was not recovered. This is an independent audit, not an author-issued correction.

## 1. The coordinate discrepancy

| Source | Right ascension (degrees) | Declination (degrees) |
|---|---:|---:|
| Yu et al., Table 1 | 144.448 | −30.6278 |
| Kramer et al. (2021), primary timing position, MJD 55045 | 114.4635338125 | −30.6613069028 |
| Primary position propagated to 2016.0 for this experiment | 114.4635284939 | −30.6613031921 |

The erroneous values appear in the [author's TeX row](sources/yu2026_table1_J0737.tex) and the [publisher's machine-readable row](sources/yu2026_table1_J0737.tsv), as well as the manuscript table. The [primary timing input and citation](trajectory/data/J0737_primary_position.json) give RA **07:37:51.248115** and Dec **−30:39:40.70485**. J2000 is the equinox, not the position epoch. Proper motion to 2016 changes the position by only about **0.021 arcseconds**.

Run the coordinate check without installing any scientific packages:

```bash
python3 scripts/check_coordinates.py
```

## 2. Changing only the coordinates changes the encounter result

The Gaussian-prior experiment uses **32,768 paired histories**. Both branches share the distance, proper motions, uncertainties, random draws, cluster inputs, Galactic potential, and encounter criterion. Only the pulsar's RA and Dec change.

| Cluster | Paper's reported fraction | Replay with printed coordinates | Same replay with corrected coordinates |
|---|---:|---:|---:|
| Theia 58 | 5.35% | 1,725 / 32,768 = **5.2643%** | **0 / 32,768** |
| OC 0450 | 13.87% | 4,514 / 32,768 = **13.7756%** | **0 / 32,768** |

For each zero-event result, the one-sided 95% binomial upper bound is **0.00914%**, conditional on that sampled model. Zero sampled encounters does **not** establish zero physical probability. An encounter fraction is also not a posterior probability that a cluster is the birthplace.

![Paired coordinate experiment: minimum past separation from Theia 58 and OC 0450](trajectory/results/coordinate_effect.png)

The independent uniform-radial-prior experiment also gives **0 / 32,768** corrected-coordinate encounters for both clusters. Its printed-coordinate Theia result differs modestly from the paper; that mismatch is retained. The complete arrays and summaries are in [trajectory/results](trajectory/results/).

## 3. Evidence beyond a table typo: the actual published figure

The audit extracts vector vertices directly from the [author-supplied Figure 2 PDF](audit/sources/j0737_oc_orbits.pdf). Under the declared Galactic potential and spatial frame, fitting its red trajectory gives a present-day RA near **144.526°**, rather than **114.464°**, for the tested caption-based timing interpretation.

Freeing the curve's temporal placement over 0–24 Myr still leaves a minimum separation of **24.252°** from the primary pulsar position. Synthetic trajectories at both the correct and printed positions, multiple fitting starts, and bounded figure-quantization perturbations test the inversion method.

This is strong evidence that the plotted numerical trajectory is inconsistent with the real J0737−3039 position under the tested model. It is **not** an exact recovery of the authors' inputs or proof of global uniqueness across every possible plotting convention. The initial nominal figure replay failed, including for the cluster curves; the frame counterchecks documenting that failure are retained.

Inspect [vector extraction](audit/results/vector_figure.json), [trajectory inversion](audit/results/vector_inversion.json), [synthetic recovery controls](audit/results/inversion_controls.json), and [frame alternatives](audit/results/frame_counterchecks.json).

## Reproduce and inspect

Python 3.12 was used. The recorded package versions are pinned in [requirements.txt](requirements.txt).

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/verify.py
.venv/bin/python scripts/reproduce.py --quick
```

The quick run recalculates paired 512-history trajectories, checks them against the corresponding archived histories, re-extracts the figure vertices, reruns the figure inversion and its controls, and checks selected trajectories with the independent ODE implementation. Recalculated files go into the ignored `reproduced/` directory; the archived evidence is preserved.

```bash
# Recalculate the full archived experiment and robustness checks:
.venv/bin/python scripts/reproduce.py --full
```

See [reproduction details](docs/REPRODUCING.md), [methods and limitations](docs/METHODS_AND_LIMITATIONS.md), and the [claim-to-evidence map](docs/EVIDENCE.md). [SHA256SUMS](SHA256SUMS) and [export provenance](provenance/EXPORT.json) identify the exact files. Third-party figure attribution and licensing are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Scientific implication and scope

The claimed J0737−3039 encounters with **Theia 58 and OC 0450 require reassessment**: the coordinate-only counterfactual removes them in the declared model, and the published figure provides additional evidence of a propagated position problem. This repository does not establish a new birthplace, revise all clusters or other binary systems, measure a natal kick, or overturn a binary-neutron-star formation channel in general.

One broader-uncertainty test produced **one** corrected-coordinate OC 0450 encounter in 4,096 histories; it is included. That result rules out describing the finding as universal impossibility. Exact execution provenance, a corrected full-catalog search, and a separately defined physical-origin/null model remain necessary for stronger conclusions.

Audit performed on 19 October 2026 with AI-assisted scientific computation. Agreement between analysis branches is treated as a check of implementations, not as independent astronomical observations. No author endorsement is implied.
