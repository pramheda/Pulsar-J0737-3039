# J0737−3039: a sky-position error in a published birthplace study

A 2026 astronomy paper tried to find where the double pulsar **J0737−3039** was born. It traced the system's motion backwards to see whether it passed close to star clusters that could have been its birthplace.

**The paper puts the pulsar in the wrong place on the sky—about 25.7° away from its measured position. In our reconstruction, correcting that position removes the close passes to both proposed clusters in the main simulations.**

This matters because those close passes are the evidence for linking the pulsar to the clusters **Theia 58** and **OC 0450**. Our results show that these proposed birthplace links need to be reassessed. They do not tell us where the pulsar was actually born.

This repository contains the source records, code, saved results, and checks behind that finding. The paper examined is [*Kinematic Evidence for Open Cluster Origins of Galactic Binary Neutron Stars*](https://doi.org/10.3847/2041-8213/ae4bd1), by Yu and colleagues, published in *The Astrophysical Journal Letters* in 2026 ([manuscript examined](https://arxiv.org/abs/2602.14145v1)).

## What we found

| Check | Result | Why it matters |
|---|---|---|
| Compare the paper's sky position with the measured position | They differ by about **25.7°**. | The published position is wrong. |
| Trace possible past paths using the paper's position | Close passes to both clusters occur at rates similar to those reported in the paper. | Our separate calculation approximately reproduces the reported result. |
| Repeat with only the sky position corrected | **No close passes** to either cluster in the main samples of 32,768 possible paths. | The position error can change the birthplace inference. |
| Examine the path drawn in the paper's own figure | Under the model tested, it points toward the wrong sky region too. | There is evidence that the problem affected the plotted calculation, beyond the table entry. |

**Certainty and limits:** The incorrect published position is directly checkable. Its effect is demonstrated in our reconstruction. The figure provides further evidence, but we do not have the authors' original calculation files, so we cannot establish exactly what went into every original run. This is an independent audit, not a correction issued by the authors.

## What changes when the position is corrected?

We simulated **32,768 possible past paths** to account for uncertain measurements and unknown motion toward or away from us. For each path, we ran the calculation twice: once with the paper's position and once with the corrected position. Everything else—including the random draws—stayed the same.

The table shows how often the pulsar passed close enough to each cluster to meet the experiment's encounter criterion. These are **fractions of simulated paths**, not probabilities that a cluster really was the birthplace.

| Possible birthplace | Paper's reported fraction | Our calculation: paper's position | Our calculation: corrected position |
|---|---:|---:|---:|
| Theia 58 | 5.35% | **5.26%** (1,725 of 32,768) | **0 of 32,768** |
| OC 0450 | 13.87% | **13.78%** (4,514 of 32,768) | **0 of 32,768** |

![How correcting the sky position changes the closest past approach to Theia 58 and OC 0450](trajectory/results/coordinate_effect.png)

These results use the paper's Gaussian-mixture assumption for the unknown radial velocity. A second experiment using a uniform radial-velocity distribution also found no corrected-position encounters in 32,768 paths for either cluster. Its result for Theia 58 with the paper's position differed modestly from the published number; that mismatch is documented.

**Zero in a simulation does not mean impossible in nature.** In a separate test allowing broader uncertainties in the cluster motions, one corrected-position path out of 4,096 passed close to OC 0450. The conclusion is that the reported links are sensitive to the position error and require reassessment, not that these clusters are ruled out under every assumption.

## The evidence in detail

### 1. The published sky position is wrong

Astronomers specify a sky position using right ascension and declination, similar to longitude and latitude.

| Source | Right ascension (degrees) | Declination (degrees) |
|---|---:|---:|
| Yu et al. (2026), Table 1 | 144.448 | −30.6278 |
| Kramer et al. (2021), measured timing position at MJD 55045 | 114.4635338125 | −30.6613069028 |
| Measured position adjusted to 2016.0 for our calculation | 114.4635284939 | −30.6613031921 |

The incorrect values appear in the [manuscript's table source](sources/yu2026_table1_J0737.tex) and the [publisher's machine-readable table](sources/yu2026_table1_J0737.tsv). The measured position comes from [Kramer et al. (2021)](https://arxiv.org/abs/2112.06795); the exact input and citation are recorded in [J0737_primary_position.json](trajectory/data/J0737_primary_position.json).

The separation on the sky is **25.719°**. The pulsar's motion between the measurement epoch and 2016 changes its position by only about **0.021 arcseconds**, far too little to explain the discrepancy. J2000 specifies the coordinate equinox, not the measurement date. We have not established how the wrong numbers entered the paper.

You can check the coordinate discrepancy with Python alone, without installing scientific packages:

```bash
python3 scripts/check_coordinates.py
```

### 2. The paper's own trajectory figure also points to a problem

We extracted the numerical drawing coordinates directly from the [published Figure 2 PDF](audit/sources/j0737_oc_orbits.pdf), then fitted an orbital model to the pulsar's drawn path.

With the tested Galactic model, coordinate frame, and interpretation of the figure's timing, the fitted present-day right ascension is about **144.526°**—near the paper's incorrect value, rather than the measured **114.464°**. Allowing the curve's time placement to vary over 0–24 million years still leaves the inferred sky position at least **24.252°** from the measured position.

This supports the interpretation that the error affected the plotted calculation. It does not recover the authors' exact inputs or exclude every possible plotting convention. We checked the fitting method using artificial trajectories with known inputs, different starting guesses, and small perturbations to the extracted figure coordinates.

We also retain checks that did not succeed: the initial attempt to reproduce the figure failed for the cluster paths as well as the pulsar path. Alternative coordinate-frame choices improved the comparison without matching it exactly. Those failures limit what can be claimed from a figure mismatch alone.

The numerical records are in [figure extraction](audit/results/vector_figure.json), [orbit fitting](audit/results/vector_inversion.json), [known-input checks](audit/results/inversion_controls.json), and [coordinate-frame checks](audit/results/frame_counterchecks.json).

### 3. What the result does—and does not—change

The finding concerns the proposed link between **J0737−3039 and two specific clusters**. It does not establish a new birthplace, test every cluster, revise other binary systems, or by itself overturn a general theory of how double neutron stars form.

For each zero-encounter result in the main experiment, the one-sided 95% binomial upper bound is **0.00914%**, conditional on the sampled model. This bounds that model's encounter fraction; it is not a bound on the physical probability of a birthplace.

Stronger astrophysical conclusions would require the authors' original calculation files, a corrected search across the full cluster catalogue, and a model that compares genuine birth associations with chance close passes. See [methods and limitations](docs/METHODS_AND_LIMITATIONS.md) for the assumptions, sensitivity tests, and unresolved questions.

## Run the checks yourself

Python 3.12 was used. Package versions are recorded in [requirements.txt](requirements.txt).

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/verify.py
.venv/bin/python scripts/reproduce.py --quick
```

`verify.py` checks the saved evidence and file checksums. The quick reproduction recalculates 512 paired paths, compares them with the archived results, repeats the figure analysis and its checks, and checks selected paths with a separately implemented orbit integrator. New results go into `reproduced/`.

To recalculate the full archived experiment and sensitivity tests:

```bash
.venv/bin/python scripts/reproduce.py --full
```

## Where to find things

| What you want to inspect | Where to look |
|---|---|
| Evidence supporting each claim | [Evidence guide](docs/EVIDENCE.md) |
| Assumptions, limitations, and unsuccessful checks | [Methods and limitations](docs/METHODS_AND_LIMITATIONS.md) |
| Instructions for repeating the calculations | [Reproduction guide](docs/REPRODUCING.md) |
| Saved path simulations and summaries | [Trajectory results](trajectory/results/) |
| File checksums and records of where the files came from | [SHA256SUMS](SHA256SUMS), [export record](provenance/EXPORT.json) |
| Attribution and licensing for third-party material | [Third-party notices](THIRD_PARTY_NOTICES.md) |

Audit performed on 19 October 2026 with AI-assisted scientific computation. Agreement between calculations checks the software; it does not constitute new independent astronomical observations. No author endorsement is implied.
