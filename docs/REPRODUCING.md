# Reproduction

Use Python 3.12 and install `requirements.txt` in a fresh virtual environment. No data from a parent workspace, API keys, paid compute, or vendored binary library is required. The tested Linux galpy wheel is version 1.12.0, SHA-256 `55c5c5ea0a53fa954dab1ccc3ba2cc3fdcb5d68c80474ba4b8630c480ecc86a0`.

## Immediate checks

```bash
python3 scripts/check_coordinates.py
.venv/bin/python scripts/verify.py
```

The first command uses only standard-library arithmetic and the retained source rows. The second verifies the checksum inventory, coordinate conversion, archived encounter counts, paired input equality, event thresholds, timestep comparisons, recorded independent ODE checks, figure-source hash, and the retained one-event uncertainty counterexample. It checks saved evidence; it does not replace fresh numerical reproduction.

## Fresh calculations

```bash
.venv/bin/python scripts/reproduce.py --quick
.venv/bin/python scripts/reproduce.py --full
```

Both commands make an isolated copy under `reproduced/<UTC timestamp>/`, preserve the archived files, and keep a log for every subprocess. Open `REPRODUCTION.json` for exact package versions, command outcomes, and comparisons with the archived results. The quick run uses 512 shared histories per coordinate branch and repeats the vector extraction, inversion, inversion controls, and independent ODE checks of selected archived histories. The full run recalculates the four 32,768-history Gaussian/uniform coordinate pairs, timestep controls, 8,192-history conditional-peculiar alternatives, and frame/uncertainty/encounter-summary checks as well.

Raw vector extraction is reproducible from the original source asset with Poppler:

```bash
pdftocairo -svg audit/sources/j0737_oc_orbits.pdf figure2-reextracted.svg
```

Different Poppler versions can serialize equivalent paths differently, so exact SVG checksums need not match. The archived SVG is included as the fixed input for the numerical audit. `figure_vector.py` explicitly calibrates each axis using five tick positions and labels; the code does not infer a position from a screenshot by eye.

## Inspecting arrays

Load NPZ files with `numpy.load(path, allow_pickle=False)`. Each baseline row stores the two encounter flags, minimum separations and lookback times, binary distance and heliocentric radial velocity, and weights used for explicitly labeled sensitivity calculations. Use the direct-prior unweighted counts for the quoted binomial bounds. A weighted zero estimate must not be assigned the unweighted bound.

The retained Gaussian pair uses exactly the same distance and radial draws. Proper-motion and cluster draws are reconstructed from the same recorded RNG seed and equations. The independent ODE code reconstructs those draws exactly and checks selected events, non-events, boundary-near histories, and radial-velocity extremes.

## Portability changes

The exported numerical code changes only directory/import/output paths and removes the dependency on an environment-specific vendor folder. Equations, physical constants, seeds, criteria, and saved result arrays are retained. The numerical code and result provenance are listed in `provenance/EXPORT.json`. The standard-library coordinate checker, verification script, reproduction wrapper, and documentation were added for this standalone repository.
