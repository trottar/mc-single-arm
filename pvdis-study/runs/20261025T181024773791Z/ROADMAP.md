# PVDIS PAC study: first resonance benchmark

## Scope and authorization

- Broader PAC study; MOLLER workshop is its first presentation milestone.
- Total beam-time budget: **30 PAC days**, shared by simultaneous HMS/SHMS measurements; sequential configurations partition that total.
- Remote `pvdis` was created from `master` commit `2a27ddb7442ab2148e47fe4f3e998cf1be0a9fa8`.
- This work creates local review material only. No commits, pushes, or hosted CI runs.
- Uploaded references, `sources/`, and the existing `mc-single-arm/` workspace are read-only in this task.

## Exact writable scope

Study-owned files beneath `pvdis-study/`: this roadmap, `README.md`, `DEPENDENCIES.md`, `data/`, `vendor/`, `src/`, `tools/`, `build/`, `results/`, immutable `runs/` snapshots, and `output/pdf/`. The user's cumulative-report request authorizes updating these study files and adding the PDF and diagnostic collection workflow.

## Ordered work

1. Pin the public model source and audit external dependencies.
2. Transcribe Table I of arXiv:1304.7741v2 with independent column checks.
3. Execute unchanged F1F2IN21 for free proton/neutron. Use deuteron electromagnetic smearing only with verified table provenance and explicit initialization.
4. Construct the light-valence, charge-symmetric gamma-Z rotation with documented coupling conventions, finite-Q2 factors, and model assumptions. Do not call it a complete weak resonance model.
5. Check algebraic limits, normalization, repeatability, and invalid-input rejection. Compare model values at published mean kinematics with data; no acceptance or radiative unfolding is inferred.
6. Produce a CSV, plots, and a short technical result with limitations. Identify dependencies still needed for PDF/rate projections and HT/CSV fits.

## Completed first pass

All six first-benchmark steps above were completed locally. The pinned F1F2IN21 code and three externally sourced CD-Bonn table blobs passed identity checks. Four mean kinematics plus 804 W-grid points executed with bounds/FP checks and byte-identical fresh-process repetition. Results, analytic checks, caveats, and missing dependencies are recorded in `README.md`, `DEPENDENCIES.md`, and `results/`. These are model comparisons, not detector/rate/weak-resonance validation.

## Gates

- No fitted coupling ellipse without explicit likelihood, nuisance templates, and correlations.
- No projected statistical errors from unverified rate normalization.
- No physical LD2 detector claim from prior nominal paddle software checks.
- Do not use an unsmeared weak numerator over a smeared electromagnetic denominator as the nominal deuteron prediction.
- No invented PDF grids or replacement definitions for the absent original Fortran helpers.

## Cumulative evidence collection

Each collected run must validate the source/model/input identities, raw repeat, numerical CSV contents, and diagnostic gates before it can enter the run index. Preserve its inputs, source, logs, data, plots, units, and machine diagnostics in a new snapshot directory. Build the cumulative PDF from indexed snapshots; never overwrite an older snapshot. Keep analysis refresh (using retained raw model output) distinct from a fresh compiled model run.
