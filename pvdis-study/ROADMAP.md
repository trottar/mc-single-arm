# PVDIS study roadmap: inelastic/resonance A_PV first

## Primary goal - user correction, October 25, 2026

The urgent scientific deliverable is **inelastic/resonance A_PV**, with preliminary MOLLER workshop material by Monday afternoon October 26 and more complete results by Tuesday October 27. Broader PAC optimization is secondary. Notebook auditing, reporting infrastructure, and coupling ellipses do not substitute for this physics artifact.

The original four-point benchmark is only a first pass. It is not completion of the resonance goal. Missing external fit files do not block the next resonance-theory calculations.

## Active order of work

1. Resolve and document the electromagnetic-to-gamma-Z normalization, axial F3 assumptions, longitudinal ratio, and resonance isospin limits.
2. Quantify those assumptions against the published deuteron resonance data. Keep arbitrary response variations distinct from calibrated uncertainties.
3. Produce fresh 11 GeV A_PV(W,Q2) reference curves and representative angle/momentum kinematics, with machine-readable diagnostics and source provenance.
4. Add credible statistical projections only after a verified LD2 rate normalization and an explicit 30-day schedule.
5. Produce concise workshop-ready inelastic/resonance figures with assumptions and outstanding physics gates visible.
6. Only then advance the broader PAC C1/C2/HT/CSV ellipse and setting optimization.

## Completion criteria for the urgent deliverable

- A documented gamma-Z prescription and its unresolved physical inputs.
- Published-data comparison with finite-Q2 kinematics and identified response sensitivities.
- Fresh 11 GeV curves and representative settings, not substituted historical figures.
- Statistical reach only where rates and exposure are verified; otherwise explicitly unavailable.
- A compact scientific result in the cumulative PDF, backed by CSV, per-point JSON, exact source/table identities and execution logs.

## Current execution checkpoint

- Executed: phase-two fresh 11 GeV scan (896 points), four benchmark means, 10,800 controlled response cases, 40 kinematic landmarks; source/raw reconstruction verified. Sealed study: `resonance_studies/20261025T193128315988Z`.
- Documented and quantified: vector/axial separation, longitudinal sensitivity, and a pure-isospin component-limit discrepancy. Scenarios remain assumptions, not fitted or calibrated uncertainties.
- Still open: scientifically justified resonance/background gamma-Z construction and axial amplitudes, acceptance/radiative consistency, rate normalization and statistical reach, workshop interpretation review.
- Broader coupling-ellipse optimization remains secondary. Do not describe this checkpoint as completion of the urgent goal.

## Scope and authorization

- Broader PAC study; MOLLER workshop is its first presentation milestone.
- Total beam-time budget: **30 PAC days**, shared by simultaneous HMS/SHMS measurements; sequential configurations partition that total.
- Remote `pvdis` was created from `master` commit `2a27ddb7442ab2148e47fe4f3e998cf1be0a9fa8`.
- Calculations and review artifacts were prepared locally. The user subsequently explicitly authorized committing and pushing the study and reproducibility bundle to pvdis for chat analysis. No merge to master, farm execution, or hosted CI spending is part of that authorization.
- Uploaded references, `sources/`, and the existing `mc-single-arm/` workspace are read-only in this task.

## Exact writable scope

Study-owned files beneath `pvdis-study/`: this roadmap, `CURRENT.md`, `README.md`, `DEPENDENCIES.md`, `data/`, `vendor/`, `src/`, `tools/`, `build/`, `results/`, immutable `runs/`, `notebook_reviews/`, and `resonance_studies/` snapshots, and `output/pdf/`. Original references and historical snapshots remain unchanged. The user's explicit memory request is recorded only as an additive note in the designated internal-memory extensions folder.

## Historical first-benchmark work (completed locally, not the full goal)

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
