# Notebook review

Source SHA-256: afc987aefb602b4310b0b8e71aa3274ddb634ac37d522d912a736117a79e8193

Cell numbers are zero-based. Full notebook execution and Hall C fits were not performed.

## N01: Joint confidence threshold overwritten

Cells: [14, 15]

delta = chi2.ppf(cl, df=2) is immediately followed by delta = 1.0 in both fit functions.

The displayed 68.3% ellipse encloses 39.35% for a two-dimensional Gaussian. A joint 68.3% contour needs delta=2.2977 and axes 1.5158 times larger.

Expose delta-chi2/coverage explicitly and retain separate one-parameter errors and joint contours.

## N02: Hall C time scaling is inverted and applied to every setting

Cells: [3, 8]

hallc_stat_loq=sqrt(split/0.5), hallc_stat_hiq=sqrt((1-split)/0.5); cell 8 applies loq to all HMS and SHMS files.

At split=0.2 the notebook multiplies errors by 0.63246 instead of 1.58114 for reduced low-angle exposure relative to the same-total 50:50 baseline. It does not enforce the user's 30-day total budget.

Use sqrt(reference_time/allocated_time) from explicit setting-pair exposure records; reject zero exposure and duplicate time allocations.

## N03: Amplitude correlation fraction used as variance fraction

Cells: [4, 14, 15]

solid11_corr=sqrt(0.004^2+0.002^2)/solid11_syst, but fit covariance uses syst_corr=syst*sqrt(f_corr).

If the listed 0.4% and 0.2% sources are the intended correlated components, f_corr should be their squared sum divided by total systematic variance. Off-diagonal variance is inflated by about 27.5% for the nominal SoLID input.

Define correlated source amplitudes directly or use the squared amplitude fraction consistently; document the intended decomposition.

## N04: Undefined setting groups prevent a clean sequential run

Cells: [22, 28]

data_sets_hms1, data_sets_shms1, data_sets_hms2, data_sets_shms2 are used without assignment; result_11gev_hallc1 and result_11gev_hallc1_all are also used without definition.

Supplying the missing files alone would not make the exported notebook execute top to bottom.

Define the setting groups from an explicit schedule and remove dependencies on hidden interactive state.

## N05: Display recentering overwrites fitted parameter records

Cells: [18, 19, 20, 23, 30]

Multiple result['params'] assignments replace the fitted centers; the four-parameter result is replaced by a two-element display vector in cell 23.

Returned params can disagree with fitted, cov_full, chi2, and residuals. Shifted ellipses illustrate precision, not measured best fits.

Preserve fitted_params and covariances; store display_center separately, with the shift and purpose recorded.

## N06: Common systematics stop at dataset boundaries

Cells: [14, 15]

Covariance is block diagonal between dataset dictionaries, with one same-sign correlated source inside each block.

Polarimetry or other common sources across HMS/SHMS settings are not represented unless all relevant rows are put into one explicitly shared source block.

Represent each physical systematic as a source vector spanning every affected measurement; validate signs and overlap.

## N07: Hadronic templates supplied, but application and authority need review

Cells: [5, 15]

A=(c1*X+c2*Y)*(1+beta_HT/((1-x)^3*Q2)+beta_CSV*x^2); beta applies only where the second coefficient is nonzero.

beta_HT has units GeV^2 and beta_CSV is dimensionless. Both terms multiply the full asymmetry. The coefficient mask is not an explicit dataset-physics label; no priors or bounds are supplied.

Retain this exact ansatz as notebook provenance; validate its scientific interpretation and use an explicit apply_hadronic flag before the PAC fit.

## N08: Derived external constraints treated as independent priors

Cells: [5, 18]

Elastic+Cs, the Qweak-derived X constraint, and projected P2 X constraint are separate independent scalar datasets with placeholder x=0.99,Q2=10000.

Independence of the derived Qweak constraint from the elastic+Cs input is not established by the notebook. P2 is a projection, not a measurement.

Recover the upstream likelihood/provenance and shared information; do not claim an independently combined world fit until this is resolved.

## N09: External tables missing and parsers silently skip rows

Cells: [6, 8, 11, 12]

28 Hall C tables, four HRS tables, solid11.dat, and pvdis_11gev_CJ15nlo.csv are referenced but not embedded. Parsers ignore wrong column counts and nonnumeric lines.

Full projections cannot be reproduced from this notebook alone; malformed numerical rows can disappear without a diagnostic.

Supply exact tables and hashes; use explicit schemas, row rejection logs, unit checks, and finite-value/rate validation.

## N10: Labeling and legacy contour geometry need cleanup

Cells: [8, 30, 36]

SHMS datasets are named HMS; the Lite legend can use nominal syst; PCA/max-extent approximates historical contour point sets.

Labels can misidentify datasets, and an enclosing display ellipse is not an inferred likelihood/covariance contour.

Correct labels and retain historical contour points with their original confidence convention and provenance.