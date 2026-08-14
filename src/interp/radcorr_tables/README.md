# 3He radiative-weight tables

`get_radcorr_3he.f` uses tables extracted from
`Newfit_20260710_fullxquad_15angles.tar.gz`.  Keep that archive in
`src/interp/radcorr_tables/` and run:

```bash
cd src
make prepare_radcorr_tables
```

The legacy `src/interp/` archive location is also accepted.

The archive expands to:

```text
Newfit_20260710_fullxquad_15angles/
  SF1_G1F1cmplt_QE95/ ... SF5_G1F1cmplt_QE95/
  radiated_data_26.0deg_short.csv ... radiated_data_33.0deg_short.csv
```

Each selected SF model supplies 15 angles in 0.5-degree increments.  The CSV
columns used are `Ebeam`, `nu`, `XSborn_unp`, and `XSrad_unp`.  `Ebeam` and
`nu` are MeV (`Ebeam=10380 MeV`); the MC converts `nu_model` from GeV once at
the reader boundary.

The stored quantity is:

```text
rad_weight_factor = XSrad_unp / XSborn_unp
```

It is bilinearly interpolated in `(theta, nu)`, with no extrapolation.  The
finite upper endpoint is 9910 MeV; the lower endpoint depends on angle.  A
lookup between two angles succeeds only over their overlapping `nu` domain:
the effective lower/upper limits are the maximum/minimum of their individual
finite limits.

Every finite grid is validated during initialization as 1 MeV spaced.  The
production reader therefore uses direct `nu` indexing after exact grid-point
and global-endpoint checks; only the selected SF model's 15 tables are loaded.

Historical A1n documentation may use `radcorr` for the inverse ratio
`XSborn_unp / XSrad_unp`.  This table reader deliberately uses the ratio above
to turn the MC's Born F1/F2 weight into a radiated weight.

For offline validation of all models and angles, run:

```bash
python3 util/check_radcorr_tables.py \
  src/interp/radcorr_tables/Newfit_20260710_fullxquad_15angles
cd util/radcorr && make
cd ../../src && ../util/radcorr/test_get_radcorr_3he
```

From `src/`, the integrated targets are `make test_radcorr_tables`,
`make test_radcorr_reader`, `make test_radcorr_born`, and `make test_radcorr`.
Use `make clean_radcorr_tables` only when a fresh archive extraction is needed;
ordinary `make clean` intentionally preserves local table data.
