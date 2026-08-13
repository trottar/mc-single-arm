# 3He SF fit tables

Place the CSV tables required by `GETSF_F1F2fit` in this directory:

- `Table_3He_F1F2_SF1.csv`
- `Table_3He_F1F2_SF2.csv`
- `Table_3He_F1F2_SF3.csv`
- `Table_3He_F1F2_SF4.csv`
- `Table_3He_F1F2_SF5.csv`

`GETSF_F1F2fit` reads these repo-local paths when `mc_single_arm` runs from
`src/`.  If the files are absent, it returns `STAT=.false.`.


## Archive preparation

Normal `make` does not require any table archive.  Before using
`SF model flag = 1`, run `make prepare_sf_tables` in `src/`.  The target checks
all five CSVs and, if needed, extracts one of:

- `interp/F1F221_3He_XZ_20250828_tables.tar.gz`
- `interp/F1F221_3He_XZ_20250828_tables2.tar.gz`
- `interp/F1F221_3He_XZ_20250828_tables.tar`
- `interp/F1F221_3He_XZ_20250828_tables.tar.*` (split tar parts)

The supplied `.tar.gz` contains an `outfiles_tmp/` prefix; the target strips it
and places the CSVs directly in `interp/sf_tables/`.
