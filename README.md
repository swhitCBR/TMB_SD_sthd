# TMB_SD_sthd

Model-fitting pipeline for the CVPAS steelhead component models (HOR and TCJ routing, HOR_TCJ, HOR_CHP and TCJ_CHP survival), with the TMB templates it needs. Split out of `CVPAS_STH_app` so the pipeline is self-contained and pinned to specific package versions.

## Layout

| Path | Contents |
|---|---|
| `model_fitting/` | Pipeline scripts; `runscript.R` runs them in order, `config.R` holds settings |
| `TMB/` | TMB templates (`HOR_CHP_global.cpp`, `HOR_CHP_global_pred.cpp`); compiled on first run, binaries are git-ignored |
| `tmp_data/` | Local copies of upstream input files (git-ignored) |
| `data/` | Intermediates written by the data-prep steps (git-ignored) |
| `output/` | Pipeline results (git-ignored) |
| `dAIC_threshold.R` | Shared dAIC threshold sourced by the fitting scripts |
| `renv.lock`, `renv/` | Locked package versions; `renv/local/` holds the TMBhelp tarball |

## Running

From the repository root (open `TMB_SD_sthd.Rproj`; renv activates automatically):

```r
renv::restore()                       # first time on a new machine
source("model_fitting/runscript.R")   # about 30 minutes with 8 workers
```

Settings (`out_dir`, `n_workers`, `skip`, `tmbhelp_version`) can be defined before sourcing to override the defaults in `model_fitting/config.R`. To test without overwriting results, set `out_dir` to a scratch folder first (see the header of `runscript.R`).

Inputs are currently read from `../CVPAS-steelhead/` (the repository must sit next to this one). Rtools is needed to compile the TMB templates on Windows.

## Pinned versions

- **TMBhelp 1.0.0** (https://github.com/swhitCBR/TMBhelp, tag `v1.0.0`), installed from `renv/local/TMBhelp_1.0.0.tar.gz`. `runscript.R` stops if the installed version differs from `tmbhelp_version` in `config.R`.
- All other packages are locked in `renv.lock` (R 4.5.1).

## Updating

1. Install the new TMBhelp release, e.g. `renv::install("path/to/TMBhelp_X.Y.Z.tar.gz")`, and/or update other packages with `renv::update()`.
2. Set `tmbhelp_version` in `model_fitting/config.R`.
3. Run the pipeline into a scratch `out_dir` and compare with the previous results.
4. If the results are acceptable, run `renv::snapshot()` and commit `renv.lock`; otherwise `renv::restore()` returns to the locked versions.
