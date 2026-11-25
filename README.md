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
source("model_fitting/runscript.R")   # about 30 minutes (workers: n_workers in config.R, default 12)
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

## FAQ

### How do I update dependencies with renv?

> [!NOTE]
> 1. **Start clean.** Restart R so no packages are loaded, and run `renv::status()` to confirm the library matches `renv.lock`. Replacing a loaded package can leave the library in a broken state.
> 2. **Update the packages.**
>    - `renv::update()` updates everything to the latest available versions.
>    - `renv::update("glmmTMB")` or `renv::install("TMB@1.9.19")` updates one package, or installs a specific version.
>    - For TMBhelp, install the new tarball with `renv::install("path/to/TMBhelp_X.Y.Z.tar.gz")` and set `tmbhelp_version` in `model_fitting/config.R`.
> 3. **Restart R, then test.** Run the pipeline into a scratch folder (e.g. `out_dir <- "C:/scratch/out"`) and compare the results with the current `output/`, using the 1e-6 tolerance from `model_fitting/README.md`.
> 4. **Record or revert.**
>    - If the results are acceptable, run `renv::snapshot()` and commit `renv.lock`. Check that `renv/.gitignore` still contains `!local/` and that `renv/local/TMBhelp_*.tar.gz` is still tracked.
>    - If not, `renv::restore()` returns the library to the locked versions.
>
> **Cautions for this project**
> - `renv::restore()` fails on the TMBhelp record. Use `renv::restore(exclude = "TMBhelp")`, then install the tarball separately.
> - TMB, glmmTMB and Rcpp are tightly coupled. Update them together and expect a rebuild from source, which needs Rtools.
> - `renv::history()` shows earlier lockfile versions in git, and `renv::revert(commit = "<hash>")` goes back to one of them.
