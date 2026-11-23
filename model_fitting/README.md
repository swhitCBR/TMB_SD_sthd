# model_fitting

Scripts for the component-model fitting pipeline (HOR, TCJ, HOR_TCJ, HOR_CHP, TCJ_CHP). Run them all with `source("model_fitting/runscript.R")` from the repository root; see the top-level README for setup and the pinned package versions.

## Verification after the move from CVPAS_STH_app (2026-11-23)

**Purpose.** Check that the pipeline still produces the same results after being moved out of `CVPAS_STH_app` into this repository, with TMBhelp pinned to v1.0.0 (previously loaded from the sibling source folder) and the other packages locked by renv.

**What was run.** All 13 steps of `runscript.R` (data prep, template compile, five candidate-fit scripts, three confidence-set scripts) in this repository, in three stages: data prep and compile, glmmTMB fits, then the parallel TMB fits and confidence sets. R 4.5.1, TMBhelp 1.0.0, TMB 1.9.18, glmmTMB 1.1.14, MuMIn 1.48.19, 8 workers. Total about 55 minutes. Inputs were read from `tmp_data/` (copies of the files in `../CVPAS-steelhead/`).

**How it was checked.** The new files in `output/` were compared with the stored results in `CVPAS_STH_app/output/` (`all.equal`, tolerance 1e-6 unless stated), looking at AIC tables, confidence sets, coefficients and fit objects.

| Output | Result |
|---|---|
| HOR, TCJ, HOR_TCJ glmmTMB results | AIC tables equal; model AIC and fixed effects agree to about 1e-10 |
| `HOR_CHP_TMB_all_mods.rds` (1024 models) | nll and AIC differ by 4e-7 at most; only optimizer iteration counts differ (up to 5) |
| `HOR_CHP_TMB_mod_fits_ls.rds` | nll and AIC differ by about 2e-7, all fits converged; parameters differ by up to 3e-4 and random-effect estimates by up to 4e-5 (optimizer stopping point on a flat likelihood) |
| `HOR_CHP_S_conf_ls.rds`, `HOR_CHP_TMB_all_mods_unscl.rds` | AIC differs by 6e-8 or less; coefficients by up to 3e-4 |
| `TCJ_CHP_TMB_all_mods.rds` | AIC differs by about 1e-10; only `TCJ_CHP_comp_ls$XX_in` differs (scaled in the current `TCJ_CHP_comp()`, as noted on 2026-11-22) |

Conclusion: the pipeline gives the same model sets, AIC rankings and confidence sets as in `CVPAS_STH_app`; remaining differences are at optimizer-tolerance level.

**Not verified.**
- The data-prep intermediates in `data/` were not compared directly (R could not read `CVPAS_STH_app/data/` in this session). They are supported only indirectly, because every downstream result matches.
- `TCJ_CHP_S_conf_ls.rds` has no stored counterpart in `CVPAS_STH_app/output/`.
- `HOR_env_obs.rds` and `TCJ_env_obs.rds` were not compared (see below).
- The edits made after this run (the `local()` change below) have not been re-run.

**Known side effects seen in this run.**
- `HOR_env_obs.rds` and `TCJ_env_obs.rds` came out at 3.1 GB each (stored versions: 130 MB and 25 MB) because the cand-mods scripts save the global environment, which still held the 3.1 GB route workspace from the data-prep step. The data-prep scripts now run their loading inside `local()` so the input workspaces no longer stay in the global environment.
- glmmTMB reports that it was built against TMB 1.9.21 while 1.9.18 is installed; the glmmTMB results still match the stored ones.

## Changes made when the scripts were moved

- `library(TMBhelp)` replaces `devtools::load_all("../TMBhelp")`; `runscript.R` stops if the installed TMBhelp differs from `tmbhelp_version` in `config.R`.
- New step `compile_templates.R` compiles the two TMB templates in `TMB/` before any parallel fit.
- `config.R` settings no longer overwrite values defined before sourcing, and it defines `input_dirs` / `input_path()`: upstream inputs are looked for in `tmp_data/`, then `../CVPAS-steelhead/tmp_data`, then `../CVPAS-steelhead`.
- Repository-root checks and comments refer to `TMB_SD_sthd`.
