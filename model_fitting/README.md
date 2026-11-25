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

**Data-prep intermediates (`data/`).** Stage 1 was run twice, before and after wrapping the data-prep loading in `local()`. The six regenerated files were `identical()` between the two runs. They were also compared with copies of the stored `CVPAS_STH_app/data/` files:

| File | Result |
|---|---|
| `HOR_CHP_mod_dat_ls.RData`, `HOR_TCJ_survmod_data.RData`, `TCJ_CHP_mod_dat_ls.RData` | identical to the stored files |
| `route_dfs_HOR_TCJ.rds`: `hor.df`, `hor.stn.df`, `d.df` | equal to the stored copy (1e-6) |
| `route_dfs_HOR_TCJ.rds`: `tcj.df`, `tcj.stn.df` | differ from the stored copy; the stored copy is stale (see below) |

The stored `data/` has no `route_dfs_HOR_TCJ.rds`. Its two route files (`route_dfs_HOR_CHP.rds` and `HOR_TCJ_route_mod_data.RData`, which holds an object named `route_dfs_HOR_CHP`) are the same older copy, in which `tcj.df` and `tcj.stn.df` are duplicates of the HOR frames (4774 rows). The regenerated TCJ frames have 1937 rows. The stored TCJ models were fit on 1937 observations, as were the new ones, and they agree to about 1e-10, so the regenerated TCJ frames are the ones that produced the stored TCJ results and the old copy is unused. The regenerated `HOR_TCJ_route_mod_data.RData` holds `route_dfs_HOR_TCJ`, not `route_dfs_HOR_CHP`.

**Not verified.**
- `HOR_CHP_mod_dat_sep_ls.rds` has no stored copy to compare with.
- `TCJ_CHP_S_conf_ls.rds` has no stored counterpart in `CVPAS_STH_app/output/`.
- `HOR_env_obs.rds` and `TCJ_env_obs.rds` were not compared (see below).
- The candidate-fit and confidence-set steps were run before the `local()` change and were not re-run afterwards; only stage 1 was re-run.

**Known side effects seen in this run.**
- `HOR_env_obs.rds` and `TCJ_env_obs.rds` came out at 3.1 GB each (stored versions: 130 MB and 25 MB) because the cand-mods scripts save the global environment, which still held the 3.1 GB route workspace from the data-prep step. The data-prep scripts now run their loading inside `local()` so the input workspaces no longer stay in the global environment.
- glmmTMB reports that it was built against TMB 1.9.21 while 1.9.18 is installed; the glmmTMB results still match the stored ones.

## Changes made when the scripts were moved

- `library(TMBhelp)` replaces `devtools::load_all("../TMBhelp")`; `runscript.R` stops if the installed TMBhelp differs from `tmbhelp_version` in `config.R`.
- New step `compile_templates.R` compiles the two TMB templates in `TMB/` before any parallel fit.
- `config.R` settings no longer overwrite values defined before sourcing, and it defines `input_dirs` / `input_path()`: upstream inputs are looked for in `tmp_data/`, then `../CVPAS-steelhead/tmp_data`, then `../CVPAS-steelhead`.
- Repository-root checks and comments refer to `TMB_SD_sthd`.

## Design matrix code and settings (2026-11-24)

The HOR_CHP and TCJ_CHP candidate fits now build their design matrices and candidate model strings (`allint_DMs`) with the functions in `R/` instead of `HOR_CHP_comp()` and `TCJ_CHP_comp()`. The old functions remain in the `*_surv_mods_data_prep.R` scripts, and their header comments still describe the design they build (the legacy design below).

| File in `R/` | Contents |
|---|---|
| `HOR_CHP_comp_ALT.R`, `TCJ_CHP_comp_ALT.R` | `*_precomp()` loads the workspace and builds the scaled covariates; `*_modelprep()` builds `TMB_data_baseline` and `allint_DMs`; `*_comp_ALT()` runs both and returns the same list as the old `*_comp()` |
| `HOR_CHP_designmat.R`, `TCJ_CHP_designmat.R` | The global model formula and the `fixed` terms; `HOR_CHP_designmat()` builds the design matrix with `model.matrix()` and the candidate strings with `MuMIn::dredge()` (used for both surveys) |
| `legacy_adapters.R` | Temporary code that relabels columns (`B_x_Temp`, `R_x_SWP`, ...), applies the settings below and compares designs (`compare_designs()`); delete it once the new column names are adopted |

Settings in `config.R` (define before sourcing `runscript.R` to override):

| Setting | Default | Meaning |
|---|---|---|
| `HOR_CHP_r_x_temp` | `"Tmsd.hor.7dadm"` | Temperature variable in the HOR_CHP route x temperature interaction (`R_x_Temp`). `"Tclc.hor.3"` is the legacy choice |
| `CHP_wyt` | `"water.yr"` | Source of the `WYT_wet` and `WYT_drought` dummies in HOR_CHP and TCJ_CHP. `"legacy"` is the release-group lookup, which assigned the wrong water-year type to 6 of 19 (HOR_CHP) and 6 of 18 (TCJ_CHP) release groups because of an indexing error |

Used by `HOR_CHP_cand_mods_run.R`, `HOR_CHP_d2_AIC_unscl_mod_ls.R`, `TCJ_CHP_cand_mods_run_init.R` and `TCJ_CHP_conf_comp.R`. To reproduce the results of the verification run above, set `HOR_CHP_r_x_temp <- "Tclc.hor.3"` and `CHP_wyt <- "legacy"` before sourcing.

Effect of the new defaults on the fits (details in the `README.md` of each `reference/` folder):

| Survey | Result |
|---|---|
| HOR_CHP, `R_x_Temp` only | The 512 models without `R_x_Temp` are unchanged; the other 512 change (median absolute change in AIC 0.3). Same best model |
| HOR_CHP, plus `WYT_drought` | All 1,024 models change (AIC worse by 11 to 22); different best model; confidence set 16 models |
| TCJ_CHP, `WYT_drought` | AIC changes by 6 on average; same best model and same 5-model confidence set |

Frozen results are in `reference/` (tracked in git, unlike `output/`), each folder with a `README.md` and `MANIFEST.csv`:

- `HOR_CHP_legacy_design_2026-11-23/`, `TCJ_CHP_legacy_design_2026-11-23/`: before the change (legacy design)
- `HOR_CHP_Rxtemp_Tmsd_2026-11-24/`: `R_x_Temp` from `Tmsd.hor.7dadm`, legacy `WYT_`
- `HOR_CHP_Rxtemp_Tmsd_WYTfixed_2026-11-24/`, `TCJ_CHP_WYTfixed_2026-11-24/`: current defaults

`reference/verify_HOR_CHP_refactor.R` and `reference/verify_TCJ_CHP_refactor.R` check that the new code in legacy mode reproduces the frozen legacy design exactly (set `run_fits <- TRUE` to refit and compare the AIC values).
