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

## Expanded candidate sets and the 2026-11-25 full run

**HOR_CHP candidate expansion.** `HOR_CHP_fixed_terms` in `R/HOR_CHP_designmat.R` was reduced to `c("WYT_", "flength")`, so the five environmental covariates (Tmsd, log.VNS, SWP, CVP, Qomt) and the `barrier.fac` / `route.fac` main effects are droppable, required only when one of their interactions is in the model (`dredge()`'s marginality). The candidate set grows from 1024 (all mains fixed, only the 10 interactions varying) to **3643 models** (32 + 2x243 + 3125, by barrier/route main-effect presence). Verified before fitting with `legacy_HOR_CHP_design(precomp, fixed = ...)` (`dredge(evaluate = FALSE)`), and by fitting two new-style candidates (dropped mains) through `TMBhelp::fit_nested_tmb_mods()`, both converging with SEs.

**Verification against the previous output.** The expanded fit (scratch run, same settings) reproduces all 1024 old models with **identical AICs** (max difference 0); 2619 models are new. The best model is unchanged and stays rank 1. The dAIC < 2 set grows from 16 to 22 models (a few with dropped mains enter). The previous 1024-model tables are kept next to the new ones as `HOR_CHP_TMB_all_mods_pre_expansion.rds` and `HOR_CHP_TMB_mod_fits_ls_pre_expansion.rds`.

**HOR_TCJ dredge conversion.** The 25 hand-written formulas were replaced by a global model + `dredge()` (582 candidates; models with both Year and WYT excluded because the WYT dummies are collinear with Year). `dev/verify_HOR_TCJ_dredge.R` confirms all 25 manual models appear in the dredged set with AIC differences below 1.1e-11, identical parameter counts and AIC ordering.

**Full pipeline run, 2026-11-25** (R 4.5.2, TMBhelp 1.0.0, 14 workers, `output/`):

| Step | Time (min) |
|---|---|
| `HOR_TCJ_separate_route_mods_data_prep.R` | 0.8 |
| `HOR_TCJ_surv_mods_data_prep.R` | 2.0 |
| `HOR_CHP_surv_mods_data_prep.R` | 0.6 |
| `TCJ_CHP_surv_mods_data_prep.R` | 0.3 |
| `compile_templates.R` | 0 |
| `HOR_cand_mods_run.R` | 1.7 |
| `TCJ_cand_mods_run.R` | 1.2 |
| `HOR_TCJ_cand_mods_run.R` (dredge, 582 models x 2) | 5.4 |
| `HOR_CHP_cand_mods_run.R` (3643 models) | 23.3 |
| `TCJ_CHP_cand_mods_run_init.R` | 0.9 |
| `TCJ_CHP_conf_comp.R` | 0.3 |
| `HOR_CHP_conf_comp.R` (22 d2 refits) | 1.8 |
| `HOR_CHP_d2_AIC_unscl_mod_ls.R` (22 refits, all converged) | 4.1 |
| Total | **42.2** |

Caveat found in this run (resolved 2026-11-25): `HOR_CHP_cand_mods_run.R`, `TCJ_CHP_cand_mods_run_init.R` and `HOR_CHP_d2_AIC_unscl_mod_ls.R` used to write their `.rds` results only `if (!file.exists(...))`, so a rerun over an existing `output/` kept stale tables (`output/HOR_CHP_TMB_all_mods_unscl.rds` still held the 16-model pre-expansion confidence set). The guards were removed; the files are now always overwritten, and `output/HOR_CHP_TMB_all_mods_unscl.rds` was replaced with the 22-model version.

## HOR_TCJ without Year (`hor_tcj_include_year`, 2026-11-25)

`pipeline_config$hor_tcj_include_year` (default `TRUE`) controls whether `Year` is in the HOR_TCJ survival global model. With `FALSE` the dredge has 388 candidates (the Year-free subset of the 582; the AICs are the same because the models are the same, only dAIC and weights are recomputed) and WYT models compete only with each other. The outputs then get the suffix `_noYear` (`HOR_TCJ_all_mods_noYear.rds`, `HOR_TCJ_d2_mods_noYear.rds`, `HOR_TCJ_conf_sets_noYear.rds`) so the standard files are not overwritten. From the stored table the Year-free d2 set has 13 models (5 with WYT); see `CVPAS_STH_app/dev/model_set_changes_TMB_SD_sthd_README.md`. The option has not been run: set it with `source("model_fitting/config.R"); pipeline_config$hor_tcj_include_year <- FALSE` before sourcing `runscript.R` (use a scratch `out_dir` first).
