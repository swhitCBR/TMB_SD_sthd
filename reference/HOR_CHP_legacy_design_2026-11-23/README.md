# HOR_CHP results with the legacy design matrix (2026-11-23)

Frozen copies of the HOR_CHP model-fitting results produced by this repository **before** the design matrix was changed. Keep them to check that refactored code, or the downstream CVPASapp repository, reproduces the same numbers. Do not edit or regenerate these files; if the design changes, add a new dated folder next to this one.

## What "legacy design" means

The 1,024 candidate models were fit with the design matrix built by `HOR_CHP_comp()` (`model_fitting/HOR_CHP_surv_mods_data_prep.R`), which has two features that were later judged to be unintended:

- **`R_x_Temp` uses `Tclc.hor.3`.** The route x temperature interaction multiplies the route indicator by the scaled `Tclc.hor.3`, while the temperature main effect and `B_x_Temp` use `Tmsd.hor.7dadm`. The intended variable is `Tmsd.hor.7dadm` in both.
- **`WYT_drought` does not follow each row's `water.yr`.** `WYT_wet` equals the row's `water.yr == "wet"`, but `WYT_drought` is built from a release-group lookup (`tmp_tab[tmp_tab[,"row"],"col"]` in the release-group code). In these data it disagrees with `water.yr == "drought"` on 1,544 of 4,611 rows. Investigation (2026-11-24): each of the 19 release groups has exactly one `water.yr`, but the lookup assigns the wrong type to 6 groups (2014 S2, 2014 S3 and 2015 S1 are `drought` but were coded `dry`; 2016 S1, S2 and S3 are `dry` but were coded `drought`), 1,544 rows in all. The cause is that `which(WYT_ref_tab != 0, arr.ind = TRUE)` lists the non-zero cells column by column (wet groups, then dry, then drought), and `tmp_tab[tmp_tab[,"row"],"col"]` then indexes that list by group number as if it were in group order. It is only correct when the groups happen to be sorted by water-year type. `WYT_wet` is unaffected.

Because `WYT_` terms are in every candidate, a change to `WYT_drought` would move all 1,024 models. A change to `R_x_Temp` moves only the 512 models that include it.

## Files

| File | Contents |
|---|---|
| `HOR_CHP_comp_ls.rds` | Output of `HOR_CHP_comp()`: `allint_DMs`, `TMB_data_baseline`, `XX_in`, `XX_in_w_int_WYT` (scaled covariates) |
| `HOR_CHP_comp_ls_unscl.rds` | Same with `z_scale_vars = FALSE` |
| `HOR_CHP_TMB_mod_fits_ls.rds` | Raw fits for all 1,024 candidates (`OPT`, `nll`, `AIC`, `k`, `converged`, `RELGRP_err_v`, `wall_secs`) |
| `HOR_CHP_TMB_all_mods.rds` | Summary list: `HOR_CHP_comp_ls`, `AIC_DF_full`, `pt_estsDF`, `par_nm_ind_ls`, `allint_DMs`, `TMB_data_baseline`, `XX_in_w_int_WYT`, `scaled_vars`, `center`, `scale`. Its `XX_in_w_int_WYT` element is `NULL`; use `HOR_CHP_comp_ls$XX_in_w_int_WYT` |
| `HOR_CHP_S_conf_ls.rds` | Confidence set: `AIC_DF_d2` (13 models) and `S_coef_confset_DF` (199 rows) |
| `HOR_CHP_TMB_all_mods_unscl.rds` | Confidence-set models refit on unscaled covariates |
| `MANIFEST.csv` | Size and md5 of every file here |

## Design matrix layout (`XX_in_w_int_WYT`, 4,611 rows x 21 columns)

Column order, which the `allint_DMs` strings index by position:

1. `(Intercept)`, `WYT_wet`, `WYT_drought`, `route.facB`, `barrier.facTRUE`, `flength`, `Tmsd.hor.7dadm`, `log.VNS.hor.5`, `SWP.hor.5`, `CVP.hor.5`, `Qomt.hor.1net` (always in, characters 1-11 of each string)
2. `B_x_Temp`, `B_x_VNS`, `B_x_SWP`, `B_x_CVP`, `B_x_OMT`, `R_x_Temp`, `R_x_VNS`, `R_x_SWP`, `R_x_CVP`, `R_x_OMT` (toggled, characters 12-21)

`allint_DMs` is a character vector of 1,024 strings of 21 characters. The first interaction (`B_x_Temp`) varies fastest: string 1 has no interactions and string 2 has only `B_x_Temp`.

## How these were produced

| Item | Value |
|---|---|
| Date and machine | 2026-11-23 19:17-19:40 local, NERVA |
| Command | `source("model_fitting/runscript.R")` from the repository root, all 13 steps, 12 workers, 26.4 minutes |
| Git | HEAD `1047562bc84beffe21af917886e7c65b287bfd31` plus uncommitted edits (worker count in `config.R`, removal of the `Survival model/` subfolder from input paths) |
| R / packages | R 4.5.2, TMB 1.9.18, glmmTMB 1.1.14, MuMIn 1.48.19, TMBhelp 1.0.0. `rlang` 1.3.0 and `Rcpp` 1.1.2 differ from `renv.lock` |
| Inputs | `tmp_data/tmp_img4FULL.Rdata` (md5 `92649ca776bdd2ef5fde919a9891943a`), `tmp_data/Multiyear v2a AllVar v3.RData` (md5 `b5916dcaeb2becc02417b356e2c12649`) |
| Intermediates | `data/HOR_CHP_mod_dat_ls.RData` (md5 `a57048a22f235f663f815f66d2e4d390`), `data/HOR_CHP_mod_dat_sep_ls.rds` (md5 `7cf6b0ec35110ca56a35e890292696f0`) |

Rows: the 4,637 rows of `x.df` minus the 26 with a missing `det.pred.qomt.rms`. All 1,024 fits converged.

## Refactor verification log (2026-11-24)

The design matrix and candidate strings are now built by `R/HOR_CHP_designmat.R` (`model.matrix()` and `MuMIn::dredge()` from a global formula and `fixed` terms), with temporary adapters in `R/legacy_adapters.R`. `reference/verify_HOR_CHP_refactor.R` repeats the design checks; set `run_fits <- TRUE` before sourcing it to repeat the fits (about 9 minutes each, 12 workers).

| Stage | Change | Result against the files in this folder |
|---|---|---|
| A | Legacy mode: same `R_x_Temp` and `WYT_drought` as before | `HOR_CHP_comp_ALT()` and the formula-based design give a list `identical()` to `HOR_CHP_comp_ls.rds`. The unscaled `TMB_data_baseline` is identical too (the unscaled `XX_in_w_int_WYT` differs only in the order of its attributes). All 1,024 refits have nll, AIC and parameters exactly equal to `HOR_CHP_TMB_mod_fits_ls.rds` |
| B | `R_x_Temp` uses `Tmsd.hor.7dadm` | Only the `R_x_Temp` column changes. The 512 models without `R_x_Temp` are exactly unchanged. The 512 with it change (median absolute change in AIC 0.32, maximum 4.34). All 1,024 converged. The best model (candidate 833) is the same; the number of models within 2 AIC of the best goes from 13 to 14; Spearman correlation of AIC with the legacy values is 0.964 |
| C | `WYT_drought` from each row's `water.yr` | Applied in `../HOR_CHP_Rxtemp_Tmsd_WYTfixed_2026-11-24/`: all 1,024 models change (worse AIC by 11 to 22), the best model differs and the confidence set grows to 16 models |

The Stage A and B fit lists are in `output/refactor_check/` (git-ignored, not frozen here).

## Using these to verify another repository (e.g. CVPASapp)

```r
ref <- "reference/HOR_CHP_legacy_design_2026-11-23"
old_fits <- readRDS(file.path(ref, "HOR_CHP_TMB_mod_fits_ls.rds"))
old_comp <- readRDS(file.path(ref, "HOR_CHP_comp_ls.rds"))

# design matrix and candidate strings
identical(new_comp$allint_DMs, old_comp$allint_DMs)
all.equal(new_comp$XX_in_w_int_WYT, old_comp$XX_in_w_int_WYT, tolerance = 1e-8)

# fits: compare nll and AIC model by model
new_aic <- vapply(new_fits, function(x) x$AIC, 0); old_aic <- vapply(old_fits, function(x) x$AIC, 0)
max(abs(new_aic - old_aic))
```

Tolerances: `model_fitting/README.md` found differences of up to about 4e-7 in nll and AIC and up to 3e-4 in parameters between machines (optimizer stopping point on a flat likelihood). Same-machine reruns of the same code should agree much more closely.
