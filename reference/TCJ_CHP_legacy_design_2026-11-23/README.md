# TCJ_CHP results with the legacy design matrix (2026-11-23)

Frozen copies of the TCJ_CHP model-fitting results produced by this repository before the design matrix code was refactored. Keep them to check that refactored code, or the downstream CVPASapp repository, reproduces the same numbers. Do not edit or regenerate these files; add a new dated folder if the design changes.

TCJ_CHP fits use the same TMB template as HOR_CHP (`HOR_CHP_global`). See also `../HOR_CHP_legacy_design_2026-11-23/`.

## Design matrix layout (`XX_in_w_int_WYT`, 1,914 rows x 13 columns)

Column order, which the `allint_DMs` strings index by position:

1. `(Intercept)`, `WYT_wet`, `WYT_drought`, `TCJroute.facF` (always in, characters 1-4 of each string)
2. `flength`, `Qomt.tcj.1net`, `log.VNS.tcj`, `CVP.tcj.4`, `SWP.tcj.4`, `Tmsd.tcj.7dadm`, `R_x_SWP`, `R_x_CVP`, `R_x_OMT` (toggled, characters 5-13)

`allint_DMs` has 216 strings of 13 characters: the 2^9 = 512 combinations of the nine toggled terms, minus the 296 with a route interaction but not its covariate main effect. The first toggled term (`flength`) varies fastest. `TCJroute.facF` is 1 for TRN fish and 0 for MAC fish. `log.VNS.tcj` is `log(VNS.tcj)`, and the scaled design replaces missing `Qomt.tcj.1net` by 0 before scaling.

## Files

| File | Contents |
|---|---|
| `TCJ_CHP_comp_ls.rds` | Output of `TCJ_CHP_comp()` (scaled): `allint_DMs`, `WYT_ref_tab`, `rel_grps_df_tmp`, `XX_in_add`, `XX_in_w_int_WYT`, `TMB_data_baseline`, `XX_in` |
| `TCJ_CHP_comp_ls_unscl.rds` | Same with `z_scale_vars = FALSE`. Its design matrix has 12 `NA` values (6 rows, in `Qomt.tcj.1net` and `R_x_OMT`), so it cannot be fit as is |
| `TCJ_CHP_TMB_mod_fits_ls.rds` | Raw fits for all 216 candidates (`OPT`, `nll`, `AIC`, `k`, `converged`, `RELGRP_err_v`, `wall_secs`) |
| `TCJ_CHP_TMB_all_mods.rds` | Summary list: `TCJ_CHP_comp_ls`, `AIC_DF_full`, `pt_estsDF`, `par_nm_ind_ls`, `allint_DMs`, `TMB_data_baseline`, `XX_in_w_int_WYT`. The last element is `NULL`; use `TCJ_CHP_comp_ls$XX_in_w_int_WYT` |
| `TCJ_CHP_S_conf_ls.rds` | Confidence set: `AIC_DF_d2` (5 models) and `S_coef_confset_DF` (44 rows) |
| `MANIFEST.csv` | Size and md5 of every file here |

## Known quirks of this design

- **`WYT_drought` is mislabelled for 6 of the 18 release groups.** Each release group has one water-year type, but the lookup (`tmp_tab[tmp_tab[,"row"],"col"]` in the release-group code) indexes the list of non-zero table cells by group number, and that list is ordered column by column (wet, dry, drought). Groups 2014 S3, 2015 S1 and 2015 S2 are `drought` but were coded `dry`; 2016 S1, S2 and S3 are `dry` but were coded `drought`. This affects 520 of 1,914 rows and every candidate model, because the `WYT_` terms are always in. `WYT_wet` is unaffected. The HOR_CHP design has the same bug (see its README).
- `n_relgrps` is hard-coded as 18 (correct for these data: 18 release groups).
- `col_S` in `TMB_data_baseline` is `0` (a single zero), a consequence of `S_col_ids <- nrow(XX_in_w_int_WYT)` in the original code. It was kept as is.

## How these were produced

| Item | Value |
|---|---|
| Date and machine | 2026-11-23 19:34-19:36 local, NERVA |
| Command | `source("model_fitting/runscript.R")` from the repository root, all 13 steps, 12 workers |
| Git | HEAD `1047562bc84beffe21af917886e7c65b287bfd31` plus uncommitted edits (worker count in `config.R`, removal of the `Survival model/` subfolder from input paths) |
| R / packages | R 4.5.2, TMB 1.9.18, glmmTMB 1.1.14, MuMIn 1.48.19, TMBhelp 1.0.0. `rlang` 1.3.0 and `Rcpp` 1.1.2 differ from `renv.lock` |
| Inputs | `tmp_data/tmp_img4FULL.Rdata` (md5 `92649ca776bdd2ef5fde919a9891943a`), `tmp_data/lvec.tcj_ls.rds` (md5 `7733fdd1458c4c844286da0097784192`) |
| Intermediate | `data/TCJ_CHP_mod_dat_ls.RData` (md5 `46b4f143d0289489ea0e9277035bb26a`) |

All 216 fits converged.

## Refactor verification log (2026-11-24)

The design matrix and candidate strings can now be built by the formula-based code (`R/TCJ_CHP_designmat.R` with `HOR_CHP_designmat()`, and the legacy adapter `legacy_TCJ_CHP_design()`), with `TCJ_CHP_precomp()` and `TCJ_CHP_modelprep()` in `R/TCJ_CHP_comp_ALT.R`. `reference/verify_TCJ_CHP_refactor.R` repeats the checks (set `run_fits <- TRUE` to refit; about 1 minute).

| Check | Result against the files in this folder |
|---|---|
| `TCJ_CHP_comp_ALT()` (split functions) | List `identical()` to `TCJ_CHP_comp_ls.rds`, scaled and unscaled |
| Formula-based design, legacy mode | `allint_DMs` (216 strings, same order), all 13 columns and the whole list `identical()` to `TCJ_CHP_comp_ls.rds`. The unscaled variant matches in `all.equal()` but is not `identical()`; the reason was not investigated |
| Refit of all 216 candidates | nll, AIC and parameters exactly equal to `TCJ_CHP_TMB_mod_fits_ls.rds` |
| Pipeline steps `TCJ_CHP_cand_mods_run_init.R` and `TCJ_CHP_conf_comp.R` with the new design | `TCJ_CHP_TMB_all_mods.rds` equal in every element, confidence set (5 models) and its coefficients and standard errors exactly equal. Outputs were written to `output/TCJ_CHP_check_legacyA/` (git-ignored, not frozen) |

The `WYT_drought` fix is applied in `../TCJ_CHP_WYTfixed_2026-11-24/` (`legacy_TCJ_CHP_design(precomp, wyt = "water.yr")`): AIC values shift by 6 on average, the best model and the 5-model confidence set are unchanged.

## Using these to verify another repository (e.g. CVPASapp)

```r
ref <- "reference/TCJ_CHP_legacy_design_2026-11-23"
old_fits <- readRDS(file.path(ref, "TCJ_CHP_TMB_mod_fits_ls.rds"))
old_comp <- readRDS(file.path(ref, "TCJ_CHP_comp_ls.rds"))
identical(new_comp$allint_DMs, old_comp$allint_DMs)
all.equal(new_comp$XX_in_w_int_WYT, old_comp$XX_in_w_int_WYT, tolerance = 1e-8)
max(abs(vapply(new_fits, function(x) x$AIC, 0) - vapply(old_fits, function(x) x$AIC, 0)))
```

Tolerances: see the HOR_CHP legacy README.
