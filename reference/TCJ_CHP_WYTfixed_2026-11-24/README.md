# TCJ_CHP results with `WYT_` from `water.yr` (2026-11-24)

Frozen output of the TCJ_CHP steps after the `WYT_wet` and `WYT_drought` dummies were changed to follow each row's own water-year type (`water.yr`). The release-group lookup used before assigned the wrong type to 6 of 18 release groups (see `../TCJ_CHP_legacy_design_2026-11-23/README.md`). Do not edit or regenerate these files; add a new dated folder if the design changes again.

TCJ_CHP fits use the same TMB template as HOR_CHP (`HOR_CHP_global`).

## What changed relative to `../TCJ_CHP_legacy_design_2026-11-23/`

Only the `WYT_drought` column. The candidate strings (`allint_DMs`, 216 strings of 13 characters), the column order and all other columns are identical. `WYT_` terms are in every candidate, so every model changes.

| Item | Legacy | This folder |
|---|---|---|
| Best model (lowest AIC) | `1111110010100`, AIC 3255.26 | same model, AIC 3262.38 |
| Change in AIC per candidate | | absolute change median 6.0, maximum 8.1 |
| Rank correlation (Spearman) of AIC across the 216 models | | 0.99 |
| Confidence set | 5 models | the same 5 models (`1111110010100`, `1111110011100`, `1111111010100`, `1111110110100`, `1111110010101`) |

All 216 fits converged. The models have the same number of parameters in both versions, so the AIC differences are differences in log likelihood. In contrast to HOR_CHP, the ranking and the confidence set are essentially unchanged.

Files and layout are as in the legacy folder: `TCJ_CHP_comp_ls.rds`, `TCJ_CHP_comp_ls_unscl.rds` (12 `NA` values in the design matrix, 6 rows with a missing `Qomt.tcj.1net`, so not fit), `TCJ_CHP_TMB_mod_fits_ls.rds`, `TCJ_CHP_TMB_all_mods.rds` (its `XX_in_w_int_WYT` element is `NULL`), `TCJ_CHP_S_conf_ls.rds`, `MANIFEST.csv`.

## How these were produced

| Item | Value |
|---|---|
| Date and machine | 2026-11-24 21:28-21:30 local, NERVA |
| Command | `out_dir <- "output/TCJ_CHP_WYTfixed"` and `source("model_fitting/runscript.R")` with every step skipped except `TCJ_CHP_cand_mods_run_init.R` and `TCJ_CHP_conf_comp.R`; 12 workers |
| Setting | `CHP_wyt <- "water.yr"` (default in `model_fitting/config.R`); `"legacy"` reproduces the release-group lookup |
| Git | HEAD `1047562` plus uncommitted edits (the `R/` folder, `model_fitting/` scripts and `config.R`) |
| R / packages | R 4.5.2, TMB 1.9.18, glmmTMB 1.1.14, MuMIn 1.48.19, TMBhelp 1.0.0 |
| Input | `data/TCJ_CHP_mod_dat_ls.RData`, md5 `46b4f143d0289489ea0e9277035bb26a` |

## Verifying another repository against this folder

```r
ref <- "reference/TCJ_CHP_WYTfixed_2026-11-24"
old_fits <- readRDS(file.path(ref, "TCJ_CHP_TMB_mod_fits_ls.rds"))
old_comp <- readRDS(file.path(ref, "TCJ_CHP_comp_ls.rds"))
identical(new_comp$allint_DMs, old_comp$allint_DMs)
all.equal(new_comp$XX_in_w_int_WYT, old_comp$XX_in_w_int_WYT, tolerance = 1e-8)
max(abs(vapply(new_fits, function(x) x$AIC, 0) - vapply(old_fits, function(x) x$AIC, 0)))
```

Use the tolerances in the HOR_CHP legacy README.
