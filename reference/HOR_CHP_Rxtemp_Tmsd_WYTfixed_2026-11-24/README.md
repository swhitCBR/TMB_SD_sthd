# HOR_CHP results with `R_x_Temp` from `Tmsd.hor.7dadm` and `WYT_` from `water.yr` (2026-11-24)

Frozen output of the HOR_CHP steps after both design changes: the route x temperature interaction uses `Tmsd.hor.7dadm`, and the `WYT_wet` and `WYT_drought` dummies follow each row's own water-year type (`water.yr`) instead of the release-group lookup, which assigned the wrong type to 6 of 19 release groups (see `../HOR_CHP_legacy_design_2026-11-23/README.md`). Do not edit or regenerate these files; add a new dated folder if the design changes again.

This is the current design. Earlier states:

| Folder | `R_x_Temp` | `WYT_` |
|---|---|---|
| `../HOR_CHP_legacy_design_2026-11-23/` | route x `Tclc.hor.3` | release-group lookup (6 groups mislabelled) |
| `../HOR_CHP_Rxtemp_Tmsd_2026-11-24/` | route x `Tmsd.hor.7dadm` | release-group lookup (6 groups mislabelled) |
| this folder | route x `Tmsd.hor.7dadm` | `water.yr` |

## What changed relative to `../HOR_CHP_Rxtemp_Tmsd_2026-11-24/`

Only the `WYT_drought` column (the `WYT_wet` column was already correct). The candidate strings (`allint_DMs`, 1,024), the column order and all other columns are identical. `WYT_` terms are in every candidate, so every model changes.

| Item | Rxtemp_Tmsd (legacy `WYT_`) | This folder |
|---|---|---|
| Best model (lowest AIC) | candidate 833, `111111111110000001011`, AIC 7345.51 | `111111111110010001011`, AIC 7364.59 |
| Change in AIC per candidate | | all 1,024 models are worse: +11.3 to +22.3, median +17.4 |
| Models within 2 AIC of the best | 14 | 16 (6 in common with the previous 14) |
| Rank correlation (Spearman) of AIC across the 1,024 models | | 0.62 |
| Confidence set (`dAIC_threshold`) | 14 models | 16 models |

All 1,024 fits converged. The models have the same number of parameters in both versions, so the AIC differences are differences in log likelihood: the release-group lookup, although wrong, fits these data better than the correct water-year types. This has not been investigated beyond reporting it.

Files and layout are as in the legacy folder (21 columns of `XX_in_w_int_WYT`, `allint_DMs` of 1,024 strings of 21 characters): `HOR_CHP_comp_ls.rds`, `HOR_CHP_comp_ls_unscl.rds`, `HOR_CHP_TMB_mod_fits_ls.rds`, `HOR_CHP_TMB_all_mods.rds` (its `XX_in_w_int_WYT` element is `NULL`), `HOR_CHP_S_conf_ls.rds`, `HOR_CHP_TMB_all_mods_unscl.rds`, `MANIFEST.csv`.

## How these were produced

| Item | Value |
|---|---|
| Date and machine | 2026-11-24 21:14-21:28 local, NERVA |
| Command | `out_dir <- "output/HOR_CHP_WYTfixed"` and `source("model_fitting/runscript.R")` with every step skipped except `HOR_CHP_cand_mods_run.R`, `HOR_CHP_conf_comp.R` and `HOR_CHP_d2_AIC_unscl_mod_ls.R`; 12 workers |
| Settings | `HOR_CHP_r_x_temp <- "Tmsd.hor.7dadm"` and `CHP_wyt <- "water.yr"` (defaults in `model_fitting/config.R`); `CHP_wyt <- "legacy"` reproduces the release-group lookup |
| Git | HEAD `1047562` plus uncommitted edits (the `R/` folder, `model_fitting/` scripts and `config.R`) |
| R / packages | R 4.5.2, TMB 1.9.18, glmmTMB 1.1.14, MuMIn 1.48.19, TMBhelp 1.0.0 |
| Input | `data/HOR_CHP_mod_dat_ls.RData`, md5 `a57048a22f235f663f815f66d2e4d390` |

The unscaled refit step prints `NA/NaN function evaluation` warnings during optimisation, as it did in the other folders.

## Verifying another repository against this folder

```r
ref <- "reference/HOR_CHP_Rxtemp_Tmsd_WYTfixed_2026-11-24"
old_fits <- readRDS(file.path(ref, "HOR_CHP_TMB_mod_fits_ls.rds"))
old_comp <- readRDS(file.path(ref, "HOR_CHP_comp_ls.rds"))
identical(new_comp$allint_DMs, old_comp$allint_DMs)
all.equal(new_comp$XX_in_w_int_WYT, old_comp$XX_in_w_int_WYT, tolerance = 1e-8)
max(abs(vapply(new_fits, function(x) x$AIC, 0) - vapply(old_fits, function(x) x$AIC, 0)))
```

Use the tolerances in the legacy README.
