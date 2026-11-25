# HOR_CHP results with `R_x_Temp` using `Tmsd.hor.7dadm` (2026-11-24)

Frozen output of the HOR_CHP steps of the pipeline after the design matrix was moved to the formula-based code (`R/HOR_CHP_designmat.R`, with `R/legacy_adapters.R`) and the route x temperature interaction was changed to use `Tmsd.hor.7dadm`, the same temperature variable as the main effect and `B_x_Temp`. Do not edit or regenerate these files; add a new dated folder if the design changes again.

Compare with `../HOR_CHP_legacy_design_2026-11-23/`, which holds the results before the change.

## What differs from the legacy folder

| Item | Legacy (2026-11-23) | This folder |
|---|---|---|
| `R_x_Temp` | route indicator x `Tclc.hor.3` | route indicator x `Tmsd.hor.7dadm` |
| `WYT_drought` | release-group lookup (still the legacy values, see the legacy README) | unchanged here; fixed in `../HOR_CHP_Rxtemp_Tmsd_WYTfixed_2026-11-24/` |
| Other 19 design columns, column order, `allint_DMs` (1,024 strings) | | identical to legacy |
| Fits without `R_x_Temp` (512 models) | | exactly the same as legacy |
| Fits with `R_x_Temp` (512 models) | | changed: median absolute change in AIC 0.32, maximum 4.34 |
| Best model | candidate 833 | same |
| Confidence set (`dAIC_threshold`) | 13 models | 14 models: the 13 legacy models plus `111111111111000111011` |

Files are the same set as in the legacy folder (see its README for contents and the design-matrix column layout): `HOR_CHP_comp_ls.rds`, `HOR_CHP_comp_ls_unscl.rds`, `HOR_CHP_TMB_mod_fits_ls.rds`, `HOR_CHP_TMB_all_mods.rds`, `HOR_CHP_S_conf_ls.rds`, `HOR_CHP_TMB_all_mods_unscl.rds`, and `MANIFEST.csv` (size and md5). As in the legacy folder, the `XX_in_w_int_WYT` element of `HOR_CHP_TMB_all_mods.rds` is `NULL`.

## How these were produced

| Item | Value |
|---|---|
| Date and machine | 2026-11-24 18:25-18:37 local, NERVA |
| Command | `out_dir <- "output/HOR_CHP_Rxtemp_Tmsd"` and `source("model_fitting/runscript.R")` with every step skipped except `HOR_CHP_cand_mods_run.R`, `HOR_CHP_conf_comp.R` and `HOR_CHP_d2_AIC_unscl_mod_ls.R`; 12 workers |
| Setting | `HOR_CHP_r_x_temp <- "Tmsd.hor.7dadm"` (default in `model_fitting/config.R`); `"Tclc.hor.3"` reproduces the legacy design |
| Git | HEAD `1047562` plus uncommitted edits (the `R/` folder, `model_fitting/` scripts and `config.R`) |
| R / packages | R 4.5.2, TMB 1.9.18, glmmTMB 1.1.14, MuMIn 1.48.19, TMBhelp 1.0.0 |
| Inputs | same data intermediates as the legacy folder (`data/HOR_CHP_mod_dat_ls.RData`, md5 `a57048a22f235f663f815f66d2e4d390`) |

All 1,024 fits converged. The fits are exactly equal to a separate in-session run of the same design (`reference/verify_HOR_CHP_refactor.R`, Stage B). The unscaled refits of the confidence set agree with the scaled AIC values to 7e-7 at most (legacy: 3e-7). The unscaled refit step prints `NA/NaN function evaluation` warnings during optimisation; I did not compare their number with the legacy run.

## Verifying another repository against this folder

```r
ref <- "reference/HOR_CHP_Rxtemp_Tmsd_2026-11-24"
old_fits <- readRDS(file.path(ref, "HOR_CHP_TMB_mod_fits_ls.rds"))
old_comp <- readRDS(file.path(ref, "HOR_CHP_comp_ls.rds"))
identical(new_comp$allint_DMs, old_comp$allint_DMs)
all.equal(new_comp$XX_in_w_int_WYT, old_comp$XX_in_w_int_WYT, tolerance = 1e-8)
max(abs(vapply(new_fits, function(x) x$AIC, 0) - vapply(old_fits, function(x) x$AIC, 0)))
```

Use the tolerances in the legacy README.
