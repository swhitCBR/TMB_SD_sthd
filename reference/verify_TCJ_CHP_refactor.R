stopifnot("run from the repository root (TMB_SD_sthd)" = file.exists("TMB_SD_sthd.Rproj"))

if (!exists("run_fits")) run_fits <- FALSE
if (!exists("n_workers")) source("model_fitting/config.R")

ref_dir <- "reference/TCJ_CHP_legacy_design_2026-11-23"
ref_comp <- readRDS(file.path(ref_dir, "TCJ_CHP_comp_ls.rds"))
ref_comp_unscl <- readRDS(file.path(ref_dir, "TCJ_CHP_comp_ls_unscl.rds"))
ref_fits <- readRDS(file.path(ref_dir, "TCJ_CHP_TMB_mod_fits_ls.rds"))

for(f in c("R/HOR_CHP_designmat.R", "R/TCJ_CHP_designmat.R", "R/legacy_adapters.R", "R/TCJ_CHP_comp_ALT.R")) source(f)

precomp <- TCJ_CHP_precomp()
precomp_unscl <- TCJ_CHP_precomp(z_scale_vars = FALSE)

stopifnot(identical(TCJ_CHP_comp_ALT(), ref_comp))
stopifnot(identical(TCJ_CHP_comp_ALT(z_scale_vars = FALSE), ref_comp_unscl))

design <- legacy_TCJ_CHP_design(precomp)
comp <- TCJ_CHP_modelprep(precomp, design = design)
stopifnot(identical(comp, ref_comp))
stopifnot(identical(design$allint_DMs, ref_comp$allint_DMs))
comp_unscl <- TCJ_CHP_modelprep(precomp_unscl, design = legacy_TCJ_CHP_design(precomp_unscl))
stopifnot(isTRUE(all.equal(comp_unscl, ref_comp_unscl, tolerance = 0)))

# WYT_ from water.yr: only the WYT_drought column may change
design_C <- legacy_TCJ_CHP_design(precomp, wyt = "water.yr")
cmp_C <- compare_designs(design_C$XX_s, design$XX_s)
stopifnot(identical(cmp_C$new_name[cmp_C$max_abs_diff > 0], "WYT_drought"))
stopifnot(identical(design_C$allint_DMs, design$allint_DMs))

if(run_fits){
  library(TMB)
  library(TMBhelp)
  fits <- fit_nested_tmb_mods_par(
    S_design_matrix_in = comp$allint_DMs,
    dir_in = file.path(getwd(), "TMB"),
    TMB_data_baseline_in = comp$TMB_data_baseline,
    calc_SEs = FALSE, cpp_name_in = "HOR_CHP_global",
    n_workers = n_workers)
  stopifnot(max(abs(vapply(fits, `[[`, 0, "AIC") - vapply(ref_fits, `[[`, 0, "AIC"))) < 1e-6)
}
message("TCJ_CHP refactor checks passed")
