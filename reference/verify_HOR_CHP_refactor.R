stopifnot("run from the repository root (TMB_SD_sthd)" = file.exists("TMB_SD_sthd.Rproj"))

if (!exists("run_fits")) run_fits <- FALSE
if (!exists("n_workers")) source("model_fitting/config.R")

ref_dir <- "reference/HOR_CHP_legacy_design_2026-11-23"
ref_comp <- readRDS(file.path(ref_dir, "HOR_CHP_comp_ls.rds"))
ref_comp_unscl <- readRDS(file.path(ref_dir, "HOR_CHP_comp_ls_unscl.rds"))
ref_fits <- readRDS(file.path(ref_dir, "HOR_CHP_TMB_mod_fits_ls.rds"))

source("R/scale_data_cols.R")
for(f in c("R/HOR_CHP_designmat.R", "R/legacy_adapters.R", "R/HOR_CHP_comp_ALT.R")) source(f)

# Stage A: refactored code in legacy mode reproduces the frozen design
pc <- suppressMessages(HOR_CHP_precomp())
pu <- suppressMessages(HOR_CHP_precomp(z_scale_vars = FALSE))
stopifnot(identical(HOR_CHP_comp_ALT(), ref_comp))
design_A <- legacy_HOR_CHP_design(pc)
comp_A <- HOR_CHP_modelprep(pc, design = design_A)
stopifnot(identical(comp_A, ref_comp))
stopifnot(identical(HOR_CHP_modelprep(pu, design = legacy_HOR_CHP_design(pu))$TMB_data_baseline, ref_comp_unscl$TMB_data_baseline))

# Stage B: R_x_Temp uses Tmsd.hor.7dadm; only that column may change
design_B <- legacy_HOR_CHP_design(pc, r_x_temp = "Tmsd.hor.7dadm")
cmp_B <- compare_designs(design_B$XX_s, ref_comp$XX_in_w_int_WYT)
stopifnot(identical(cmp_B$new_name[cmp_B$max_abs_diff > 0], "R_x_Temp"))
stopifnot(identical(design_B$allint_DMs, ref_comp$allint_DMs))
comp_B <- HOR_CHP_modelprep(pc, design = design_B)

# Stage C: WYT_ from water.yr; only the WYT_drought column may change
design_C <- legacy_HOR_CHP_design(pc, r_x_temp = "Tmsd.hor.7dadm", wyt = "water.yr")
cmp_C <- compare_designs(design_C$XX_s, design_B$XX_s)
stopifnot(identical(cmp_C$new_name[cmp_C$max_abs_diff > 0], "WYT_drought"))
stopifnot(identical(design_C$allint_DMs, design_B$allint_DMs))

if(run_fits){
  library(TMB)
  library(TMBhelp)
  fit_with <- function(comp) fit_nested_tmb_mods_par(
    S_design_matrix_in = comp$allint_DMs,
    dir_in = file.path(getwd(), "TMB"),
    TMB_data_baseline_in = comp$TMB_data_baseline,
    calc_SEs = FALSE, cpp_name_in = "HOR_CHP_global",
    n_workers = n_workers)
  fits_A <- fit_with(comp_A)
  stopifnot(max(abs(vapply(fits_A, `[[`, 0, "AIC") - vapply(ref_fits, `[[`, 0, "AIC"))) < 1e-6)
  fits_B <- fit_with(comp_B)
  dAIC_B <- vapply(fits_B, `[[`, 0, "AIC") - vapply(ref_fits, `[[`, 0, "AIC")
  has_R_x_Temp <- substr(comp_B$allint_DMs, 17, 17) == "1"
  stopifnot(max(abs(dAIC_B[!has_R_x_Temp])) < 1e-6)
}
message("HOR_CHP refactor checks passed")
