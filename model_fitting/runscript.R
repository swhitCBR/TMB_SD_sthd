# model_fitting/runscript.R
#
# Runs the whole model-fitting pipeline in order: data preparation, candidate model
# fits, then confidence sets / standard errors. It first sources model_fitting/config.R.
#
# This script is NOT meant to be run casually:
#   - the data preparation steps read inputs from tmp_data/, else ../CVPAS-steelhead/ (see input_dirs in
#     config.R) and overwrite files in data/
#   - the candidate-fit steps overwrite results in `out_dir` (default "output/"); the
#     glmmTMB scripts (HOR, TCJ, HOR_TCJ) save unconditionally
#   - the HOR_CHP and TCJ_CHP fits run in parallel (the workers use the installed TMBhelp,
#     which must be the release pinned by tmbhelp_version in config.R)
#   - the complete run takes on the order of 30 minutes or more
# To test, point the results at a scratch folder before sourcing, e.g.
#   out_dir <- "C:/some/scratch/folder"; dir.create(out_dir)
# (the confidence-set scripts read their inputs from out_dir, so copy
# TCJ_CHP_TMB_all_mods.rds and HOR_CHP_S_conf_ls.rds there first).
#
# Usage, from the repository root:
#   source("model_fitting/runscript.R")
#
# Settings (define before sourcing to override):
#   out_dir    folder the steps write to / read from      (default "output")
#   n_workers  worker processes for the parallel fits      (default 12, set in config.R)
#   skip       character vector of step file names to skip (default none)
#
# The steps share one R session and one global environment, as the scripts were
# written to be run that way. Each step stops the run if it fails, because later steps
# depend on earlier ones.

stopifnot("run from the repository root (TMB_SD_sthd)" = file.exists("TMB_SD_sthd.Rproj"))

source("model_fitting/config.R")

# the pipeline is pinned to one TMBhelp release (the parallel workers use the installed copy)
tmbhelp_installed <- if (requireNamespace("TMBhelp", quietly = TRUE)) as.character(utils::packageVersion("TMBhelp")) else "none"
if (tmbhelp_installed != tmbhelp_version) {
  stop("TMBhelp ", tmbhelp_version, " is required (installed: ", tmbhelp_installed,
       "). Install with renv::restore(), or renv::install(\"swhitCBR/TMBhelp@v", tmbhelp_version, "\").")
}

# ---- steps, in run order ---------------------------------------------------- #
steps <- c(
  # 1. data preparation (the *_data_prep.R scripts also define the *_comp() functions;
  #    run_data_prep is set to TRUE for these steps so the preparation itself runs)
  "HOR_TCJ_separate_route_mods_data_prep.R",
  "HOR_TCJ_surv_mods_data_prep.R",
  "HOR_CHP_surv_mods_data_prep.R",
  "TCJ_CHP_surv_mods_data_prep.R",
  # compile the TMB templates (before any parallel cluster starts)
  "compile_templates.R",
  # 2. candidate model fits (all candidates, not only confidence-set members)
  "HOR_cand_mods_run.R",
  "TCJ_cand_mods_run.R",
  "HOR_TCJ_cand_mods_run.R",
  "HOR_CHP_cand_mods_run.R",
  "TCJ_CHP_cand_mods_run_init.R",
  # 3. confidence sets and standard errors (read the results of step 2 from out_dir)
  "TCJ_CHP_conf_comp.R",
  "HOR_CHP_conf_comp.R",
  "HOR_CHP_d2_AIC_unscl_mod_ls.R"
)

# Scripts in this folder that are not steps: they define functions or are sourced by the
# steps themselves (TCJ_CHP_cust_fxns.R), or are this script and its config.
not_steps <- c("runscript.R", "config.R", "TCJ_CHP_cust_fxns.R")

# Every script in the folder must be accounted for, so a new script cannot be left out silently.
in_folder <- list.files("model_fitting", pattern = "\\.[Rr]$")
unlisted <- setdiff(in_folder, c(steps, not_steps))
if (length(unlisted)) {
  stop("scripts in model_fitting/ that are not in `steps` or `not_steps`: ",
       paste(unlisted, collapse = ", "), ". Add them to runscript.R.")
}
missing_steps <- setdiff(steps, in_folder)
if (length(missing_steps)) {
  stop("steps listed in runscript.R but not found in model_fitting/: ",
       paste(missing_steps, collapse = ", "))
}

# ---- helper: the confidence-set steps fit the TMB templates in this R session ------- #
# TMB::MakeADFun() fails with "Multiple TMB models loaded" when more than one model DLL is
# loaded. glmmTMB's own DLL (loaded by the glmmTMB steps) and the templates loaded by an earlier
# step (HOR_CHP_global, HOR_CHP_global_pred) therefore have to be unloaded before each of these
# steps. Nothing after the glmmTMB steps uses glmmTMB.
tmb_session_steps <- c("TCJ_CHP_conf_comp.R", "HOR_CHP_conf_comp.R", "HOR_CHP_d2_AIC_unscl_mod_ls.R")
unload_tmb_model_dlls <- function() {
  loaded <- getLoadedDLLs()
  for (nm in intersect(c("glmmTMB", "HOR_CHP_global", "HOR_CHP_global_pred", "HOR_CHP_predmod"), names(loaded))) {
    dyn.unload(loaded[[nm]][["path"]])
  }
}

# ---- run -------------------------------------------------------------------- #
run_t0 <- Sys.time()
for (step in setdiff(steps, skip)) {
  message(format(Sys.time(), "%H:%M:%S"), "  running ", step)
  step_t0 <- Sys.time()
  # the data_prep scripts run their preparation only when run_data_prep is TRUE; the
  # fitting scripts set it to FALSE (functions only) when they source a data_prep script
  if (grepl("_data_prep\\.R$", step)) run_data_prep <- TRUE
  if (step %in% tmb_session_steps) unload_tmb_model_dlls()
  source(file.path("model_fitting", step), local = FALSE)
  message("    finished ", step, " in ",
          round(as.numeric(difftime(Sys.time(), step_t0, units = "mins")), 1), " min")
}
message("all steps finished in ",
        round(as.numeric(difftime(Sys.time(), run_t0, units = "mins")), 1), " min")
