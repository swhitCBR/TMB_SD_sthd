# model_fitting/config.R
#
# Settings for the model-fitting pipeline, as a named list. Sourced at the top of
# model_fitting/runscript.R, which copies the list's elements into the global
# environment (the steps share one session and its global environment).
#
# To override a setting for a run, source this file, edit pipeline_config, then
# source runscript.R:
#   source("model_fitting/config.R")
#   pipeline_config$out_dir <- "C:/some/scratch/folder"; dir.create(pipeline_config$out_dir)
#   source("model_fitting/runscript.R")
# (the confidence-set scripts read their inputs from out_dir, so copy
#  TCJ_CHP_TMB_all_mods.rds and HOR_CHP_S_conf_ls.rds there first).
#
# input_path(), which locates the upstream input files using input_dirs below,
# is now defined in R/input_path.R.

pipeline_config <- list(
  # folder the steps write results to / read results from
  out_dir = "output",

  # worker processes for every parallel step (MuMIn::dredge clusters in the glmmTMB scripts and
  # fit_nested_tmb_mods_par() in the HOR_CHP / TCJ_CHP fits)
  n_workers = 14,

  # variable multiplied with the route indicator in the HOR_CHP R_x_Temp interaction ("Tclc.hor.3" reproduces the legacy design)
  HOR_CHP_r_x_temp = "Tmsd.hor.7dadm",

  # water-year type behind the WYT_ dummies in the HOR_CHP and TCJ_CHP designs: "water.yr" (each row's own type) or "legacy" (release-group lookup with an indexing bug)
  CHP_wyt = "water.yr",

  # Is Year (a categorical study-year effect) a term of the HOR_TCJ survival global model?
  # TRUE (default): 582 candidates; Year models fit best and fill the d2 set, so predictions only
  #   exist for the study years. FALSE: 388 candidates (the Year-free subset of the 582, same AICs;
  #   WYT models compete only with each other). With FALSE the HOR_TCJ outputs get the suffix
  #   "_noYear" (HOR_TCJ_all_mods_noYear.rds, ...) so the standard files are not overwritten.
  hor_tcj_include_year = TRUE,

  # step file names to skip (character vector)
  skip = character(),

  # TMBhelp release the pipeline was built and checked against (checked by runscript.R).
  # To change it: install the new release, set this value, run the pipeline and compare
  # the results with the previous output/ before committing (see README.md).
  tmbhelp_version = "1.0.0",

  # Folders searched, in order, for upstream input files. The first folder containing the
  # requested file is used, so a local copy in tmp_data/ takes precedence over the
  # CVPAS-steelhead repository next to this one.
  input_dirs = c("tmp_data", "../CVPAS-steelhead/tmp_data", "../CVPAS-steelhead"),

  # Confidence-set definitions for the additional *_conf_sets.rds outputs. Each entry is
  # one rule for keeping a subset of the candidate models of an AIC table:
  #   "d<N>"   all models with delta AIC < N ("d2" = the classic set, i.e. what the
  #            "*_d2_*" objects and output files elsewhere in this repository mean)
  #   "wt<P>"  the smallest set of models, best AIC first, whose cumulative AIC weights
  #            reach P percent
  # Currently used by the glmmTMB steps (HOR, TCJ, HOR_TCJ); the TMB steps still use
  # the single dAIC_threshold.R. The functions behind the definitions are in R/conf_sets.R.
  conf_set_defs = c("d2", "wt25"),

  # Apply the "nesting rule" (Richards 2008; Richards, Whittingham & Stephens 2011) to
  # the candidate sets BEFORE the confidence-set definitions above: models that are more
  # complex versions of a lower-AIC model are excluded first, and dAIC and AIC weights are
  # then computed over the survivors (which can shrink the confidence sets a lot; on the
  # saved candidate tables the rule keeps 6 of 1024 HOR_CHP and 16 of 216 TCJ_CHP models).
  # Currently only supported for candidate tables with a "dm" column of 0/1 design
  # strings (the HOR_CHP / TCJ_CHP TMB tables); the glmmTMB dredge tables would need
  # MuMIn::nested(). The rule itself is in R/nest_rule_mods.R.
  nesting_rule = FALSE
)
