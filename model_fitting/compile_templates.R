# model_fitting/compile_templates.R
# Compiles the TMB templates the pipeline uses (TMB/HOR_CHP_global.cpp and
# TMB/HOR_CHP_global_pred.cpp) if no compiled library exists. fit_nested_tmb_mods_par()
# needs pre-compiled templates, and it must be compiled before any cluster starts.
# Flags are those used by TMBhelp::TMB_model_setup(). Delete the .dll/.so/.o files in
# TMB/ to force a recompile.

library(TMB)

for (tmpl in c("HOR_CHP_global", "HOR_CHP_global_pred")) {
  dll <- file.path("TMB", paste0(tmpl, .Platform$dynlib.ext))
  if (file.exists(dll)) {
    message("template already compiled: ", dll)
  } else {
    message("compiling ", tmpl)
    old_wd <- setwd("TMB")
    tryCatch(TMB::compile(paste0(tmpl, ".cpp"), "-O1 -g", DLLFLAGS = ""),
             finally = setwd(old_wd))
    stopifnot("compilation did not produce a library" = file.exists(dll))
  }
}
