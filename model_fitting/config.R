# model_fitting/config.R
#
# Settings for the model-fitting pipeline; sourced first by model_fitting/runscript.R.
# Each setting is only defined if it does not already exist, so values defined before
# sourcing runscript.R (e.g. out_dir <- "C:/scratch") take precedence.

# folder the steps write results to / read results from
if (!exists("out_dir")) out_dir <- "output"

# worker processes for the parallel TMB fits
if (!exists("n_workers")) n_workers <- 8

# step file names to skip (character vector)
if (!exists("skip")) skip <- character()

# TMBhelp release the pipeline was built and checked against (checked by runscript.R).
# To change it: install the new release, set this value, run the pipeline and compare
# the results with the previous output/ before committing (see README.md).
if (!exists("tmbhelp_version")) tmbhelp_version <- "1.0.0"

# Folders searched, in order, for upstream input files. The first folder containing the
# requested file is used, so a local copy in tmp_data/ takes precedence over the
# CVPAS-steelhead repository next to this one. Define input_dirs before sourcing to override.
if (!exists("input_dirs")) {
  input_dirs <- c("tmp_data", "../CVPAS-steelhead/tmp_data", "../CVPAS-steelhead")
}

# input_path("Survival model/Multiyear v2a AllVar v3.RData") returns the path of the first
# existing match in input_dirs (and says which one it chose); it stops with the list of
# folders searched if the file is in none of them.
input_path <- function(rel_file, dirs = input_dirs) {
  cand <- file.path(dirs, rel_file)
  hit <- cand[file.exists(cand)]
  if (!length(hit)) {
    stop("input file '", rel_file, "' not found in: ", paste(dirs, collapse = ", "), call. = FALSE)
  }
  message("input: ", normalizePath(hit[1], winslash = "/"))
  hit[1]
}
