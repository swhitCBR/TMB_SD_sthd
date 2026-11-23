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
