# R/input_path.R
#
# Locate upstream input files for the model-fitting pipeline. Moved out of
# model_fitting/config.R (which now holds settings only) so config.R is purely data.

#' Path of an upstream input file
#'
#' Returns the path of the first existing match of `rel_file` in `dirs`, and
#' messages which folder it chose. Stops with the list of folders searched if
#' the file is in none of them.
#'
#' @param rel_file File name (or relative path) of the input file, e.g.
#'   `"Multiyear v2a AllVar v3.RData"`.
#' @param dirs Folders searched, in order; the first folder containing the file
#'   is used. Defaults to `pipeline_config$input_dirs` (model_fitting/config.R),
#'   so a local copy in `tmp_data/` takes precedence over the CVPAS-steelhead
#'   repository next to this one. Falls back to those same folders when
#'   `pipeline_config` is not in the global environment.
#'
#' @return Character scalar: the path of the first existing match.
#' @export
input_path <- function(rel_file,
                       dirs = if (exists("pipeline_config", envir = globalenv(), inherits = FALSE)) {
                         globalenv()$pipeline_config$input_dirs
                       } else {
                         c("tmp_data", "../CVPAS-steelhead/tmp_data", "../CVPAS-steelhead")
                       }) {
  cand <- file.path(dirs, rel_file)
  hit <- cand[file.exists(cand)]
  if (!length(hit)) {
    stop("input file '", rel_file, "' not found in: ", paste(dirs, collapse = ", "), call. = FALSE)
  }
  message("input: ", normalizePath(hit[1], winslash = "/"))
  hit[1]
}
