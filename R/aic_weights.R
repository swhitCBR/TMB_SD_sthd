# R/aic_weights.R
#
# AIC weights for the candidate-model tables. Deduplicated from
# HOR_cand_mods_run.R / TCJ_cand_mods_run.R / HOR_TCJ_cand_mods_run.R,
# which each defined an identical copy.

#' AIC weights (Akaike weights)
#'
#' Computes Akaike weights for the models in `df` from the AIC values in
#' `aic_column`: relative likelihoods `exp(-0.5 * delta AIC)` normalized to sum
#' to 1.
#'
#' @param df Data frame with one row per model.
#' @param aic_column Name of the column holding the AIC values.
#'
#' @return Numeric vector of weights, aligned with the rows of `df`.
#' @export
get_aic_wts <- function(df, aic_column) {
  delta_aic <- df[[aic_column]] - min(df[[aic_column]], na.rm = TRUE)
  rel_likelihood <- exp(-0.5 * delta_aic)
  rel_likelihood / sum(rel_likelihood, na.rm = TRUE)
}
