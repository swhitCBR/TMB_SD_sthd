# R/nest_rule_mods.R
#
# The "nesting rule" (Richards 2008; Richards, Whittingham & Stephens 2011): a model
# that is a more complex version of a model with a lower AIC has non-informative extra
# parameters and is excluded from model selection / averaging. Copied from
# CVPAS_STH_app/R/nest_rule_mods.R; controlled in this pipeline by the `nesting_rule`
# setting in model_fitting/config.R and applied by compile_conf_sets() (R/conf_sets.R)
# before the confidence-set definitions.

#' Apply the nesting rule to a table of candidate models
#'
#' Removes every model that is a more complex version of a model with a lower
#' AIC (the "nesting rule"). A model is dropped when some other model in the
#' table (i) contains only terms that are also in it and (ii) has a strictly
#' lower AIC. Optionally applies a delta AIC cutoff.
#'
#' @details
#' Nesting is determined from the 0/1 inclusion string of each model (one
#' character per term in the full design matrix, "1" = term included), so it is
#' valid when candidate models differ only in which design-matrix columns they
#' use, with the same random effects and error structure. A model is nested in
#' another when its set of included terms is a subset of the other's.
#'
#' The best-ranked model is always retained. Ties in AIC are not treated as
#' "better", so tied models are both kept. Akaike weights in `nest_wt` are
#' recomputed over the retained models only, so they no longer sum to one over
#' the original candidate set. Memory use grows with the square of the number of
#' models (a few thousand models is fine).
#'
#' For tables produced by `MuMIn::dredge()`, `MuMIn::nested()` implements the
#' same rule from the model formulas, e.g. `subset(dredge_table, !nested(.))`;
#' this function is for models identified by inclusion strings.
#'
#' @param aic_tab Data frame with one row per candidate model.
#' @param design_col Name of the column holding the 0/1 inclusion strings; all
#'   strings must have the same length.
#' @param aic_col Name of the column holding AIC values. Rows with missing AIC
#'   are removed with a warning.
#' @param delta_max Optional cutoff on delta AIC (relative to the best model in
#'   `aic_tab`) applied in addition to the nesting rule. The default `Inf`
#'   applies the nesting rule only.
#'
#' @return The retained rows of `aic_tab`, ordered by AIC, with added columns
#'   `delta_aic` (relative to the best model) and `nest_wt` (Akaike weights
#'   renormalised over the retained models). Attributes `n_input` and
#'   `n_removed_by_nesting` record the number of models supplied and the number
#'   removed by the nesting rule (before any delta AIC cutoff).
#'
#' @references
#' Richards, S. A. (2008). Dealing with overdispersed count data in applied
#' ecology. Journal of Applied Ecology, 45, 218-227.
#'
#' Richards, S. A., Whittingham, M. J., & Stephens, P. A. (2011). Model
#' selection and model averaging in behavioural ecology: the utility of the
#' IT-AIC framework. Behavioral Ecology and Sociobiology, 65, 77-89.
#'
#' @seealso [MuMIn::nested()] for the equivalent filter on `dredge()` tables.
#'
#' @examples
#' tab <- data.frame(design = c("110", "111", "100", "101"),
#'                   AIC = c(10, 10.4, 12, 11))
#' # "111" is dropped: it is a more complex version of "110", which has lower AIC
#' apply_nest_rule(tab)
#' apply_nest_rule(tab, delta_max = 1.5)
apply_nest_rule <- function(aic_tab, design_col = "design", aic_col = "AIC",
                            delta_max = Inf) {
  stopifnot(is.data.frame(aic_tab),
            all(c(design_col, aic_col) %in% names(aic_tab)),
            is.numeric(delta_max), length(delta_max) == 1, delta_max >= 0)

  design <- as.character(aic_tab[[design_col]])
  aic <- aic_tab[[aic_col]]
  if (anyNA(design) || !all(grepl("^[01]+$", design)) ||
      length(unique(nchar(design))) != 1) {
    stop("`", design_col, "` must hold equal-length strings of 0s and 1s")
  }

  has_aic <- !is.na(aic)
  if (!all(has_aic)) {
    warning(sum(!has_aic), " model(s) with missing AIC were removed")
    aic_tab <- aic_tab[has_aic, , drop = FALSE]
    design <- design[has_aic]
    aic <- aic[has_aic]
  }

  ord <- order(aic)
  aic_tab <- aic_tab[ord, , drop = FALSE]
  design <- design[ord]
  aic <- aic[ord]

  X <- (do.call(rbind, strsplit(design, "")) == "1") * 1L
  # nested_in[j, i]: every term of model j is also in model i (j != i)
  nested_in <- (X %*% t(1L - X)) == 0
  diag(nested_in) <- FALSE
  # lower_aic[j, i]: model j has a strictly lower AIC than model i
  lower_aic <- outer(aic, aic, "<")
  dominated <- colSums(nested_in & lower_aic) > 0

  delta <- aic - min(aic)
  keep <- !dominated & delta <= delta_max

  out <- aic_tab[keep, , drop = FALSE]
  out$delta_aic <- delta[keep]
  w <- exp(-0.5 * out$delta_aic)
  out$nest_wt <- w / sum(w)
  rownames(out) <- NULL
  attr(out, "n_input") <- length(aic)
  attr(out, "n_removed_by_nesting") <- sum(dominated)
  out
}
