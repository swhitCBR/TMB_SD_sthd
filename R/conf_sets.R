# R/conf_sets.R
#
# Confidence-set definitions for the model-fitting pipeline. A confidence set is
# a rule for keeping a subset of the candidate models of an AIC table; the rules
# are character strings taken from pipeline_config$conf_set_defs
# (model_fitting/config.R):
#   "d<N>"   all models with delta AIC < N. "d2" is the classic delta AIC < 2 set
#            (the "*_d2_*" objects and output files elsewhere in this repository
#            refer to this definition).
#   "wt<P>"  the smallest set of models, best AIC first, whose cumulative AIC
#            weights reach P percent (e.g. "wt25" = cumulative weight >= 0.25).
#
# Optionally the nesting rule (R/nest_rule_mods.R) is applied to the candidates
# first (pipeline_config$nesting_rule), so that models with non-informative extra
# parameters are excluded before the definitions are evaluated.
#
# Used by the glmmTMB steps (HOR_cand_mods_run.R, TCJ_cand_mods_run.R,
# HOR_TCJ_cand_mods_run.R) and the TMB candidate-fit steps (HOR_CHP_cand_mods_run.R,
# TCJ_CHP_cand_mods_run_init.R) to write the *_conf_sets.rds outputs.

if (!exists("get_aic_wts")) source("R/aic_weights.R")

#' Parse one confidence-set definition string
#'
#' @param def Definition string, e.g. `"d2"`, `"d6"`, `"wt25"`, `"wt90"`.
#'
#' @return List with elements `type` (`"d"` or `"wt"`) and `value` (numeric:
#'   the delta-AIC threshold, or the cumulative-weight target as a fraction).
#' @export
parse_conf_set_def <- function(def) {
  m <- regmatches(def, regexec("^d([0-9.]+)$", def))[[1]]
  if (length(m)) return(list(type = "d", value = as.numeric(m[2])))
  m <- regmatches(def, regexec("^wt([0-9.]+)$", def))[[1]]
  if (length(m)) return(list(type = "wt", value = as.numeric(m[2]) / 100))
  stop("unrecognized confidence-set definition: '", def,
       "' (expected 'd<N>' or 'wt<P>', e.g. 'd2', 'wt25')", call. = FALSE)
}

#' Membership mask of one confidence set
#'
#' @param aic_tab AIC table (one row per model), sorted by increasing AIC.
#' @param def Definition string (see [parse_conf_set_def()]).
#'
#' @return Named logical vector aligned with the rows of `aic_tab`.
#' @export
conf_set_mask <- function(aic_tab, def) {
  if (is.unsorted(aic_tab$AIC)) {
    stop("conf_set_mask() expects the AIC table sorted by increasing AIC", call. = FALSE)
  }
  parsed <- parse_conf_set_def(def)
  if (parsed$type == "d") {
    return(aic_tab$dAIC < parsed$value)
  }
  # wt: keep models, best AIC first, until the cumulative AIC weight reaches the
  # target (the crossing model is included, so the set's weight >= target)
  wt_cumul <- cumsum(get_aic_wts(aic_tab, "AIC"))
  hit <- which(wt_cumul >= parsed$value)
  seq_len(nrow(aic_tab)) <= (if (length(hit)) min(hit) else nrow(aic_tab))
}

#' Compile all confidence sets of one model set
#'
#' Subsets an AIC table and its model list by every definition in `defs` and
#' builds the meta-confidence set (the models meeting at least one definition,
#' with one logical column per definition showing which ones it meets).
#'
#' When `nesting_rule` is TRUE, the nesting rule ([apply_nest_rule()], R/nest_rule_mods.R)
#' is applied to the candidate table first, so that models which are more complex versions
#' of a lower-AIC model are excluded before the definitions are evaluated; dAIC and AIC
#' weights are then computed over the surviving models. This needs a `design_col` column
#' of 0/1 design strings (as in the HOR_CHP / TCJ_CHP TMB candidate tables).
#'
#' @param aic_tab AIC table (one row per model), with rownames matching the
#'   names of `mods_ls` (when supplied). It is sorted by increasing AIC internally.
#' @param mods_ls Named list of fitted models, aligned with the rownames of
#'   `aic_tab`; or NULL for a table-only compilation (the TMB steps keep their
#'   fitted models separately).
#' @param defs Character vector of definition strings (see
#'   [parse_conf_set_def()]). Defaults to `pipeline_config$conf_set_defs`.
#' @param nesting_rule Apply the nesting rule before the definitions? Defaults
#'   to `pipeline_config$nesting_rule` when it exists, else FALSE.
#' @param design_col Column holding the 0/1 design strings used by the nesting rule.
#'
#' @return List with elements:
#'   \item{defs}{the definitions used, in order}
#'   \item{sets}{named list with one entry per definition, each a list with
#'     `aictab` (the subset AIC table) and, when `mods_ls` was supplied, `mods`
#'     (the matching models)}
#'   \item{meta}{data frame of all models in at least one set, with a logical
#'     column per definition and the number of definitions met (`n_defs`)}
#'   \item{n_removed_by_nesting}{number of models removed by the nesting rule
#'     (0 when `nesting_rule` is FALSE)}
#' @export
compile_conf_sets <- function(aic_tab, mods_ls = NULL,
                              defs = if (exists("pipeline_config") && !is.null(pipeline_config$conf_set_defs)) pipeline_config$conf_set_defs else "d2",
                              nesting_rule = if (exists("pipeline_config") && !is.null(pipeline_config$nesting_rule)) pipeline_config$nesting_rule else FALSE,
                              design_col = "dm") {
  n_removed <- 0L
  if (nesting_rule) {
    if (!exists("apply_nest_rule")) source("R/nest_rule_mods.R")
    if (!design_col %in% names(aic_tab)) {
      stop("compile_conf_sets(): nesting_rule = TRUE needs a '", design_col,
           "' column of 0/1 design strings (as in the TMB candidate tables); ",
           "for dredge() tables use MuMIn::nested() instead", call. = FALSE)
    }
    n_before <- nrow(aic_tab)
    aic_tab <- apply_nest_rule(aic_tab, design_col = design_col, aic_col = "AIC")
    n_removed <- attr(aic_tab, "n_removed_by_nesting")
    message("nesting rule: removed ", n_removed, " of ", n_before, " candidate models")
    # apply_nest_rule() drops the rownames; re-key the rows by the design strings
    rownames(aic_tab) <- aic_tab[[design_col]]
  }
  aic_tab <- aic_tab[order(aic_tab$AIC), , drop = FALSE]
  if (is.null(rownames(aic_tab)) && !is.null(mods_ls)) {
    stop("compile_conf_sets(): aic_tab needs rownames matching names(mods_ls)", call. = FALSE)
  }

  masks <- lapply(defs, conf_set_mask, aic_tab = aic_tab)
  names(masks) <- defs

  in_any <- Reduce("|", masks)
  if (!is.null(mods_ls)) {
    mod_nms <- rownames(aic_tab)[in_any]
    missing <- setdiff(mod_nms, names(mods_ls))
    if (length(missing)) {
      stop("compile_conf_sets(): no model in mods_ls for row(s): ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
  }

  sets <- lapply(seq_along(masks), function(i) {
    out <- list(aictab = aic_tab[masks[[i]], , drop = FALSE])
    if (!is.null(mods_ls)) out$mods <- mods_ls[rownames(aic_tab)[masks[[i]]]]
    out
  })
  names(sets) <- defs

  meta <- cbind(aic_tab[in_any, , drop = FALSE],
                do.call(cbind, lapply(masks, function(m) m[in_any])))
  meta$n_defs <- rowSums(as.matrix(meta[defs]))

  list(defs = defs, sets = sets, meta = meta, n_removed_by_nesting = n_removed)
}
