# R/glmmTMB_par.R
#
# Parallel versions of the glmmTMB steps in the model_fitting/*_cand_mods_run.R scripts:
#   - dredge_par_glmmTMB()  replaces MuMIn::dredge() on a glmmTMB global model
#   - get_models_par()    replaces MuMIn::get.models(<dredge table>, subset = NA)
#   - fit_glmmTMB_par()   replaces an lapply() over a set of glmmTMB formulas
# glmm_cluster() starts the worker cluster the three functions use.
#
# Each function starts its own cluster (and stops it again) unless one is passed in `cl`.
# The workers refit the models from their calls, so the data frame a model was fitted to
# is copied to every worker under the same name it has in the calling session. That name
# has to be a variable in `envir` (default: the caller); dredge_par_glmmTMB() also needs it to be a
# global variable, as it is in the scripts.
#
# Fits are deterministic, so the models equal those of the serial code (AIC and fixed
# effects agree to numerical precision); only the order in which they are computed differs.
#
# Usage:  source("R/glmmTMB_par.R")

#' Name of the data variable used by a glmmTMB model or call
#'
#' @param x A fitted glmmTMB model, or a call to `glmmTMB()`.
#' @return The `data` argument as a character string.
#' @noRd
data_name_of <- function(x) {
  cl <- if (is.call(x)) x else x$call
  d <- cl$data
  if (!is.symbol(d)) {
    stop("the model's `data` argument must be a variable name (e.g. data = my.df), not an expression",
         call. = FALSE)
  }
  as.character(d)
}

#' Start a worker cluster for fitting glmmTMB models
#'
#' Starts a PSOCK cluster, gives the workers the main session's library paths and
#' `na.action` option, loads glmmTMB on them, and copies the named data frames into their
#' global environments, so calls such as `glmmTMB(y ~ x, data = my.df)` can be evaluated
#' there.
#'
#' @param data_names Character vector of names of the objects to copy to the workers.
#' @param n_workers Number of worker processes. Capped at the number of logical cores
#'   (with a warning).
#' @param envir Environment the objects named in `data_names` are looked up in.
#'
#' @return A `cluster` object. Stop it with `parallel::stopCluster()` when finished.
#' @export
#' @examples
#' \dontrun{
#' cl <- glmm_cluster("hor.df", n_workers = 12)
#' parallel::stopCluster(cl)
#' }
glmm_cluster <- function(data_names = character(), n_workers = 12L, envir = parent.frame()) {
  n_workers <- as.integer(n_workers)
  stopifnot(length(n_workers) == 1L, !is.na(n_workers), n_workers >= 1L)
  n_cores <- parallel::detectCores(logical = TRUE)
  if (!is.na(n_cores) && n_workers > n_cores) {
    warning("n_workers (", n_workers, ") exceeds the number of logical cores; using ", n_cores,
            call. = FALSE)
    n_workers <- n_cores
  }

  # --no-init-file: workers started in the project folder would each run .Rprofile and
  # activate renv (about 20 s for 12 workers); the library paths are set explicitly below
  cl <- parallel::makeCluster(n_workers, type = "PSOCK", rscript_args = "--no-init-file")
  ready <- FALSE
  on.exit(if (!ready) parallel::stopCluster(cl), add = TRUE)

  init_worker <- function(lib_paths, na_action) {
    .libPaths(lib_paths)
    options(na.action = na_action)
    suppressPackageStartupMessages(library(glmmTMB))
    invisible(NULL)
  }
  # a closure is sent with its environment; the global environment is not serialised, whereas
  # this function's frame would drag the caller's objects (e.g. the dredge table) to every worker
  environment(init_worker) <- globalenv()
  parallel::clusterCall(cl, init_worker, .libPaths(), getOption("na.action"))

  for (nm in data_names) {
    value <- get(nm, envir = envir, inherits = TRUE)
    parallel::clusterCall(cl, assign, nm, value, envir = globalenv())
  }

  ready <- TRUE
  cl
}

#' Parallel `MuMIn::dredge()` for a glmmTMB global model
#'
#' Runs `MuMIn::dredge()` with its `cluster` argument on workers that already have the
#' model's data. The result is the same `model.selection` table `dredge()` returns without
#' a cluster, with the same model numbering (row names).
#'
#' `dredge()` evaluates the data variable from inside the MuMIn namespace when it runs on a
#' cluster, so the model's data frame must be a global variable.
#'
#' @param global_model A fitted glmmTMB model whose call has `data = <variable name>`.
#' @param ... Further arguments for `MuMIn::dredge()` (`fixed`, `beta`, `rank`, ...).
#' @inheritParams get_models_par
#'
#' @return A `model.selection` object, as from `MuMIn::dredge()`.
#' @export
#' @examples
#' \dontrun{
#' m <- glmmTMB(y ~ a * b + (1 | g), family = "binomial", data = my.df)
#' dd <- dredge_par_glmmTMB(m, fixed = "a", beta = "none")
#' }
dredge_par_glmmTMB <- function(global_model, ..., cl = NULL, n_workers = 12L, envir = parent.frame()) {
  model_expr <- substitute(global_model)
  dots <- as.list(substitute(list(...)))[-1L]
  if (is.null(cl)) {
    cl <- glmm_cluster(data_name_of(global_model), n_workers = n_workers, envir = envir)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  # dredge() looks up its model and cluster arguments by name from inside MuMIn, so they are
  # passed as values (a cluster, and the model, are objects it can use directly)
  dredge_call <- as.call(c(list(quote(MuMIn::dredge), global_model), dots, list(cluster = cl)))
  res <- eval(dredge_call, envir)

  # record the call as it would read in the script, not with the whole model inlined
  attr(res, "call") <- as.call(c(list(quote(MuMIn::dredge), global.model = model_expr), dots))

  # MuMIn numbers the models of a cluster run from 0 and a serial run from 1. Use the serial
  # numbering so the row names (which the scripts use to pick out models) do not change.
  idx <- suppressWarnings(as.integer(rownames(res)))
  if (length(idx) && !anyNA(idx) && min(idx) == 0L) {
    new_names <- as.character(idx + 1L)
    attr(res, "row.names") <- new_names
    model_calls <- attr(res, "model.calls")
    names(model_calls) <- new_names
    attr(res, "model.calls") <- model_calls
  }
  res
}

#' Parallel `MuMIn::get.models()` for a dredge table
#'
#' Refits the models of a `model.selection` table from `MuMIn::dredge()` on the workers; the
#' parallel equivalent of `MuMIn::get.models(object, subset = NA)`. The returned list has
#' the table's row names as names, in the table's order, so `models[rownames(object)]`
#' lines up with the table.
#'
#' @param object A `model.selection` table from `MuMIn::dredge()` on a glmmTMB model whose
#'   call has `data = <variable name>`.
#' @param rows Optional rows of `object` to refit (row positions or row names). Default: all.
#' @param cl A cluster from [glmm_cluster()] that already holds the model's data, or `NULL`
#'   to start (and afterwards stop) one with `n_workers` workers.
#' @param n_workers Number of workers when `cl` is `NULL`.
#' @param envir Environment the data variable is looked up in when `cl` is `NULL`.
#'
#' @return A named list of fitted models with the `rank` and `beta` attributes `get.models()`
#'   sets.
#' @export
#' @examples
#' \dontrun{
#' dd <- MuMIn::dredge(global_model, beta = "none")
#' models <- get_models_par(dd)
#' }
get_models_par <- function(object, rows = NULL, cl = NULL, n_workers = 12L, envir = parent.frame()) {
  calls <- attr(object, "model.calls")
  if (!inherits(object, "model.selection") || is.null(calls)) {
    stop("`object` must be a model.selection table from dredge() (it has no model calls)", call. = FALSE)
  }
  if (!is.null(rows)) {
    if (is.character(rows)) rows <- match(rows, rownames(object))
    if (anyNA(rows)) stop("`rows` not found in rownames(object)", call. = FALSE)
    calls <- calls[rows]
    nms <- rownames(object)[rows]
  } else {
    nms <- rownames(object)
  }

  if (is.null(cl)) {
    cl <- glmm_cluster(data_name_of(attr(object, "global.call")), n_workers = n_workers, envir = envir)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  # clusterApplyLB hands the next model to whichever worker is free
  models <- parallel::clusterApplyLB(cl, calls, eval, envir = globalenv())
  names(models) <- nms
  for (a in c("rank", "beta")) attr(models, a) <- attr(object, a)
  models
}

#' Fit a set of glmmTMB models in parallel
#'
#' Fits `glmmTMB(<response> ~ <rhs> [+ <random>], data = <data_name>, ...)` once per element
#' of `rhs`, spread over the workers. The calls are evaluated in the global environment of
#' each worker, so the fitted models do not carry this function's environment (and its data).
#'
#' @param rhs Named or unnamed character vector (or list) of formula right-hand sides, e.g.
#'   `c(year = "Year", VNS = "log.VNS.hor.5")`. The result takes its names.
#' @param data_name Name of the data frame (a variable in `envir`) to fit to.
#' @param response Left-hand side of the formulas, e.g. `"detd"`.
#' @param random Optional random-effect term appended to every right-hand side, e.g.
#'   `"(1|YrRel)"`.
#' @param ... Further arguments passed to `glmmTMB()` (`family`, `na.action`, ...). They are
#'   not evaluated here: like the formula and data, they are written into each model's call
#'   and evaluated on the worker, so they must not refer to variables that exist only in the
#'   calling session.
#' @inheritParams get_models_par
#'
#' @return A list of fitted glmmTMB models, one per element of `rhs`, with the names of `rhs`.
#' @export
#' @examples
#' \dontrun{
#' fits <- fit_glmmTMB_par(c(null = "1", year = "Year"), "my.df", response = "detd",
#'                         random = "(1|YrRel)", family = binomial(link = "logit"))
#' }
fit_glmmTMB_par <- function(rhs, data_name, response, random = NULL, ...,
                            cl = NULL, n_workers = 12L, envir = parent.frame()) {
  rhs <- unlist(rhs)
  formulas <- paste(response, "~", rhs, if (!is.null(random)) paste("+", random))
  names(formulas) <- names(rhs)
  extra_args <- match.call(expand.dots = FALSE)$...

  calls <- lapply(formulas, function(formula_text) {
    f <- stats::as.formula(formula_text, env = globalenv())
    as.call(c(list(as.name("glmmTMB"), f, data = as.name(data_name)), extra_args))
  })

  if (is.null(cl)) {
    cl <- glmm_cluster(data_name, n_workers = n_workers, envir = envir)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  # clusterApplyLB hands the next model to whichever worker is free
  models <- parallel::clusterApplyLB(cl, calls, eval, envir = globalenv())
  names(models) <- names(formulas)
  models
}
