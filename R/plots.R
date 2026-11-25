# R/plots.R
#
# Diagnostic and results plots, wrapped as functions so the pipeline steps stay
# readable and the plots can be reused with other model sets. Each function
# prints its plots (matching the old top-level behavior of the scripts) and
# returns them invisibly.

#' HOR-TCJ covariate diagnostics
#'
#' Histograms of the HOR-TCJ survival-model covariates by water-year type (and
#' by release group for the Z-scaled data) and pairwise plots of the flow and
#' temperature covariates. Moved from HOR_TCJ_cand_mods_run.R.
#'
#' @param kf.df Unscaled HOR-TCJ survival model data frame.
#' @param kf.dfZ Z-scaled version of the same data.
#'
#' @return Invisibly, a list of the plots.
#' @export
plot_HOR_TCJ_diagnostics <- function(kf.df, kf.dfZ) {
  plots <- list(
    ggplot2::ggplot(data = kf.df, ggplot2::aes(x = VNS.hor.5)) +
      ggplot2::facet_wrap(~WYT) + ggplot2::geom_histogram(),
    ggplot2::ggplot(data = kf.df, ggplot2::aes(x = log.VNS.hor.5)) +
      ggplot2::facet_wrap(~WYT) + ggplot2::geom_histogram(),
    ggplot2::ggplot(data = kf.dfZ, ggplot2::aes(x = log.VNS.hor.5, fill = rel.grp)) +
      ggplot2::facet_wrap(~WYT) + ggplot2::geom_histogram(),
    ggplot2::ggplot(data = kf.df, ggplot2::aes(x = Tmsd.hor.7dadm)) +
      ggplot2::facet_wrap(~WYT) + ggplot2::geom_histogram(),
    ggplot2::ggplot(data = kf.df, ggplot2::aes(x = Tmsd.hor.7dadm, fill = rel.grp)) +
      ggplot2::facet_wrap(~WYT) + ggplot2::geom_histogram(),
    ggplot2::ggplot(data = kf.dfZ, ggplot2::aes(x = Tmsd.hor.7dadm, fill = rel.grp)) +
      ggplot2::facet_wrap(~WYT) + ggplot2::geom_histogram())
  invisible(lapply(plots, print))
  plot(kf.df$VNS.hor, kf.df$Tmsd.hor.7dadm)
  plot(kf.df$log.VNS.hor.5, kf.df$Tmsd.hor.7dadm)
  plot(kf.df$Tmsd.hor.7dadm, kf.df$log.VNS.hor.5)
  invisible(plots)
}

#' TCJ export-covariate correlation plot
#'
#' Old exploratory plot of the two export indicators. Moved from
#' TCJ_cand_mods_run.R (the accompanying `cor()` call stays in the script).
#'
#' @param tcj.df TCJ routing model data frame.
#' @export
plot_TCJ_exit_correlation <- function(tcj.df) {
  plot(tcj.df$SWP.exit.tcj, tcj.df$CVP.exit.tcj)
}

#' Confidence-set coefficient plot
#'
#' Estimates with 95% intervals of the confidence-set models, point size scaled
#' by AIC weight and color by model. Moved from TCJ_CHP_conf_comp.R.
#'
#' @param coef_df Data frame of coefficient estimates with columns `S_parID`,
#'   `par_nm`, `estimate`, `SE`, `AICwt` and `dm` (one row per model x
#'   coefficient), as built by `get_SEcoef_tab()` in TCJ_CHP_conf_comp.R.
#'
#' @return Invisibly, the plot.
#' @export
plot_confset_coefs <- function(coef_df) {
  p <- ggplot2::ggplot(data = coef_df,
         ggplot2::aes(x = factor(S_parID),
             y = estimate, ymin = estimate - SE * 1.96, ymax = estimate + SE * 1.96,
             size = AICwt * 2.5,
             fill = factor(dm))) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", colour = "gray40") +
    ggplot2::geom_linerange(width = 0, position = ggplot2::position_dodge(width = 1), size = 0.5) +
    ggplot2::geom_point(shape = 21, position = ggplot2::position_dodge(width = 1)) +
    ggplot2::scale_x_discrete(limits = rev, labels = rev(unique(coef_df$par_nm))) +
    ggplot2::coord_flip()
  print(p)
  invisible(p)
}
