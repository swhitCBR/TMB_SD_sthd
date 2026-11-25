legacy_col_names <- c(
  "barrier.facTRUE:Tmsd.hor.7dadm"="B_x_Temp",
  "barrier.facTRUE:log.VNS.hor.5"="B_x_VNS",
  "barrier.facTRUE:SWP.hor.5"="B_x_SWP",
  "barrier.facTRUE:CVP.hor.5"="B_x_CVP",
  "barrier.facTRUE:Qomt.hor.1net"="B_x_OMT",
  "route.facB:Tmsd.hor.7dadm"="R_x_Temp",
  "route.facB:log.VNS.hor.5"="R_x_VNS",
  "route.facB:SWP.hor.5"="R_x_SWP",
  "route.facB:CVP.hor.5"="R_x_CVP",
  "route.facB:Qomt.hor.1net"="R_x_OMT")

legacy_HOR_CHP_design <- function(precomp,
                                  formula = HOR_CHP_global_formula,
                                  fixed = HOR_CHP_fixed_terms,
                                  r_x_temp = c("Tclc.hor.3", "Tmsd.hor.7dadm"),
                                  wyt = c("legacy", "water.yr")){

  r_x_temp <- match.arg(r_x_temp)
  wyt <- match.arg(wyt)

  dat <- precomp$XX_in
  dat$barrier.fac <- factor(dat$barrier, levels = c(FALSE, TRUE))

  if(wyt == "legacy"){
    wyt_mat <- precomp$XX_in_w_int_WYT[, c("WYT_wet", "WYT_drought")]
    stopifnot(all(rowSums(wyt_mat) <= 1))
    dat$WYT_ <- factor(c("dry", "wet", "drought")[1 + wyt_mat$WYT_wet + 2*wyt_mat$WYT_drought],
                       levels = c("dry", "wet", "drought"))
  } else {
    dat$WYT_ <- factor(dat$water.yr, levels = c("dry", "wet", "drought"))
  }

  col_overrides <- NULL
  if(r_x_temp == "Tclc.hor.3"){
    col_overrides <- list("route.facB:Tmsd.hor.7dadm" = as.numeric(dat$route.fac == "B") * as.vector(dat$Tclc.hor.3))
  }

  HOR_CHP_designmat(dat, formula = formula, fixed = fixed,
                    col_overrides = col_overrides, col_names = legacy_col_names)
}

legacy_TCJ_col_names <- c(
  "TCJroute.facF:SWP.tcj.4"="R_x_SWP",
  "TCJroute.facF:CVP.tcj.4"="R_x_CVP",
  "TCJroute.facF:Qomt.tcj.1net"="R_x_OMT")

legacy_TCJ_CHP_design <- function(precomp,
                                  formula = TCJ_CHP_global_formula,
                                  fixed = TCJ_CHP_fixed_terms,
                                  wyt = c("legacy", "water.yr")){

  wyt <- match.arg(wyt)

  dat <- precomp$XX_in
  dat$TCJroute.fac <- factor(c("A", "F")[1 + dat$TCJroute.facF], levels = c("A", "F"))

  if(wyt == "legacy"){
    wyt_mat <- precomp$XX_in_w_int_WYT[, c("WYT_wet", "WYT_drought")]
    stopifnot(all(rowSums(wyt_mat) <= 1))
    dat$WYT_ <- factor(c("dry", "wet", "drought")[1 + wyt_mat$WYT_wet + 2*wyt_mat$WYT_drought],
                       levels = c("dry", "wet", "drought"))
  } else {
    dat$WYT_ <- factor(dat$water.yr, levels = c("dry", "wet", "drought"))
  }

  out <- HOR_CHP_designmat(dat, formula = formula, fixed = fixed, col_names = legacy_TCJ_col_names)
  out$XX_s <- legacy_copy_col_attrs(out$XX_s, precomp$XX_in_w_int_WYT)
  out
}

legacy_copy_col_attrs <- function(new, legacy){
  for(nm in intersect(names(new), names(legacy))){
    if(is.matrix(legacy[[nm]])) attributes(new[[nm]]) <- attributes(legacy[[nm]])
  }
  new
}

compare_designs <- function(new, legacy){
  new <- as.matrix(new)
  legacy <- as.matrix(legacy)
  stopifnot(identical(dim(new), dim(legacy)))
  data.frame(position = seq_len(ncol(new)),
             new_name = colnames(new),
             legacy_name = colnames(legacy),
             same_name = colnames(new) == colnames(legacy),
             max_abs_diff = vapply(seq_len(ncol(new)), function(j) max(abs(new[, j] - legacy[, j])), 0))
}
