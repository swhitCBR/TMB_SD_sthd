HOR_CHP_global_formula <- ~ WYT_ + route.fac + barrier.fac + flength + Tmsd.hor.7dadm + log.VNS.hor.5 + SWP.hor.5 + CVP.hor.5 + Qomt.hor.1net +
  barrier.fac*(Tmsd.hor.7dadm + log.VNS.hor.5 + SWP.hor.5 + CVP.hor.5 + Qomt.hor.1net) +
  route.fac*(Tmsd.hor.7dadm + log.VNS.hor.5 + SWP.hor.5 + CVP.hor.5 + Qomt.hor.1net)

HOR_CHP_fixed_terms <- c("WYT_", "route.fac", "barrier.fac", "flength", "Tmsd.hor.7dadm", "log.VNS.hor.5", "SWP.hor.5", "CVP.hor.5", "Qomt.hor.1net")

sort_term <- function(x) vapply(strsplit(x, ":", fixed = TRUE), function(p) paste(sort(p), collapse = ":"), "")

HOR_CHP_designmat <- function(data, formula = HOR_CHP_global_formula, fixed = HOR_CHP_fixed_terms,
                              col_overrides = NULL, col_names = NULL){

  mm <- model.matrix(formula, model.frame(formula, data, na.action = na.pass))
  assign_ids <- attr(mm, "assign")
  term_labels <- attr(terms(formula), "term.labels")

  dat_tmp <- data[complete.cases(data[, all.vars(formula)]), ]
  dat_tmp$DUMMY_Y <- rep(c(0, 1), length.out = nrow(dat_tmp))
  glob_tmp <- lm(update(formula, DUMMY_Y ~ .), data = dat_tmp, na.action = na.fail)

  cand_calls <- MuMIn::dredge(glob_tmp, fixed = fixed, evaluate = FALSE)
  cand_terms <- lapply(cand_calls, function(cl) attr(terms(eval(cl$formula)), "term.labels"))
  term_in <- t(sapply(cand_terms, function(m) as.integer(sort_term(term_labels) %in% sort_term(m))))

  # first varying term changes fastest
  varying <- setdiff(seq_along(term_labels), match(sort_term(fixed), sort_term(term_labels)))
  cand_order <- do.call(order, rev(lapply(varying, function(j) term_in[, j])))

  col_in <- term_in[, pmax(assign_ids, 1), drop = FALSE]
  col_in[, assign_ids == 0] <- 1L
  allint_DMs <- unname(apply(col_in[cand_order, , drop = FALSE], 1, paste, collapse = ""))

  XX_s <- as.data.frame(mm, check.names = FALSE)
  for(nm in names(col_overrides)){
    stopifnot(nm %in% names(XX_s))
    XX_s[[nm]] <- col_overrides[[nm]]
  }
  if(!is.null(col_names)){
    stopifnot(all(names(col_names) %in% names(XX_s)))
    names(XX_s)[match(names(col_names), names(XX_s))] <- col_names
  }

  list("XX_s"=XX_s, "allint_DMs"=allint_DMs)
}
