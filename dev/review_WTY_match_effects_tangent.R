library(ggplot2)

ref_leg <- "reference/HOR_CHP_Rxtemp_Tmsd_2026-11-24"
ref_fix <- "reference/HOR_CHP_Rxtemp_Tmsd_WYTfixed_2026-11-24"

fits_leg <- readRDS(file.path(ref_leg, "HOR_CHP_TMB_mod_fits_ls.rds"))
fits_fix <- readRDS(file.path(ref_fix, "HOR_CHP_TMB_mod_fits_ls.rds"))
comp_leg <- readRDS(file.path(ref_leg, "HOR_CHP_comp_ls.rds"))
comp_fix <- readRDS(file.path(ref_fix, "HOR_CHP_comp_ls.rds"))

dms <- comp_leg$allint_DMs
stopifnot(identical(dms, comp_fix$allint_DMs))

# release groups: label, true water-year type, legacy (lookup) water-year type
xx <- comp_leg$XX_in
wyt_cols <- comp_leg$XX_in_w_int_WYT
grp_rows <- data.frame(
  idx = comp_leg$TMB_data_baseline$relgrps_ind,
  label = paste(xx$year, xx$rel.grp),
  water.yr = as.character(xx$water.yr),
  wyt_legacy = ifelse(wyt_cols$WYT_wet == 1, "wet", ifelse(wyt_cols$WYT_drought == 1, "drought", "dry")))
grp <- unique(grp_rows)
grp <- grp[order(grp$idx), ]
grp$n_rows <- as.vector(table(grp_rows$idx)[as.character(grp$idx)])
grp$mislabelled <- grp$water.yr != grp$wyt_legacy
stopifnot(nrow(grp) == 19, !anyDuplicated(grp$idx))

pull_model <- function(fit, id){
  par <- fit[[id]]$OPT$par
  s <- par[names(par) == "S_pars"]
  list(re = fit[[id]]$RELGRP_err_v,
       logSD = unname(par["logSD_RELGRP"]),
       b_wet = unname(s[2]), b_drought = unname(s[3]),
       nll = fit[[id]]$nll)
}

# candidate models to compare (same candidate string in both fits)
id_best_leg <- which.min(vapply(fits_leg, function(x) x$AIC, 0))
id_best_fix <- which.min(vapply(fits_fix, function(x) x$AIC, 0))
ids <- c(base = 1, best_legacy = id_best_leg, best_fixed = id_best_fix, global = length(dms))

model_summary <- do.call(rbind, lapply(names(ids), function(nm){
  a <- pull_model(fits_leg, ids[[nm]]); b <- pull_model(fits_fix, ids[[nm]])
  data.frame(model = nm, candidate = ids[[nm]], dm = dms[ids[[nm]]],
             nll_legacy = a$nll, nll_fixed = b$nll,
             sd_RE_legacy = sd(a$re), sd_RE_fixed = sd(b$re),
             logSD_legacy = a$logSD, logSD_fixed = b$logSD,
             b_wet_legacy = a$b_wet, b_wet_fixed = b$b_wet,
             b_drought_legacy = a$b_drought, b_drought_fixed = b$b_drought)
}))
model_summary

# group-level detail for the best legacy candidate
detail_for <- function(id){
  a <- pull_model(fits_leg, id); b <- pull_model(fits_fix, id)
  out <- grp
  out$re_legacy <- a$re
  out$re_fixed <- b$re
  out$wyt_effect_legacy <- a$b_wet*(out$wyt_legacy == "wet") + a$b_drought*(out$wyt_legacy == "drought")
  out$wyt_effect_fixed <- b$b_wet*(out$water.yr == "wet") + b$b_drought*(out$water.yr == "drought")
  out$total_legacy <- out$wyt_effect_legacy + out$re_legacy
  out$total_fixed <- out$wyt_effect_fixed + out$re_fixed
  out$re_change <- out$re_fixed - out$re_legacy
  out
}
detail <- detail_for(id_best_leg)
detail[, c("label", "n_rows", "water.yr", "wyt_legacy", "re_legacy", "re_fixed", "re_change", "wyt_effect_legacy", "wyt_effect_fixed", "total_legacy", "total_fixed")]

# mean random effect by water-year type, legacy and fixed fits
aggregate(cbind(re_legacy, re_fixed) ~ wyt_legacy, data = detail, FUN = mean)
aggregate(cbind(re_legacy, re_fixed) ~ water.yr, data = detail, FUN = mean)
aggregate(cbind(re_legacy, re_fixed, re_change) ~ mislabelled, data = detail, FUN = function(x) mean(abs(x)))

# are the fixed-fit random effects aligned with the legacy label?
legacy_drought <- as.numeric(detail$wyt_legacy == "drought")
legacy_wet <- as.numeric(detail$wyt_legacy == "wet")
c(cor_RE_fixed_with_legacy_drought = cor(detail$re_fixed, legacy_drought),
  cor_RE_fixed_with_legacy_wet = cor(detail$re_fixed, legacy_wet),
  cor_RE_fixed_with_true_drought = cor(detail$re_fixed, as.numeric(detail$water.yr == "drought")),
  cor_RE_legacy_with_RE_fixed = cor(detail$re_legacy, detail$re_fixed),
  cor_total_legacy_with_total_fixed = cor(detail$total_legacy, detail$total_fixed))

# same across every candidate: does the fixed fit need larger random effects?
sd_re <- data.frame(
  legacy = vapply(fits_leg, function(x) sd(x$RELGRP_err_v), 0),
  fixed = vapply(fits_fix, function(x) sd(x$RELGRP_err_v), 0))
summary(sd_re$fixed / sd_re$legacy)
summary(vapply(seq_along(fits_fix), function(i) fits_fix[[i]]$nll - fits_leg[[i]]$nll, 0))

ggplot(detail, aes(x = re_legacy, y = reorder(label, idx))) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(aes(xend = re_fixed, yend = reorder(label, idx)), colour = "grey60") +
  geom_point(aes(colour = "legacy WYT")) +
  geom_point(aes(x = re_fixed, colour = "fixed WYT")) +
  facet_grid(mislabelled ~ ., scales = "free_y", space = "free_y", labeller = label_both) +
  labs(x = "release-group random effect (best legacy candidate)", y = NULL, colour = NULL)

ggplot(detail, aes(x = total_legacy, y = total_fixed, label = label)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(shape = mislabelled)) +
  geom_text(vjust = -0.7, size = 3) +
  labs(x = "WYT effect + random effect, legacy fit", y = "WYT effect + random effect, fixed fit")
