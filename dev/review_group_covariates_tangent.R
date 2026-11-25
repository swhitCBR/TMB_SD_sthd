library(ggplot2)

ref_leg <- "reference/HOR_CHP_Rxtemp_Tmsd_2026-11-24"
ref_fix <- "reference/HOR_CHP_Rxtemp_Tmsd_WYTfixed_2026-11-24"

fits_leg <- readRDS(file.path(ref_leg, "HOR_CHP_TMB_mod_fits_ls.rds"))
fits_fix <- readRDS(file.path(ref_fix, "HOR_CHP_TMB_mod_fits_ls.rds"))
comp <- readRDS(file.path(ref_leg, "HOR_CHP_comp_ls.rds"))

# group-level survival effect from the fixed fit of the best legacy candidate: WYT effect + random effect
id_best_leg <- which.min(vapply(fits_leg, function(x) x$AIC, 0))
g <- unique(data.frame(idx = comp$TMB_data_baseline$relgrps_ind,
                       label = paste(comp$XX_in$year, comp$XX_in$rel.grp),
                       water.yr = as.character(comp$XX_in$water.yr)))
g <- g[order(g$idx), ]
g$year <- as.numeric(sub(" .*", "", g$label))
g$S <- as.numeric(sub(".*S", "", g$label))
par_fix <- fits_fix[[id_best_leg]]$OPT$par
s_fix <- par_fix[names(par_fix) == "S_pars"]
g$re_fixed <- fits_fix[[id_best_leg]]$RELGRP_err_v
g$total_fixed <- g$re_fixed + s_fix[2]*(g$water.yr == "wet") + s_fix[3]*(g$water.yr == "drought")

# release-group covariates from x.df (same rows as the design matrix)
e <- new.env()
load("data/HOR_CHP_mod_dat_ls.RData", envir = e)
x <- subset(e$x.df, !is.na(det.pred.qomt.rms))
x$label <- paste(x$year, x$rel.grp)
x$doy <- as.numeric(format(as.POSIXct(x$rel.date, tz = "UTC"), "%j"))

g$n <- as.vector(table(x$label)[g$label])
g$doy <- tapply(x$doy, x$label, mean)[g$label]
g$flength <- tapply(x$flength, x$label, mean)[g$label]
g[, c("label", "n", "water.yr", "doy", "flength", "re_fixed", "total_fixed")]

# pre-specified: release date, release order within the season, fish length, year
sapply(c("doy", "S", "flength", "year"), function(v)
  c(pearson = cor(g[[v]], g$total_fixed), spearman = cor(g[[v]], g$total_fixed, method = "spearman")))

# share of the group effect that lies between years
c(r2_year_factor = summary(lm(total_fixed ~ factor(year), data = g))$r.squared,
  r2_year_linear = summary(lm(total_fixed ~ year, data = g))$r.squared)

# scan of numeric x.df columns (group means); outcome-type columns are excluded
excl <- "^(tag$|tt\\.|tag\\.surv|CHP|BBR|TCJ|MAC|TRN|HOR\\.|det|form|year|day\\.|dawn\\.|dusk\\.|crep\\.|night\\.|censor|use\\.|line3|barrier|barPos|YrR|water|drought|xtrm|route|rel\\.grp|Y_|S$|X\\.)"
num <- names(x)[vapply(x, is.numeric, NA)]
cand <- num[!grepl(excl, num)]
gmean <- sapply(cand, function(v) tapply(x[[v]], x$label, mean, na.rm = TRUE)[g$label])
keep <- apply(gmean, 2, function(z) sum(is.finite(z)) >= 17 && sd(z, na.rm = TRUE) > 0)
gmean <- gmean[, keep]
c(n_candidates = length(cand), n_scanned = ncol(gmean))

rho <- apply(gmean, 2, function(z) cor(z, g$total_fixed, method = "spearman", use = "complete.obs"))
top <- head(sort(abs(rho), decreasing = TRUE), 15)
data.frame(covariate = names(top), spearman = rho[names(top)])

# permutation test of the largest |Spearman correlation| over all scanned columns
set.seed(5718)
perm_max <- replicate(3000, {
  y <- sample(g$total_fixed)
  max(abs(apply(gmean, 2, function(z) cor(z, y, method = "spearman", use = "complete.obs"))))
})
c(observed_max = max(abs(rho)), p_value = mean(perm_max >= max(abs(rho))), perm_95th = unname(quantile(perm_max, 0.95)))

ggplot(g, aes(x = year, y = total_fixed, label = paste0("S", S))) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_text() +
  labs(x = "release year", y = "group effect (WYT effect + random effect), fixed fit")
