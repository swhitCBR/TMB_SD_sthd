# dev/verify_HOR_TCJ_dredge.R
#
# Checks that the dredge()-based HOR-TCJ candidate fit (model_fitting/HOR_TCJ_cand_mods_run.R)
# reproduces the former hand-written candidate models (the hor_tcj_rhs list, kept at the
# top of that script): every manual model must appear in the dredged candidate set, with
# the same AIC and the same number of parameters.
#
# Usage, from the repository root, after running HOR_TCJ_cand_mods_run.R (typically with
# out_dir pointed at a scratch folder so output/ is not overwritten):
#   new_dir <- "C:/some/scratch/folder"   # where the new HOR_TCJ_all_mods.rds was written
#   source("dev/verify_HOR_TCJ_dredge.R")
# Leave new_dir unset to check the results in output/ itself.

if (!exists("new_dir")) new_dir <- "output"

old <- readRDS("output/HOR_TCJ_all_mods.rds")             # manual candidate fits (24 models)
new <- readRDS(file.path(new_dir, "HOR_TCJ_all_mods.rds")) # dredge() candidate fits

# match models by their term sets (labels of the fixed-effects terms of the formula).
# The components of an interaction label are sorted, because glmmTMB orders them
# differently depending on how the formula was constructed ("barrier:x" vs "x:barrier").
canonical_label <- function(l) paste(sort(strsplit(l, ":", fixed = TRUE)[[1]]), collapse = ":")
term_key <- function(f) paste(sort(vapply(attr(terms(f), "term.labels"),
                                          canonical_label, character(1))), collapse = " + ")
old_keys <- vapply(old$HOR_TCJ_glmm_ls, function(m) term_key(formula(m)), character(1))
new_keys <- vapply(new$HOR_TCJ_glmm_ls, function(m) term_key(formula(m)), character(1))

old_aic <- vapply(old$HOR_TCJ_glmm_ls, AIC, numeric(1))
old_k <- vapply(old$HOR_TCJ_glmm_ls, function(m) length(m$fit$par), numeric(1))
new_aic <- vapply(new$HOR_TCJ_glmm_ls, AIC, numeric(1))
new_k <- vapply(new$HOR_TCJ_glmm_ls, function(m) length(m$fit$par), numeric(1))

matched <- match(old_keys, new_keys)
if (anyNA(matched)) {
  cat("manual models NOT found in the dredged candidate set:\n")
  print(names(old_keys)[is.na(matched)])
} else {
  cat("all", length(old_keys), "manual models found in the dredged candidate set\n")
}

cmp <- data.frame(
  manual = names(old_keys),
  AIC_manual = old_aic,
  AIC_dredge = new_aic[matched],
  k_manual = old_k,
  k_dredge = new_k[matched],
  row.names = NULL)
cmp$AIC_diff <- cmp$AIC_dredge - cmp$AIC_manual
print(cmp, digits = 8)

cat("\nlargest |AIC difference|:", max(abs(cmp$AIC_diff), na.rm = TRUE), "\n")
cat("k (parameter count) differences:", sum(cmp$k_manual != cmp$k_dredge, na.rm = TRUE), "\n")
cat("manual AIC ranks preserved:",
    identical(order(cmp$AIC_manual), order(cmp$AIC_dredge)), "\n")
cat("dredged candidate set size:", length(new_keys),
    "(manual list:", length(old_keys), ")\n")
