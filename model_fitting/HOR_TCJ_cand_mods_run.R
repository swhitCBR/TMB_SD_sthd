# model_fitting/HOR_TCJ_cand_mods_run.R  (copy of Preliminaries/HOR_TCJ_all_mod_run.R, renamed for the model_fitting pipeline)
# Settings: model_fitting/config.R (pipeline_config), copied into the global environment.
if (!exists("pipeline_config")) source("model_fitting/config.R")
list2env(pipeline_config, envir = globalenv())  # out_dir, n_workers (dredge cluster size), conf_set_defs

library(glmmTMB)
library(ggplot2)
library(MuMIn)
source("dAIC_threshold.R")
# dredge() needs na.fail (there are no NAs in the model covariates)
options(na.action = "na.fail")

# adapted from
# C:/repos/CVPAS_STH_app/Back end/model_averaging/HOR_TJC_all_mods_glmmTMB.R

load("data/HOR_TCJ_survmod_data.RData")

# Z scaling data

c("detd","Year","barrier","YrRel","log.VNS.hor.5","flength","Tmsd.hor.7dadm","CVP.hor.5","SWP.hor.5") %in% names(kf.df)
# (subsetting, scaling numeric values (retainng center[mean] and scale[sd] values))
kd.df_unscl_num_mat <- kf.df[,c("log.VNS.hor.5","flength","Tmsd.hor.7dadm","CVP.hor.5","SWP.hor.5")]
kd.df_scl_cols <- base::scale(kd.df_unscl_num_mat,center = T,scale = T)
kf.dfZ <- cbind(kf.df[,-which(colnames(kf.df) %in% colnames(kd.df_scl_cols))],kd.df_scl_cols)
# attributes(kd.df_scl_cols) # shape and scale information

kf.dfZ[1:6,c(1:5,524:529)]
dim(kf.dfZ)

kd.df_unscl_num_mat

# Former hand-written candidate list (right-hand sides), kept for reference and for
# dev/verify_HOR_TCJ_dredge.R, which checks that the dredge() candidate set below
# reproduces these models. Every term and interaction of these formulas is in the
# global model, so each of them is among the dredged candidates.
hor_tcj_rhs <- c(
  "int_only" = "1",
  "year" = "Year",
  "VNS" = "log.VNS.hor.5",
  "barrier" = "barrier",
  "WYT" = "WYT",
  "WYT+barrier" = "WYT + barrier",
  "WYT+barrier+temp" = "WYT + barrier + Tmsd.hor.7dadm",
  "WYT+barrier+VNS" = "WYT + barrier + log.VNS.hor.5",
  "WYT+flength+barrier" = "WYT + flength + barrier",
  "WYT+flength+barrier+VNS" = "WYT + flength + barrier + log.VNS.hor.5",
  "WYT+flength+barrier+temp" = "WYT + flength + barrier + Tmsd.hor.7dadm",
  "WYT+flength+barrier+VNS+CVP" = "WYT + flength + barrier + CVP.hor.5 + log.VNS.hor.5",
  "WYT+flength+barrier+VNS+SWP" = "WYT + flength + barrier + SWP.hor.5 + log.VNS.hor.5",
  "WYT+flength+barrier+VNS+temp" = "WYT + flength + barrier + log.VNS.hor.5 + Tmsd.hor.7dadm",
  "WYT+flength+barrier+temp+CVP" = "WYT + flength + barrier + Tmsd.hor.7dadm + CVP.hor.5",
  "WYT+flength+barrier+temp+SWP" = "WYT + flength + barrier + Tmsd.hor.7dadm + SWP.hor.5",
  "WYT+flength+barrier*(temp)" = "WYT + flength + barrier + Tmsd.hor.7dadm + barrier:Tmsd.hor.7dadm",
  "WYT+flength+barrier*(VNS)" = "WYT + flength + log.VNS.hor.5 + barrier + barrier:log.VNS.hor.5",
  "WYT+flength+barrier*(temp) + VNS" = "WYT + flength + barrier + log.VNS.hor.5 + Tmsd.hor.7dadm + barrier:Tmsd.hor.7dadm",
  "WYT+flength+barrier*(VNS) + temp" = "WYT + flength + log.VNS.hor.5 + Tmsd.hor.7dadm + barrier + barrier:log.VNS.hor.5",
  "WYT+flength+barrier*(VNS+temp)" = "WYT + flength + log.VNS.hor.5 + Tmsd.hor.7dadm + barrier + barrier:log.VNS.hor.5 + barrier:Tmsd.hor.7dadm",
  "WYT+flength+barrier*(temp)+CVP" = "WYT + flength + barrier + Tmsd.hor.7dadm + barrier:Tmsd.hor.7dadm + CVP.hor.5",
  "WYT+flength+barrier*(VNS)+SWP" = "WYT + flength + log.VNS.hor.5 + barrier + barrier:log.VNS.hor.5 + SWP.hor.5",
  "WYT+flength+barrier*(temp+CVP)" = "WYT + flength + barrier + Tmsd.hor.7dadm + CVP.hor.5 + barrier:Tmsd.hor.7dadm + barrier:CVP.hor.5",
  "WYT+flength+barrier*(temp+SWP)" = "WYT + flength + barrier + Tmsd.hor.7dadm + SWP.hor.5 + barrier:Tmsd.hor.7dadm + barrier:SWP.hor.5")

#### old manual fits (replaced by the dredge() implementation below) ####
# source("R/glmmTMB_par.R")
# HOR_TCJ_glmm_ls <- fit_glmmTMB_par(hor_tcj_rhs, "kf.df", response = "detd", random = "(1|YrRel)",
#                                         family = binomial(link = "logit"), na.action = na.exclude, n_workers = n_workers)
# HOR_TCJ_glmm_lsZ <- fit_glmmTMB_par(hor_tcj_rhs, "kf.dfZ", response = "detd", random = "(1|YrRel)",
#                                          family = binomial(link = "logit"), na.action = na.exclude, n_workers = n_workers)

#### new dredge() implementation ####
# The candidate set is all subsets of a global model (as for the HOR and TCJ routing
# models), instead of the hand-written hor_tcj_rhs list. The global model contains every
# main covariate and every barrier interaction of the manual list; dredge() keeps an
# interaction only together with its main effects, so the candidate set is a superset of
# the manual list. Models with both Year and WYT are excluded (subset = !(`cond(Year)` && `cond(WYT)`)):
# WYT is collinear with Year (its dummies drop from rank-deficient fits) and the manual
# list also never mixes the two.
source("R/glmmTMB_par.R")
source("R/aic_weights.R")
HOR_TCJ_glob_form <- detd ~ Year + WYT + barrier + flength + log.VNS.hor.5 +
  Tmsd.hor.7dadm + CVP.hor.5 + SWP.hor.5 +
  barrier:Tmsd.hor.7dadm + barrier:log.VNS.hor.5 + barrier:CVP.hor.5 + barrier:SWP.hor.5 +
  (1|YrRel)

bt <- proc.time()
HOR_TCJ_glmm_glob <- glmmTMB(HOR_TCJ_glob_form, family = binomial(link = "logit"),
                             na.action = na.fail, data = kf.df)
HOR_TCJ_glob.mod.dd <- dredge_par_glmmTMB(HOR_TCJ_glmm_glob, beta = "none",
                                          subset = !(`cond(Year)` && `cond(WYT)`), n_workers = n_workers)
HOR_TCJ_glmm_ls <- get_models_par(HOR_TCJ_glob.mod.dd, n_workers = n_workers)
proc.time() - bt

bt <- proc.time()
HOR_TCJ_glmm_globZ <- glmmTMB(HOR_TCJ_glob_form, family = binomial(link = "logit"),
                              na.action = na.fail, data = kf.dfZ)
HOR_TCJ_globZ.mod.dd <- dredge_par_glmmTMB(HOR_TCJ_glmm_globZ, beta = "none",
                                           subset = !(`cond(Year)` && `cond(WYT)`), n_workers = n_workers)
HOR_TCJ_glmm_lsZ <- get_models_par(HOR_TCJ_globZ.mod.dd, n_workers = n_workers)
proc.time() - bt

HOR_TCJ_AIC_tabZ_unsort <- data.frame(
  AIC=sapply(HOR_TCJ_glmm_lsZ,AIC),
  nll=sapply(HOR_TCJ_glmm_lsZ,logLik)*-1,
  dAIC=sapply(HOR_TCJ_glmm_lsZ,AIC)-min(sapply(HOR_TCJ_glmm_lsZ,AIC)),
  k=sapply(HOR_TCJ_glmm_lsZ, function(xx) length(xx$fit$par)))

HOR_TCJ_AIC_tabZ <- HOR_TCJ_AIC_tabZ_unsort[order(HOR_TCJ_AIC_tabZ_unsort$AIC),]


HOR_TCJ_AIC_tabZ$wt <- get_aic_wts(HOR_TCJ_AIC_tabZ,"AIC")
HOR_TCJ_AIC_tabZ$wt_cumul <- cumsum(HOR_TCJ_AIC_tabZ$wt)


########################## #
# Alternative criteria
########################## #

# 95% of the weight of the full set
HOR_TCJ_AIC_tabZ[1:min(which(HOR_TCJ_AIC_tabZ$wt_cumul >0.95)),]
# within 4 AIC
HOR_TCJ_AIC_tabZ[which(HOR_TCJ_AIC_tabZ$dAIC <4),]
# within 6 AIC
HOR_TCJ_AIC_tabZ[which(HOR_TCJ_AIC_tabZ$dAIC <6),]

summary(HOR_TCJ_glmm_ls[[1]])

# covariate diagnostics (moved to R/plots.R)
source("R/plots.R")
plot_HOR_TCJ_diagnostics(kf.df, kf.dfZ)

length(HOR_TCJ_glmm_ls[[1]]$fit$par)
coef(HOR_TCJ_glmm_ls[[1]])
summary(HOR_TCJ_glmm_ls[[1]])$k


### AIC table of the unscaled fits
HOR_TCJ_aictab <- data.frame(
  AIC=sapply(HOR_TCJ_glmm_ls,AIC),
  nll=sapply(HOR_TCJ_glmm_ls,logLik)*-1,
  dAIC=sapply(HOR_TCJ_glmm_ls,AIC)-min(sapply(HOR_TCJ_glmm_ls,AIC)),
  k=sapply(HOR_TCJ_glmm_ls, function(xx) length(xx$fit$par)))

HOR_TCJ_aictab<- HOR_TCJ_aictab[order(HOR_TCJ_aictab$AIC),]
HOR_TCJ_aictab

HOR_TCJ_mod_all_ls <- list(
  "HOR_TCJ_aictab"=HOR_TCJ_aictab,
  "HOR_TCJ_glmm_ls"=HOR_TCJ_glmm_ls)

saveRDS(HOR_TCJ_mod_all_ls,file.path(out_dir, "HOR_TCJ_all_mods.rds"))

HOR_TCJ_aictab_d2 <- HOR_TCJ_aictab |> dplyr::filter(dAIC<dAIC_threshold)
rownames(HOR_TCJ_aictab_d2)

HOR_TCJ_d2_mods <- HOR_TCJ_glmm_ls[rownames(HOR_TCJ_aictab_d2)]

HOR_TCJ_mod_d2_ls <- list(
  "HOR_TCJ_aictab"=HOR_TCJ_aictab,
  "HOR_TCJ_d2_mods"=HOR_TCJ_d2_mods)


# HOR_TCJ_mod_d2_ls <- readRDS("../CVPAS_beta/src/HOR_TCJ_d2_mods.rds")
# saveRDS(HOR_TCJ_mod_d2_ls,"../CVPAS_beta/src/HOR_TCJ_d2_mods.rds")
saveRDS(HOR_TCJ_mod_d2_ls,file.path(out_dir, "HOR_TCJ_d2_mods.rds"))

# ---- confidence sets for the definitions in pipeline_config$conf_set_defs -------- #
# Additional output: per-definition confidence sets and the meta set (the models
# meeting at least one definition), for the unscaled ("unscl") and Z-scaled ("scl") models.
source("R/conf_sets.R")
HOR_TCJ_conf_sets <- list(
  "unscl" = compile_conf_sets(HOR_TCJ_aictab, HOR_TCJ_glmm_ls, pipeline_config$conf_set_defs),
  "scl"   = compile_conf_sets(HOR_TCJ_AIC_tabZ, HOR_TCJ_glmm_lsZ, pipeline_config$conf_set_defs)
)
saveRDS(HOR_TCJ_conf_sets, file.path(out_dir, "HOR_TCJ_conf_sets.rds"))
