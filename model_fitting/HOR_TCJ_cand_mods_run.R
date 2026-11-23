# model_fitting/HOR_TCJ_cand_mods_run.R  (copy of Preliminaries/HOR_TCJ_all_mod_run.R, renamed for the model_fitting pipeline)
# Setting (define before sourcing to override): folder results are written to / read from.
if (!exists("out_dir")) out_dir <- "output"

library(glmmTMB)
library(ggplot2)
source("dAIC_threshold.R")
# library(dplyr)

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
# Model formulas (right-hand sides) used for both the unscaled and the Z-scaled model sets
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

# Fit the same formulas (plus a release-group random intercept) to the unscaled and the Z-scaled data. The calls are
# evaluated in the global environment so the fitted models do not carry this function's environment (and its data)
fit_HOR_TCJ_models <- function(data_name) {
  lapply(hor_tcj_rhs, function(rhs) {
    f <- as.formula(paste("detd ~", rhs, "+ (1|YrRel)"), env = globalenv())
    eval(bquote(glmmTMB(.(f), data = .(as.name(data_name)), family = binomial(link = "logit"), na.action = na.exclude)),
         globalenv())
  })
}

HOR_TCJ_glmm_ls <- fit_HOR_TCJ_models("kf.df")
HOR_TCJ_glmm_lsZ <- fit_HOR_TCJ_models("kf.dfZ")

HOR_TCJ_AIC_tabZ_unsort <- data.frame(
  AIC=sapply(HOR_TCJ_glmm_lsZ,AIC),
  nll=sapply(HOR_TCJ_glmm_lsZ,logLik)*-1,
  dAIC=sapply(HOR_TCJ_glmm_lsZ,AIC)-min(sapply(HOR_TCJ_glmm_lsZ,AIC)),
  k=sapply(HOR_TCJ_glmm_lsZ, function(xx) length(xx$fit$par)))

HOR_TCJ_AIC_tabZ <- HOR_TCJ_AIC_tabZ_unsort[order(HOR_TCJ_AIC_tabZ_unsort$AIC),]


get_aic_wts <- function(df, aic_column) {
  delta_aic <- df[[aic_column]] - min(df[[aic_column]], na.rm = TRUE)
  rel_likelihood <- exp(-0.5 * delta_aic)
  weights <- rel_likelihood / sum(rel_likelihood, na.rm = TRUE)
  return(weights)
}

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

summary(HOR_TCJ_glmm_ls[["WYT+flength+barrier*(VNS)"]])
summary(HOR_TCJ_glmm_ls[["WYT+flength+barrier+VNS"]])

ggplot(data=kf.df,aes(x=VNS.hor.5)) + facet_wrap(~WYT) + geom_histogram()
ggplot(data=kf.df,aes(x=log.VNS.hor.5)) + facet_wrap(~WYT) + geom_histogram()# $WYT
ggplot(data=kf.dfZ,aes(x=log.VNS.hor.5,fill=rel.grp)) + facet_wrap(~WYT) + geom_histogram()# $WYT
ggplot(data=kf.df,aes(x=Tmsd.hor.7dadm)) + facet_wrap(~WYT) + geom_histogram()# $WYT
ggplot(data=kf.df,aes(x=Tmsd.hor.7dadm,fill=rel.grp)) + facet_wrap(~WYT) + geom_histogram()# $WYT
ggplot(data=kf.dfZ,aes(x=Tmsd.hor.7dadm,fill=rel.grp)) + facet_wrap(~WYT) + geom_histogram()# $WYT

plot(kf.df$VNS.hor,kf.df$Tmsd.hor.7dadm)
plot(kf.df$log.VNS.hor.5,kf.df$Tmsd.hor.7dadm)
plot(kf.df$Tmsd.hor.7dadm,kf.df$log.VNS.hor.5)

length(HOR_TCJ_glmm_ls[[1]]$fit$par)
coef(HOR_TCJ_glmm_ls[[1]])
summary(HOR_TCJ_glmm_ls[[1]])$k


HOR_TCJ_glmm_glob=glmmTMB(detd~Year + barrier + Year:barrier+  (1|YrRel),
                          family=binomial(link="logit"),na.action=na.exclude,data=kf.df)

HOR_TCJ_glmm_glob=glmmTMB(detd~Year + barrier +  (1|YrRel),
                          family=binomial(link="logit"),na.action=na.exclude,data=kf.df)

### Newly added code

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


