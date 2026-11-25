# model_fitting/HOR_cand_mods_run.R  (copy of Preliminaries/HOR_all_mod_run.R, renamed for the model_fitting pipeline)
# Settings: model_fitting/config.R (pipeline_config), copied into the global environment.
if (!exists("pipeline_config")) source("model_fitting/config.R")
list2env(pipeline_config, envir = globalenv())  # out_dir, n_workers (dredge cluster size), conf_set_defs

library(glmmTMB)
library(ggplot2)
library(MuMIn)
require(parallel)
# require("corrplot")
require("car")
options(na.action="na.fail")
source("dAIC_threshold.R")

# library(dplyr)
# adapted from
# C:/repos/CVPAS_STH_app/Back end/model_averaging/HOR_TJC_all_mods_glmmTMB.R

# load("data/HOR_TCJ_survmod_data.RData")
route_dfs_HOR_TCJ <- readRDS("data/route_dfs_HOR_TCJ.rds")
hor.df <- route_dfs_HOR_TCJ$"hor.df"
hor.stn.df <- route_dfs_HOR_TCJ$"hor.stn.df"

unscl_ind <- match(c("log.VNS.ent.hor","CVP.exit.hor","SWP.exit.hor","flength"),names(hor.df))
scl_ind <- match(c("log.VNS.ent.hor","CVP.exit.hor","SWP.exit.hor","flength"),names(hor.stn.df))

head(hor.df)
data.frame(
  attributes(scale(hor.df[,unscl_ind]))["scaled:center"],
  attributes(scale(hor.df[,unscl_ind]))["scaled:scale"])

bt=proc.time()
HOR_glmm_globZ <-glmmTMB(HOR.route.fac ~ barrier.exit.hor.fac*(log.VNS.ent.hor + CVP.exit.hor + SWP.exit.hor) + flength + (1|YrRel) ,
                         family="binomial", data=hor.stn.df)
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
# clusterExport(clust, "hor.stn.df")
# clusterEvalQ(clust, library(glmmTMB))
# HOR_glmm_globZ.mod.dd<-MuMIn::dredge(HOR_glmm_globZ, fixed=c("barrier.exit.hor.fac"), beta="none")#, cluster=clust, trace=2) # M model with barrier interaction effects
# stopCluster(clust)
source("R/glmmTMB_par.R")
source("R/aic_weights.R")
HOR_glmm_globZ.mod.dd <- dredge_par_glmmTMB(HOR_glmm_globZ, fixed=c("cond(barrier.exit.hor.fac)"), beta="none", n_workers = n_workers)

# not standardizing outputs
HOR_glmm_glob <-glmmTMB(HOR.route.fac ~ barrier.exit.hor.fac*(log.VNS.ent.hor + CVP.exit.hor + SWP.exit.hor ) +flength + (1|YrRel),
                        family="binomial", data=hor.df)
proc.time()-bt 

bt=proc.time()
# cluster doesn't seem to speed it up much work
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
# clusterExport(clust, "hor.df")
# clusterEvalQ(clust, library(glmmTMB))
# HOR_glmm_glob.mod.dd<-MuMIn::dredge(HOR_glmm_glob, fixed=c("barrier.exit.hor.fac"), beta="none")#, cluster=clust, trace=2) # M model with barrier interaction effects
# stopCluster(clust)
HOR_glmm_glob.mod.dd <- dredge_par_glmmTMB(HOR_glmm_glob, fixed=c("cond(barrier.exit.hor.fac)"), beta="none", n_workers = n_workers)


proc.time()-bt 

# HOR_modsZ_d2_ls <- get.models(HOR_glmm_globZ.mod.dd, subset = delta < 2)
# HOR_modsZ_sum_d2_ls <- lapply(HOR_modsZ_d2_ls, summary)
# 
# HOR_mods_d2_ls <- get.models(HOR_glmm_glob.mod.dd, subset = delta < 2)
# HOR_mods_sum_d2_ls <- lapply(HOR_mods_d2_ls, summary)

############# #
# all models
############# #
bt=proc.time()
# HOR_modsZ_ls <- get.models(HOR_glmm_globZ.mod.dd, subset = NA)
HOR_modsZ_ls <- get_models_par(HOR_glmm_globZ.mod.dd, n_workers = n_workers)
# HOR_modsZ_sum_ls <- lapply(HOR_modsZ_ls, summary)
proc.time()-bt # 50s

bt=proc.time()
# HOR_mods_ls <- get.models(HOR_glmm_glob.mod.dd, subset = NA)
HOR_mods_ls <- get_models_par(HOR_glmm_glob.mod.dd, n_workers = n_workers)
# HOR_mods_sum_ls <- lapply(HOR_mods_ls, summary)
proc.time()-bt #1.5 mins

################## #
# standardized
################## #

HOR_AIC_tabZ_unsort <- data.frame(
  AIC=sapply(HOR_modsZ_ls,AIC),
  nll=sapply(HOR_modsZ_ls,logLik)*-1,
  dAIC=sapply(HOR_modsZ_ls,AIC)-min(sapply(HOR_modsZ_ls,AIC)),
  k=sapply(HOR_modsZ_ls, function(xx) length(xx$fit$par)))

HOR_AIC_tabZ <- HOR_AIC_tabZ_unsort[order(HOR_AIC_tabZ_unsort$AIC),]
HOR_AIC_tabZ$wt <- get_aic_wts(HOR_AIC_tabZ,"AIC")
HOR_AIC_tabZ$wt_cumul <- cumsum(HOR_AIC_tabZ$wt)

################## #
# unstandardized
################## #

HOR_AIC_tab_unsort <- data.frame(
  AIC=sapply(HOR_mods_ls,AIC),
  nll=sapply(HOR_mods_ls,logLik)*-1,
  dAIC=sapply(HOR_mods_ls,AIC)-min(sapply(HOR_mods_ls,AIC)),
  k=sapply(HOR_mods_ls, function(xx) length(xx$fit$par)))

HOR_AIC_tab <- HOR_AIC_tab_unsort[order(HOR_AIC_tab_unsort$AIC),]
HOR_AIC_tab$wt <- get_aic_wts(HOR_AIC_tab,"AIC")
HOR_AIC_tab$wt_cumul <- cumsum(HOR_AIC_tab$wt)

saveRDS(as.list(environment()), file = file.path(out_dir, "HOR_env_obs.rds"))
# HOR_env_obs <- readRDS(file.path(out_dir, "HOR_env_obs.rds"))
# list2env(HOR_env_obs,envir = .GlobalEnv)

HOR_aictab_d2 <- HOR_AIC_tab |> dplyr::filter(dAIC<dAIC_threshold)
rownames(HOR_aictab_d2)

HOR_d2_mods <- HOR_mods_ls[rownames(HOR_aictab_d2)]
length(HOR_d2_mods)

# HOR_mod_d2_ls <- readRDS("../CVPAS_beta/src/HOR_d2_mods.rds")
# saveRDS(HOR_mod_d2_ls,"../CVPAS_beta/src/HOR_d2_mods.rds")
HOR_aictabZ_d2 <- HOR_AIC_tabZ_unsort |> dplyr::filter(dAIC<dAIC_threshold)
# HOR_mod_d2_ls <- lapply(HOR_modsZ_ls[which(HOR_AIC_tab$dAIC<2)],function(x) names(x$frame[,-1]))
HOR_d2_modZ_ls <- HOR_modsZ_ls[which(HOR_AIC_tab$dAIC<dAIC_threshold)]
HOR_d2_mod_ls <- HOR_mods_ls[which(HOR_AIC_tab$dAIC<dAIC_threshold)]

# HOR_mod_d2_ls <- list(
# 
#   )

HOR_mod_d2_ls <- list(
  "HOR_d2_aictab"=HOR_aictabZ_d2,
  "HOR_d2_mods"=HOR_d2_mod_ls,
  "HOR_d2_aictabZ"=HOR_aictab_d2,
  "HOR_d2_modsZ"=HOR_d2_modZ_ls
)



saveRDS(HOR_mod_d2_ls,file.path(out_dir, "HOR_d2_mods_ls.rds"))

# saveRDS(HOR_mod_d2_ls, file = file.path(out_dir, "HOR_mod_d2_ls.rds"))
saveRDS(as.list(environment()), file = file.path(out_dir, "HOR_env_obs.rds"))

# ---- confidence sets for the definitions in pipeline_config$conf_set_defs -------- #
# Additional output: per-definition confidence sets and the meta set (the models
# meeting at least one definition), for the unscaled ("unscl") and Z-scaled ("scl") models.
source("R/conf_sets.R")
HOR_conf_sets <- list(
  "unscl" = compile_conf_sets(HOR_AIC_tab, HOR_mods_ls, pipeline_config$conf_set_defs),
  "scl"   = compile_conf_sets(HOR_AIC_tabZ, HOR_modsZ_ls, pipeline_config$conf_set_defs)
)
saveRDS(HOR_conf_sets, file.path(out_dir, "HOR_conf_sets.rds"))
