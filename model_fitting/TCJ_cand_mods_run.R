# model_fitting/TCJ_cand_mods_run.R  (copy of Preliminaries/TCJ_all_mod_run.R, renamed for the model_fitting pipeline)
# Setting (define before sourcing to override): folder results are written to / read from.
if (!exists("out_dir")) out_dir <- "output"

library(glmmTMB)
library(ggplot2)
library(MuMIn)
require(parallel)
# require("corrplot")
require("car")
options(na.action="na.fail")
source("dAIC_threshold.R")

route_dfs_HOR_TCJ <- readRDS("data/route_dfs_HOR_TCJ.rds")
tcj.df <- route_dfs_HOR_TCJ$"tcj.df"
tcj.stn.df <- route_dfs_HOR_TCJ$"tcj.stn.df"

head(tcj.stn.df,2)

# taken from line 565
# barrier.exit.tcj.fac + log.VNS.ent.tcj.cms + CVP.exit.tcj.cms + SWP.exit.tcj.cms + flength + tod.ent.tcj.fac

all_vars <- c("barrier.exit.tcj.fac","log.VNS.ent.tcj","CVP.exit.tcj","SWP.exit.tcj","flength")
all_num_vars <- c("log.VNS.ent.tcj","CVP.exit.tcj","SWP.exit.tcj","flength")


unscl_ind <- match(all_num_vars,names(tcj.df))
scl_ind <- match(all_num_vars,names(tcj.stn.df))

head(tcj.df)
data.frame(
  attributes(scale(tcj.df[,unscl_ind]))["scaled:center"],
  attributes(scale(tcj.df[,unscl_ind]))["scaled:scale"])

plot(tcj.df$SWP.exit.tcj,tcj.df$CVP.exit.tcj)
cor(tcj.df$SWP.exit.tcj,tcj.df$CVP.exit.tcj)

TCJ_glmm_glob <-glmmTMB(TCJ.route.fac ~ barrier.exit.tcj.fac*(log.VNS.ent.tcj + CVP.exit.tcj + SWP.exit.tcj) + flength + (1|YrRel) ,
                         family="binomial", data=tcj.df)
bt=proc.time()
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
# clusterExport(clust, "tcj.df")
# clusterEvalQ(clust, library(glmmTMB))
# TCJ_glmm_glob.mod.dd<-MuMIn::dredge(TCJ_glmm_glob, fixed=c("barrier.exit.tcj.fac"), beta="none")#, cluster=clust, trace=2) # M model with barrier interaction effects
# stopCluster(clust)
source("R/glmmTMB_par.R")
TCJ_glmm_glob.mod.dd <- dredge_par_glmmTMB(TCJ_glmm_glob, fixed=c("cond(barrier.exit.tcj.fac)"), beta="none")

bt=proc.time()
# TCJ_mods_ls <- get.models(TCJ_glmm_glob.mod.dd, subset = NA)
TCJ_mods_ls <- get_models_par(TCJ_glmm_glob.mod.dd)
proc.time()-bt # 50s

TCJ_glmm_globZ <-glmmTMB(TCJ.route.fac ~ barrier.exit.tcj.fac*(log.VNS.ent.tcj + CVP.exit.tcj + SWP.exit.tcj) + flength + (1|YrRel) ,
                         family="binomial", data=tcj.stn.df)
bt=proc.time()
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
# clusterExport(clust, "tcj.stn.df")
# clusterEvalQ(clust, library(glmmTMB))
# TCJ_glmm_globZ.mod.dd<-MuMIn::dredge(TCJ_glmm_globZ, fixed=c("barrier.exit.tcj.fac"), beta="none")#, cluster=clust, trace=2) # M model with barrier interaction effects
# stopCluster(clust)
TCJ_glmm_globZ.mod.dd <- dredge_par_glmmTMB(TCJ_glmm_globZ, fixed=c("cond(barrier.exit.tcj.fac)"), beta="none")


bt=proc.time()
# TCJ_modsZ_ls <- get.models(TCJ_glmm_globZ.mod.dd, subset = NA)
TCJ_modsZ_ls <- get_models_par(TCJ_glmm_globZ.mod.dd)
proc.time()-bt # 50s

# length(TCJ_modsZ_ls)
# summary(TCJ_modsZ_ls[[5]])

get_aic_wts <- function(df, aic_column) {
  delta_aic <- df[[aic_column]] - min(df[[aic_column]], na.rm = TRUE)
  rel_likelihood <- exp(-0.5 * delta_aic)
  weights <- rel_likelihood / sum(rel_likelihood, na.rm = TRUE)
  return(weights)
}

formula_vZ=sapply(TCJ_modsZ_ls,function(xx) as.character(formula(xx)[3]))

TCJ_AIC_tabZ_unsort <- data.frame(
  formula=formula_vZ,
  AIC=sapply(TCJ_modsZ_ls,AIC),
  nll=sapply(TCJ_modsZ_ls,logLik)*-1,
  dAIC=sapply(TCJ_modsZ_ls,AIC)-min(sapply(TCJ_modsZ_ls,AIC)),
  k=sapply(TCJ_modsZ_ls, function(xx) length(xx$fit$par)))

TCJ_AIC_tabZ <- TCJ_AIC_tabZ_unsort[order(TCJ_AIC_tabZ_unsort$AIC),]
TCJ_AIC_tabZ$wt <- get_aic_wts(TCJ_AIC_tabZ,"AIC")
TCJ_AIC_tabZ$wt_cumul <- cumsum(TCJ_AIC_tabZ$wt)

saveRDS(as.list(environment()), file = file.path(out_dir, "TCJ_env_obs.rds"))

TCJ_aictabZ_d2 <- TCJ_AIC_tabZ |> dplyr::filter(dAIC<dAIC_threshold)
rownames(TCJ_aictabZ_d2)

TCJ_aictabZ_d2 <- TCJ_AIC_tabZ_unsort |> dplyr::filter(dAIC<dAIC_threshold)
TCJ_d2_modZ_ls <- TCJ_modsZ_ls[which(TCJ_AIC_tabZ$dAIC<dAIC_threshold)]





formula_v=sapply(TCJ_mods_ls,function(xx) as.character(formula(xx)[3]))

TCJ_AIC_tab_unsort <- data.frame(
  formula=formula_v,
  AIC=sapply(TCJ_mods_ls,AIC),
  nll=sapply(TCJ_mods_ls,logLik)*-1,
  dAIC=sapply(TCJ_mods_ls,AIC)-min(sapply(TCJ_mods_ls,AIC)),
  k=sapply(TCJ_mods_ls, function(xx) length(xx$fit$par)))

TCJ_AIC_tab <- TCJ_AIC_tab_unsort[order(TCJ_AIC_tab_unsort$AIC),]
TCJ_AIC_tab$wt <- get_aic_wts(TCJ_AIC_tab,"AIC")
TCJ_AIC_tab$wt_cumul <- cumsum(TCJ_AIC_tab$wt)

saveRDS(as.list(environment()), file = file.path(out_dir, "TCJ_env_obs.rds"))

TCJ_aictab_d2 <- TCJ_AIC_tab |> dplyr::filter(dAIC<dAIC_threshold)
rownames(TCJ_aictab_d2)

TCJ_aictab_d2 <- TCJ_AIC_tab_unsort |> dplyr::filter(dAIC<dAIC_threshold)
TCJ_d2_mod_ls <- TCJ_mods_ls[which(TCJ_AIC_tab$dAIC<dAIC_threshold)]


TCJ_mod_d2_ls <- list(
  "TCJ_d2_aictabZ"=TCJ_aictabZ_d2,
  "TCJ_d2_modsZ"=TCJ_d2_modZ_ls,
  "TCJ_d2_aictab"=TCJ_aictab_d2,
  "TCJ_d2_mods"=TCJ_d2_mod_ls
)

saveRDS(TCJ_mod_d2_ls,file.path(out_dir, "TCJ_d2_mods.rds"))
# TCJ_mod_d2_ls <- readRDS(file.path(out_dir, "TCJ_d2_mods.rds"))


# TCJ_mod_d2_ls <- readRDS("../CVPAS_beta/src/TCJ_d2_mods.rds")
