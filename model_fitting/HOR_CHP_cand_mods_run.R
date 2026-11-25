# model_fitting/HOR_CHP_cand_mods_run.R  (copy of Preliminaries/HOR_CHP_all_mod_run_REMOTE_posit.R, renamed for the model_fitting pipeline)
# Settings: model_fitting/config.R (pipeline_config), copied into the global environment.
if (!exists("pipeline_config")) source("model_fitting/config.R")
list2env(pipeline_config, envir = globalenv())  # n_workers (worker processes used by fit_nested_tmb_mods_par()), HOR_CHP_r_x_temp, CHP_wyt

library(TMB)
library(dplyr)
library(TMBhelp)

load("data/HOR_CHP_mod_dat_ls.RData")

# comp functions are defined in the data_prep script (run_data_prep <- FALSE loads the functions only)
run_data_prep <- FALSE
source("model_fitting/HOR_CHP_surv_mods_data_prep.R")

# # compile data from HOR-CHP survival model
# HOR_CHP_comp_ls <- HOR_CHP_comp()
# HOR_CHP_comp_ls_unscl <- HOR_CHP_comp(z_scale_vars = F)

source("R/HOR_CHP_designmat.R")
source("R/legacy_adapters.R")
source("R/HOR_CHP_comp_ALT.R")
HOR_CHP_precomp_ls <- HOR_CHP_precomp()
HOR_CHP_comp_ls <- HOR_CHP_modelprep(HOR_CHP_precomp_ls,
                                     design = legacy_HOR_CHP_design(HOR_CHP_precomp_ls, r_x_temp = HOR_CHP_r_x_temp, wyt = CHP_wyt))
# HOR_CHP_comp_ls_ALT <- HOR_CHP_comp_ALT()
# HOR_CHP_comp_ls_unscl <- HOR_CHP_comp_ALT(z_scale_vars = F)

# useful test
# identical(HOR_CHP_comp_ls_ALT,HOR_CHP_comp_ls)

# 
HOR_CHP_mod_dat_sep_ls <- readRDS("data/HOR_CHP_mod_dat_sep_ls.rds")
str(HOR_CHP_mod_dat_sep_ls$x.df)

TMB_data_baseline <- HOR_CHP_comp_ls$TMB_data_baseline
XX_in_w_int_WYT <- HOR_CHP_comp_ls$XX_in_w_int_WYT

nrow(TMB_data_baseline$Dat)

# HOR_CHP_comp_ls
head(XX_in_w_int_WYT)
allint_DMs <- HOR_CHP_comp_ls$allint_DMs

# superseded by HOR_CHP_designmat() in R/HOR_CHP_designmat.R
################################## #
# Creating design matrix indices
################################## #

# this is an extremely wide dataframe with all potential covariates
# xdf_tmp <- HOR_CHP_mod_dat_sep_ls$"x.df"
# xdf_tmp$WYT_ <- factor(xdf_tmp$water.yr,levels=c("dry","wet","drought"))

# xdf_tmp$Y_DUMMY = rbinom(n=nrow(xdf_tmp),size = 1,prob=0.5)

# drop rows with missing model variables (26 rows, all in Qomt.hor.1net) so the na.fail glm below works
# mains <- c("WYT_", "route.fac", "barrier.fac", "flength", "Tmsd.hor.7dadm", "log.VNS.hor.5", "SWP.hor.5", "CVP.hor.5", "Qomt.hor.1net")
# xdf_tmp <- xdf_tmp[complete.cases(xdf_tmp[, mains]), ]
# mm_tmp_raw <- model.matrix(~WYT_ + route.fac + barrier.fac + flength + Tmsd.hor.7dadm + log.VNS.hor.5 +  SWP.hor.5 +  CVP.hor.5 + Qomt.hor.1net +
#                              barrier.fac*(Tmsd.hor.7dadm + log.VNS.hor.5 +  SWP.hor.5 +  CVP.hor.5 + Qomt.hor.1net) +
#                              route.fac*(Tmsd.hor.7dadm + log.VNS.hor.5 +  SWP.hor.5 +  CVP.hor.5 + Qomt.hor.1net),
#                            data=xdf_tmp,method="model.frame")

# glm_tmp_glob <- glm(Y_DUMMY~WYT_ + route.fac + barrier.fac + flength + Tmsd.hor.7dadm + log.VNS.hor.5 +  SWP.hor.5 +  CVP.hor.5 + Qomt.hor.1net +
#                              barrier.fac*(Tmsd.hor.7dadm + log.VNS.hor.5 +  SWP.hor.5 +  CVP.hor.5 + Qomt.hor.1net) +
#                              route.fac*(Tmsd.hor.7dadm + log.VNS.hor.5 +  SWP.hor.5 +  CVP.hor.5 + Qomt.hor.1net),
#                            data=xdf_tmp,family="binomial", control = glm.control(maxit = 1), na.action = na.fail)#, control = list(iters=1))

# Corrected: fix all main effects (not only flength) so only the 10 interactions vary, and prepend intercept + 2 WYT_ dummies (1024 strings, nchar 21)
# model_tb_tmp <- MuMIn::dredge(glm_tmp_glob, fixed = mains, evaluate = FALSE)
# terms_list <- lapply(model_tb_tmp, function(call) attr(terms(eval(call$formula)), "term.labels"))
# all_terms <- attr(terms(glm_tmp_glob), "term.labels")
# dredge sometimes writes interactions in reversed order (e.g. "SWP.hor.5:route.fac"), so compare terms with sorted components
# sort_term <- function(x) vapply(strsplit(x, ":", fixed = TRUE), function(p) paste(sort(p), collapse = ":"), "")
# indicator_mat_raw <- t(sapply(terms_list, function(m_terms) as.integer(sort_term(all_terms) %in% sort_term(m_terms))))
# allint_DMs_ALT <- apply(cbind(1, indicator_mat_raw[, 1], indicator_mat_raw[, 1], indicator_mat_raw[, -1]),
#                         MARGIN = 1, FUN = paste, collapse = "")
# same model set as the stored design matrices (order differs)
# stopifnot(identical(sort(unname(allint_DMs_ALT)), sort(allint_DMs)))

# model_tb_tmp <- MuMIn::dredge(glm_tmp_glob,fixed=c("WYT_","flength"), evaluate = FALSE)
# model_tb_tmp <- MuMIn::dredge(glm_tmp_glob,fixed=c("flength"), evaluate = FALSE)
#
# terms_list <- lapply(model_tb_tmp, function(call) attr(terms(eval(call$formula)), "term.labels"))
# all_terms <- attr(terms(glm_tmp_glob), "term.labels")
# all_terms_mod <- gsub(all_terms,pattern="barrier.fac",replacement="B")
# all_terms_mod <- gsub(all_terms_mod,pattern="route.fac",replacement="R")
# all_terms_mod <- gsub(all_terms_mod,pattern=":",replacement=" x ")
#
# indicator_mat_raw <- t(sapply(terms_list, function(m_terms) as.integer(all_terms %in% m_terms)))
# indicator_mat_raw
# dim(indicator_mat_raw)
#
# indicator_mat <- indicator_mat_raw
# colnames(indicator_mat) <- all_terms
# dim(indicator_mat)
# indicator_mat
#
# keep_rows_ind <- which(indicator_mat_raw[,1]==1)
# allint_DMs_ALT <- as.character(apply(indicator_mat_raw[keep_rows_ind,],MARGIN=1,paste,collapse=""))

# allint_DMs
# allint_DMs_ALT

# str(model_tb_tmp)

# unique(mm_tmp_raw)
# unique((mm_tmp_raw>0)+0)

# attributes(HOR_CHP_comp_ls$XX_in)
# attributes(HOR_CHP_comp_ls_unscl$XX_in)
# 
# # extracting a few relevant vars
# allint_DMs <- HOR_CHP_comp_ls$allint_DMs

# ###################################################################################################### #
# Fitting all models
# ###################################################################################################### ## 

bt <- proc.time()
HOR_CHP_TMB_mod_fits_ls <- fit_nested_tmb_mods_par(
  S_design_matrix_in = allint_DMs,
  dir_in =file.path(getwd(), "TMB"),
  TMB_data_baseline_in=TMB_data_baseline,
  calc_SEs=FALSE,cpp_name_in = "HOR_CHP_global",
  n_workers = n_workers)
proc.time()-bt # 54.3 minutes with the sequential fit_nested_tmb_mods()
stopifnot(all(vapply(HOR_CHP_TMB_mod_fits_ls, function(x) isTRUE(x$converged), NA)))  # every fit must have converged

AIC_DF_tmp <- get_summ_tab_from_mod_ls(HOR_CHP_TMB_mod_fits_ls)

if(!file.exists(file.path(out_dir, "HOR_CHP_TMB_mod_fits_ls.rds"))){
  saveRDS(HOR_CHP_TMB_mod_fits_ls,file.path(out_dir, "HOR_CHP_TMB_mod_fits_ls.rds"))}


TMB_data_baseline <- HOR_CHP_comp_ls$TMB_data_baseline
XX_in_w_int_WYT <- HOR_CHP_comp_ls$XX_in_w_int_WYT


# assignig parameter names associated w each model parameter names
par_init_v_ls <- lapply(1:length(HOR_CHP_TMB_mod_fits_ls),function(ii) {
  HOR_CHP_TMB_mod_fits_ls[[ii]]$OPT$par})

par_mat1=do.call(rbind,strsplit(allint_DMs,split=""))=="1"
colnames(par_mat1)=colnames(XX_in_w_int_WYT)
par_nm_ind_ls <- apply(par_mat1,1,which,simplify = F)

HOR_CHP_TMB_mod_fits_ls[[1]]$OPT$par
names(HOR_CHP_TMB_mod_fits_ls[[1]]$OPT$par)

all_mod_parm_coef_ls <- lapply(1:length(HOR_CHP_TMB_mod_fits_ls),function(ii){
  aa=HOR_CHP_TMB_mod_fits_ls[[ii]]$OPT$par
  S_coef=aa[names(aa)=="S_pars"]
  other_pars=aa[names(aa)!="S_pars"]
  coeff_param=c(S_coef,other_pars)
  par_nm_tmp=c(names(par_nm_ind_ls[[ii]]),names(aa)[names(aa)!="S_pars"])
  c=data.frame(dmID=ii,
               dm=allint_DMs[ii],
               # par_nm=names(par_nm_ind_ls[[ii]]),
               par_nm=par_nm_tmp,
               S_parID=c(par_nm_ind_ls[[ii]],rep(NA,length(other_pars))),
               coeff_param)
  rownames(c)=NULL
  return(c)
  })

all_mod_parm_coefDF <- do.call(rbind,all_mod_parm_coef_ls)


HOR_CHP_TMB_all_mods <- list(
   "HOR_CHP_comp_ls"=HOR_CHP_comp_ls,
   "AIC_DF_full"=AIC_DF_tmp,
   "pt_estsDF"=all_mod_parm_coefDF,
   #elevating some objects for convenience
   "par_nm_ind_ls"=par_nm_ind_ls,
   "allint_DMs"=HOR_CHP_comp_ls$"allint_DMs",
   "TMB_data_baseline"=HOR_CHP_comp_ls$"TMB_data_baseline",
   "XX_in_w_int_WYT"=HOR_CHP_TMB_mod_fits_ls$"XX_in_w_int_WYT",
   "scaled_vars"=attributes(HOR_CHP_comp_ls$XX_in)$"scaled_vars",
   "center"=attributes(HOR_CHP_comp_ls$XX_in)$"center",
   "scale"=attributes(HOR_CHP_comp_ls$XX_in)$"scale")

HOR_CHP_TMB_all_mods


if(!file.exists(file.path(out_dir, "HOR_CHP_TMB_all_mods.rds"))){
  saveRDS(HOR_CHP_TMB_all_mods,file.path(out_dir, "HOR_CHP_TMB_all_mods.rds"))
}

# ---- confidence sets for the definitions in pipeline_config$conf_set_defs -------- #
# Table-level output (the fitted models stay in HOR_CHP_TMB_all_mods): per-definition
# sets and the meta set (models meeting at least one definition), optionally after
# the nesting rule (pipeline_config$nesting_rule).
source("R/conf_sets.R")
HOR_CHP_conf_sets <- compile_conf_sets(HOR_CHP_TMB_all_mods$"AIC_DF_full",
  mods_ls = NULL, defs = pipeline_config$conf_set_defs,
  nesting_rule = pipeline_config$nesting_rule)
saveRDS(HOR_CHP_conf_sets, file.path(out_dir, "HOR_CHP_conf_sets.rds"))



# start footer

#
# sink_TMB_R(cpp_name="HOR_CHP_global", dir_in=file.path(getwd(),"TMB"), silent = T, random = c("RELGRP_err"))
# OBJ_tmp1 <- TMB_model_setup(
#   cpp_name = "HOR_CHP_global",
#   dir_in =  file.path(getwd(),"TMB"),
#   data=HOR_CHP_TMB_all_mods$"TMB_data_baseline",
#   parameters=list(
#     "S_pars"=c(S_pars_tmp_g[1:21]),
#     "RELGRP_err"=rep(0,19),
#     "logSD_RELGRP"=S_pars_tmp_g[22],
#     "P_pars"=matrix(S_pars_tmp_g[23:30],nrow=length(HOR_CHP_TMB_all_mods$TMB_data_baseline$col_P11),
#                     ncol=HOR_CHP_TMB_all_mods$TMB_data_baseline$P_lc_n)),
#   random=c("RELGRP_err"),
#   recompile=F,
#   overwrite=F,
#   check_Eigen=F,
#   check_differen = T,
#   silent = T)
# 
# OBJ$par


# sink_TMB_R <- function(cpp_name,dir_in,silent,random){
#   sink( file.path(getwd(),dir_in,paste0(cpp_name,".R") ))
#   cat("
#       library( TMB )
#       setwd('",file.path(getwd(),dir_in),"')
#       dyn.load( TMB::dynlib('",cpp_name,"') )
#       load( '",cpp_name,"_inputs.RData' )
#       OBJ = MakeADFun(data=",cpp_name,"_inputs[['data']], parameters=",cpp_name,"_inputs[['parameters']],silent=",silent,",random='",random,"',inner.control = list(maxit=20000))
#       saveRDS(OBJ,'",cpp_name,"_OBJ.rds')\n",fill=FALSE,sep="")
#   sink()
# }
# 
# sink_TMB_R(cpp_name="HOR_CHP_global", dir_in="TMB", silent = T, random = c("RELGRP_err"))
# 
# source("TMB/HOR_CHP_global.R")
# sink_TMB_R
# setwd("..")
