# model_fitting/TCJ_CHP_cand_mods_run_init.R  (copy of Preliminaries/TCJ_CHP_all_mod_run_init.R, renamed for the model_fitting pipeline)
# Setting (define before sourcing to override): folder results are written to / read from.
if (!exists("out_dir")) out_dir <- "output"
if (!exists("n_workers") || !exists("CHP_wyt")) source("model_fitting/config.R")   # n_workers (set in config.R): worker processes used by fit_nested_tmb_mods_par()

####################################### #
# data prep example starts here
####################################### #

# taken from "C:/repos/CVPAS_STH_app/Preliminaries/TCJ_CHP_match_likelihood_for_toward_RE_model2.R"

# library(TMB)
# load("tmp_data/tmp_img4FULL.Rdata")
# lvec.tcj_ls <- readRDS("tmp_data/lvec.tcj_ls.rds")
# devtools::load_all("../CVhelp")
# 
# objs_for_surv_TCJ_CHP_mods <- c("x.df","lvec.tcj_ls","sdat.det.common2")
# save(list=objs_for_surv_TCJ_CHP_mods,file="C:/repos/CVPAS_STH_app/data/TCJ_CHP_mod_dat_ls.RData")

###################### #
# cand_mods_run code starts here
###################### #

library(TMB)
library(dplyr)

library(TMBhelp)
# devtools::load_all("../CVhelp")  # no longer needed: TCJ_CHP_comp() is in TCJ_CHP_surv_mods_data_prep.R

# lvec.tcj_ls <- readRDS("C:/repos/CVPAS-steelhead/tmp_data/lvec.tcj_ls.rds")
load("data/TCJ_CHP_mod_dat_ls.RData")

source("R/HOR_CHP_designmat.R")
source("R/TCJ_CHP_designmat.R")
source("R/legacy_adapters.R")
source("R/TCJ_CHP_comp_ALT.R")
TCJ_CHP_precomp_ls <- TCJ_CHP_precomp()
TCJ_CHP_comp_ls <- TCJ_CHP_modelprep(TCJ_CHP_precomp_ls, design = legacy_TCJ_CHP_design(TCJ_CHP_precomp_ls, wyt = CHP_wyt))

# comp function is defined in the data_prep script (run_data_prep <- FALSE loads the function only)
# run_data_prep <- FALSE
# source("model_fitting/TCJ_CHP_surv_mods_data_prep.R")

# TCJ_CHP_comp_ls <- TCJ_CHP_comp(
#   # XX_in
#   # ,
#   sdat.det.in = sdat.det.common2,
#   lvec.tcj_ls_in = lvec.tcj_ls,
#   x.df = x.df
# )

head(TCJ_CHP_comp_ls$XX_in_add)
head(TCJ_CHP_comp_ls$XX_in_w_int_WYT)
head(TCJ_CHP_comp_ls$TMB_data_baseline$XX_s)

dim(TCJ_CHP_comp_ls$XX_in_w_int_WYT)

# # compile data from HOR-CHP survival model
# HOR_CHP_comp_ls <- CVhelp::HOR_CHP_comp()
# HOR_CHP_comp_ls_unscl <- CVhelp::HOR_CHP_comp(z_scale_vars = F)

TMB_data_baseline <- TCJ_CHP_comp_ls$TMB_data_baseline
XX_in_w_int_WYT <- TCJ_CHP_comp_ls$XX_in_w_int_WYT
allint_DMs <- TCJ_CHP_comp_ls$allint_DMs

# TCJ_CHP_comp_ls$TMB_data_baseline$s
# P_mat_base <- matrix(0, ncol = 4)
# 
# OBJ <- TMB_model_setup(cpp_name = "HOR_CHP_global",
#                        dir_in = file.path(getwd(),"TMB"),
#                        data=TCJ_CHP_comp_ls$TMB_data_baseline,
#                        parameters=list(
#                          "S_pars"=c(1,rep(0,ncol(TCJ_CHP_comp_ls$TMB_data_baseline$XX_s)-1)),
#                          "logSD_RELGRP"=-2,
#                          "RELGRP_err"=rep(0,TCJ_CHP_comp_ls$TMB_data_baseline$n_relgrps),
#                          "P_pars"=P_mat_base
#                        ),
#                        random="RELGRP_err",
#                        recompile=F,
#                        overwrite=F,
#                        check_Eigen=F,
#                        check_differen = T,
#                        silent = T,hessian = T)
# OBJ$fn()
# OBJ$par

# ###################################################################################################### #
# Fitting all models
# ###################################################################################################### ##

bt <- proc.time()
TCJ_CHP_TMB_mod_fits_ls <- fit_nested_tmb_mods_par(
  S_design_matrix_in = allint_DMs,
  dir_in = file.path(getwd(), "TMB"),
  TMB_data_baseline_in = TMB_data_baseline,
  calc_SEs = FALSE,
  cpp_name_in = "HOR_CHP_global",
  n_workers = n_workers
)
proc.time() - bt # 40 seconds with the sequential fit_nested_tmb_mods()
stopifnot(all(vapply(TCJ_CHP_TMB_mod_fits_ls, function(x) isTRUE(x$converged), NA)))  # every fit must have converged

AIC_DF_tmp <- get_summ_tab_from_mod_ls(TCJ_CHP_TMB_mod_fits_ls)

if (!file.exists(file.path(out_dir, "TCJ_CHP_TMB_mod_fits_ls.rds"))) {
  saveRDS(TCJ_CHP_TMB_mod_fits_ls,
          file.path(out_dir, "TCJ_CHP_TMB_mod_fits_ls.rds"))
}

# TCJ_CHP_TMB_mod_fits_ls <- readRDS(file.path(out_dir, "TCJ_CHP_TMB_mod_fits_ls.rds"))


TMB_data_baseline <- TCJ_CHP_comp_ls$TMB_data_baseline
XX_in_w_int_WYT <- TCJ_CHP_comp_ls$XX_in_w_int_WYT

# assigning parameter names associated w each model parameter names
par_init_v_ls <- lapply(1:length(TCJ_CHP_TMB_mod_fits_ls), function(ii) {
  TCJ_CHP_TMB_mod_fits_ls[[ii]]$OPT$par})

par_mat1 = do.call(rbind, strsplit(allint_DMs, split = "")) == "1"
colnames(par_mat1) = colnames(XX_in_w_int_WYT)
par_nm_ind_ls <- apply(par_mat1, 1, which, simplify = F)

TCJ_CHP_TMB_mod_fits_ls[[1]]$OPT$par
names(TCJ_CHP_TMB_mod_fits_ls[[1]]$OPT$par)

ii=1

all_mod_parm_coef_ls <- lapply(1:length(TCJ_CHP_TMB_mod_fits_ls), function(ii) {
  aa = TCJ_CHP_TMB_mod_fits_ls[[ii]]$OPT$par
  S_coef = aa[names(aa) == "S_pars"]
  other_pars = aa[names(aa) != "S_pars"]
  coeff_param = c(S_coef, other_pars)
  par_nm_tmp = c(names(par_nm_ind_ls[[ii]]), names(aa)[names(aa) != "S_pars"])
  c = data.frame(
    dmID = ii,
    dm = allint_DMs[ii],
    # par_nm=names(par_nm_ind_ls[[ii]]),
    par_nm = par_nm_tmp,
    S_parID = c(par_nm_ind_ls[[ii]], rep(NA, length(other_pars))),
    coeff_param
  )
  rownames(c) = NULL
  return(c)
})

S_coef_ls <- lapply(1:length(TCJ_CHP_TMB_mod_fits_ls), function(ii) {
  # aa=HOR_CHP_TMB_mod_fits_ls[[ii]]$"est_tab"
  # S_coef=aa$"Estimate_SE"[aa$"Parameter"=="S_pars"]
  aa = TCJ_CHP_TMB_mod_fits_ls[[ii]]$OPT$par
  S_coef = aa[names(aa) == "S_pars"]
  c = data.frame(
    dmID = ii,
    dm = allint_DMs[ii],
    par_nm = names(par_nm_ind_ls[[ii]]),
    S_parID = par_nm_ind_ls[[ii]],
    S_coef
  )
  rownames(c) = NULL
  return(c)
})



all_mod_parm_coefDF <- do.call(rbind,S_coef_ls)


TCJ_CHP_TMB_all_mods <- list(
  "TCJ_CHP_comp_ls" = TCJ_CHP_comp_ls,
  "AIC_DF_full" = AIC_DF_tmp,
  "pt_estsDF" = all_mod_parm_coefDF,
  #elevating some objects for convenience
  "par_nm_ind_ls" = par_nm_ind_ls,
  "allint_DMs" =allint_DMs,# TCJ_CHP_comp_ls$"allint_DMs", # DIFFERENCE
  "TMB_data_baseline" = TCJ_CHP_comp_ls$"TMB_data_baseline",
  "XX_in_w_int_WYT" = TCJ_CHP_TMB_mod_fits_ls$"XX_in_w_int_WYT"#,
  # "scaled_vars"=attributes(TCJ_CHP_comp_ls$XX_in)$"scaled_vars",
  # "center"=attributes(TCJ_CHP_comp_ls$XX_in)$"center",
  # "scale"=attributes(TCJ_CHP_comp_ls$XX_in)$"scale"
)

TCJ_CHP_TMB_all_mods$AIC_DF_full


if (!file.exists(file.path(out_dir, "TCJ_CHP_TMB_all_mods.rds"))) {
  saveRDS(TCJ_CHP_TMB_all_mods, file.path(out_dir, "TCJ_CHP_TMB_all_mods.rds"))
}



# ##################################################
# Demo with first 3 models
# ##################################################
# 
# allint_DMs <- apply(cbind(rep(c("11111111111"), 2^5), expand.grid(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1)))
#                     , 1, paste0, collapse = "")
# DM_mat<- do.call(rbind,lapply(1:length(allint_DMs),function(x) as.numeric(strsplit(allint_DMs[x],split="")[[1]]) ))
# colnames(DM_mat) <- colnames(TMB_data_baseline$XX_s)
# 
# head(TMB_data_baseline$XX_s)

# allint_DMs <- apply(cbind(rep(c("1111111111"), 2^3), expand.grid(c(0, 1), c(0, 1), c(0, 1)))
#                     , 1, paste0, collapse = "")
# allint_DMs<- c("1111000000000",allint_DMs)
# allint_DMs<- c("1111000000000",allint_DMs)

# 2^9=512
# allint_DMs_raw <- apply(cbind(rep(c("1111"), 2^9), expand.grid(c(0, 1), c(0, 1), c(0, 1),c(0, 1), c(0, 1), c(0, 1),c(0, 1), c(0, 1), c(0, 1)))
#                     , 1, paste0, collapse = "")
# DM_mat_raw<- do.call(rbind,lapply(1:length(allint_DMs_raw),function(x) as.numeric(strsplit(allint_DMs_raw[x],split="")[[1]]) ))
# colnames(DM_mat_raw) <- colnames(TMB_data_baseline$XX_s)
# 
# bad_ind <- which(DM_mat_raw[,"Qomt.tcj.1net"]==0 & DM_mat_raw[,"R_x_OMT"]==1 | DM_mat_raw[,"CVP.tcj.4"]==0 & DM_mat_raw[,"R_x_CVP"]==1 | DM_mat_raw[,"SWP.tcj.4"]==0 & DM_mat_raw[,"R_x_SWP"]==1)
# allint_DMs_raw[bad_ind]
# DM_mat <- DM_mat_raw[-c(bad_ind),]
# allint_DMs <- allint_DMs_raw[-c(bad_ind)]

# which(is.na(TMB_data_baseline$XX_s),arr.ind = T)

# TMB:::getUserDLL()
# 
# TCJ_CHP_frst3_mod_ls <- fit_nested_tmb_mods(
#   S_design_matrix_in = allint_DMs[c(1:3)],
#   dir_in = file.path(getwd(), "TMB"),
#   TMB_data_baseline_in = TMB_data_baseline,
#   calc_SEs = FALSE,
#   cpp_name_in = "HOR_CHP_global"
# )
# 
# # TCJ_CHP_frst3_mod_ls <- fit_nested_tmb_mods(
# #   S_design_matrix_in = allint_DMs[c(1:3)],
# #   dir_in =file.path(getwd(), "TMB"),
# #   TMB_data_baseline_in=HOR_CHP_comp_ls_unscl$"TMB_data_baseline",
# #   calc_SEs=FALSE,cpp_name_in = "HOR_CHP_global")
# 
# TCJ_CHP_frst3_mod_ls_SE <- fit_nested_tmb_mods(
#   S_design_matrix_in = allint_DMs[c(1:3)],
#   dir_in = file.path(getwd(), "TMB"),
#   TMB_data_baseline_in = TMB_data_baseline,
#   calc_SEs = TRUE,
#   cpp_name_in = "HOR_CHP_global"
# )
# 
# devtools::load_all("../TMBhelp")
# AIC_DF_frst3_DF <- get_summ_tab_from_mod_ls(TCJ_CHP_frst3_mod_ls)
