# model_fitting/HOR_CHP_d2_AIC_unscl_mod_ls.R  (copy of Preliminaries/HOR_CHP_d2_AIC_unscl_mod_ls.R, renamed for the model_fitting pipeline)
# Settings: model_fitting/config.R (pipeline_config), copied into the global environment.
if (!exists("pipeline_config")) source("model_fitting/config.R")
list2env(pipeline_config, envir = globalenv())  # HOR_CHP_r_x_temp, CHP_wyt

library(TMB)
library(dplyr)

library(TMBhelp)
# devtools::load_all("../CVhelp")  # no longer needed: HOR_CHP_comp() is in HOR_CHP_surv_mods_data_prep.R

# load("data/HOR_CHP_mod_dat_ls.RData")

# compile data from HOR-CHP survival model
# HOR_CHP_comp_ls <- CVhelp::HOR_CHP_comp()
source("R/HOR_CHP_designmat.R")
source("R/legacy_adapters.R")
source("R/HOR_CHP_comp_ALT.R")
HOR_CHP_precomp_ls_unscl <- HOR_CHP_precomp(z_scale_vars = F)
HOR_CHP_comp_ls_unscl <- HOR_CHP_modelprep(HOR_CHP_precomp_ls_unscl,
                                           design = legacy_HOR_CHP_design(HOR_CHP_precomp_ls_unscl, r_x_temp = HOR_CHP_r_x_temp, wyt = CHP_wyt))
# HOR_CHP_comp_ls_unscl <- HOR_CHP_comp(z_scale_vars = F)

# attributes(HOR_CHP_comp_ls$XX_in)
# attributes(HOR_CHP_comp_ls_unscl$XX_in)

# extracting a few relevant vars
# allint_DMs <- HOR_CHP_comp_ls$allint_DMs
# TMB_data_baseline <- HOR_CHP_comp_ls$TMB_data_baseline
XX_in_w_int_WYT <- HOR_CHP_comp_ls_unscl$XX_in_w_int_WYT

S_conf_ls <- readRDS(file.path(out_dir, "HOR_CHP_S_conf_ls.rds"))
S_conf_ls$AIC_DF_d2$dm

# placeholder when fitting
HOR_CHP_comp_ls_unscl$"TMB_data_baseline"$"XX_pred_mat" <- HOR_CHP_comp_ls_unscl$TMB_data_baseline$XX_s[1:10,]


bt=proc.time()
HOR_CHP_d2_mod_ls_unscl <- fit_nested_tmb_mods(
  S_design_matrix_in = S_conf_ls$AIC_DF_d2$"dm",
  dir_in =file.path(getwd(), "TMB"),
  TMB_data_baseline_in=HOR_CHP_comp_ls_unscl$"TMB_data_baseline",
  calc_SEs=TRUE,
  # cpp_name_in = "HOR_CHP_global", # old version that doesn't have "XX_pred_mat" for easy prediction
  cpp_name_in = "HOR_CHP_global_pred"
  )
proc.time()-bt

HOR_CHP_TMB_mod_fits_unscl_ls <- HOR_CHP_d2_mod_ls_unscl
# HOR_CHP_TMB_mod_fits_unscl_ls$`111111111110100001011`$est_tab

# bt <- proc.time()
# HOR_CHP_TMB_mod_fits_unscl_ls <- fit_nested_tmb_mods(
#   S_design_matrix_in = allint_DMs,
#   dir_in =file.path(getwd(), "TMB"),
#   TMB_data_baseline_in=TMB_data_baseline,
#   calc_SEs=FALSE,cpp_name_in = "HOR_CHP_global")
# proc.time()-bt # 54.3 minutes

AIC_DF_tmp <- get_summ_tab_from_mod_ls(HOR_CHP_TMB_mod_fits_unscl_ls)


# TMB:::getUserDLL()
# dyn.unload(dynlib("HOR_CHP_global"))


# XX_pred_mat

# HOR_CHP_d2_AIC_unscl <- get_summ_tab_from_mod_ls(HOR_CHP_d2_mod_ls_unscl)

# if(!file.exists(file.path(out_dir, "HOR_CHP_TMB_mod_fits_unscl_ls.rds"))){
#   saveRDS(HOR_CHP_TMB_mod_fits_unscl_ls,file.path(out_dir, "HOR_CHP_TMB_mod_fits_unscl_ls.rds"))}

# TMB_data_baseline <- HOR_CHP_comp_ls$TMB_data_baseline
# XX_in_w_int_WYT <- HOR_CHP_comp_ls$XX_in_w_int_WYT


# assignig parameter names associated w each model parameter names
par_init_v_ls <- lapply(1:length(HOR_CHP_TMB_mod_fits_unscl_ls),function(ii) {
  HOR_CHP_TMB_mod_fits_unscl_ls[[ii]]$OPT$par})

par_mat1=do.call(rbind,strsplit(S_conf_ls$AIC_DF_d2$"dm",split=""))=="1"
# par_mat1=do.call(rbind,strsplit(allint_DMs,split=""))=="1"
colnames(par_mat1)=colnames(XX_in_w_int_WYT)
par_nm_ind_ls <- apply(par_mat1,1,which,simplify = F)

HOR_CHP_TMB_mod_fits_unscl_ls[[1]]$OPT$par
names(HOR_CHP_TMB_mod_fits_unscl_ls[[1]]$OPT$par)

all_mod_parm_coef_ls <- lapply(1:length(HOR_CHP_TMB_mod_fits_unscl_ls),function(ii){
  aa=HOR_CHP_TMB_mod_fits_unscl_ls[[ii]]$OPT$par
  S_coef=aa[names(aa)=="S_pars"]
  other_pars=aa[names(aa)!="S_pars"]
  coeff_param=c(S_coef,other_pars)
  par_nm_tmp=c(names(par_nm_ind_ls[[ii]]),names(aa)[names(aa)!="S_pars"])
  c=data.frame(dmID=ii,
               dm=S_conf_ls$AIC_DF_d2$"dm"[ii],
               # par_nm=names(par_nm_ind_ls[[ii]]),
               par_nm=par_nm_tmp,
               S_parID=c(par_nm_ind_ls[[ii]],rep(NA,length(other_pars))),
               coeff_param)
  rownames(c)=NULL
  return(c)
})

# HOR_CHP_d2_mod_ls_unscl$


ii=1
S_coef_ls <- lapply(1:length(HOR_CHP_d2_mod_ls_unscl),function(ii){
  aa=HOR_CHP_d2_mod_ls_unscl[[ii]]$"est_tab"
  # aa=HOR_CHP_d2_mod_ls_unscl[,ii]$"est_tab"
  S_coef=aa$"Estimate_SE"[aa$"Parameter"=="S_pars"]
  c=data.frame(dmID=ii,dm=S_conf_ls$AIC_DF_d2$"dm"[ii],par_nm=names(par_nm_ind_ls[[ii]]),S_parID=par_nm_ind_ls[[ii]],S_coef)
  rownames(c)=NULL
  return(c)})

all_mod_parm_coefDF <- do.call(rbind,all_mod_parm_coef_ls)  # all parameters (S, detection, random-effect SD), as in the saved object


HOR_CHP_TMB_all_mods_unscl <- list(
  "HOR_CHP_comp_ls"=HOR_CHP_TMB_mod_fits_unscl_ls,
  "AIC_DF_d2_full"=AIC_DF_tmp,
  "pt_estsDF"=all_mod_parm_coefDF,
  #elevating some objects for convenience
  "par_nm_ind_ls"=par_nm_ind_ls,
  # "allint_DMs"=HOR_CHP_comp_ls$"allint_DMs",
  "TMB_data_baseline"=HOR_CHP_comp_ls_unscl$"TMB_data_baseline",
  "XX_in_w_int_WYT"=HOR_CHP_TMB_mod_fits_unscl_ls$"XX_in_w_int_WYT",
  "scaled_vars"=attributes(HOR_CHP_TMB_mod_fits_unscl_ls$XX_in)$"scaled_vars",
  "center"=attributes(HOR_CHP_TMB_mod_fits_unscl_ls$XX_in)$"center",
  "scale"=attributes(HOR_CHP_TMB_mod_fits_unscl_ls$XX_in)$"scale")

# HOR_CHP_TMB_all_mods


# always overwritten (the guard that was here left a stale 16-model file after the HOR_CHP expansion)
saveRDS(HOR_CHP_TMB_all_mods_unscl,file.path(out_dir, "HOR_CHP_TMB_all_mods_unscl.rds"))



# HOR_CHP_d2_unscl_supp_ls<- list(
#   "S_conf_d2_DMs"=S_conf_ls$AIC_DF_d2$"dm",
#   "TMB_data_baseline"=HOR_CHP_comp_ls_unscl$"TMB_data_baseline",
#   "HOR_CHP_d2_mod_ls_unscl"=HOR_CHP_d2_mod_ls_unscl,
#   "HOR_CHP_d2_AIC_unscl"=HOR_CHP_d2_AIC_unscl)
# 
# HOR_CHP_comp_ls_unscl$TMB_data_baseline$xx
# 
# saveRDS(HOR_CHP_d2_unscl_supp_ls,file.path(out_dir, "HOR_CHP_d2_unscl_supp_ls.rds"))


# str(HOR_CHP_d2_mod_ls_unscl$`111111111110100001011`)
# HOR_CHP_d2_mod_ls_unscl$`111111111110100001011`$AIC


