# model_fitting/HOR_CHP_cand_mods_run.R  (copy of Preliminaries/HOR_CHP_all_mod_run_REMOTE_posit.R, renamed for the model_fitting pipeline)
# Setting (define before sourcing to override): folder results are written to / read from.
if (!exists("out_dir")) out_dir <- "output"
if (!exists("n_workers")) n_workers <- 8   # worker processes used by fit_nested_tmb_mods_par()

library(TMB)
library(dplyr)

# devtools::load_all("../TMBhelp")
# devtools::load_all("../CVhelp")

load("data/HOR_CHP_mod_dat_ls.RData")

# library(CVhelp)  # no longer needed: comp functions are in HOR_CHP_surv_mods_data_prep.R
library(TMBhelp)
# comp functions are defined in the data_prep script (run_data_prep <- FALSE loads the functions only)
run_data_prep <- FALSE
source("model_fitting/HOR_CHP_surv_mods_data_prep.R")

# compile data from HOR-CHP survival model
HOR_CHP_comp_ls <- HOR_CHP_comp()
HOR_CHP_comp_ls_unscl <- HOR_CHP_comp(z_scale_vars = F)

attributes(HOR_CHP_comp_ls$XX_in)
attributes(HOR_CHP_comp_ls_unscl$XX_in)

# extracting a few relevant vars
allint_DMs <- HOR_CHP_comp_ls$allint_DMs
TMB_data_baseline <- HOR_CHP_comp_ls$TMB_data_baseline
# XX_in_w_int_WYT <- HOR_CHP_TMB_mod_fits_ls$XX_in_w_int_WYT

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
