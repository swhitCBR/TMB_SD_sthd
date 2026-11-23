# model_fitting/HOR_CHP_conf_comp.R
#
# Confidence set for the HOR_CHP survival models. Refits the candidate models with
# dAIC < dAIC_threshold (using the HOR_CHP_global template) and writes the design
# matrices, estimates and standard errors of the confidence set to
# HOR_CHP_S_conf_ls.rds, which HOR_CHP_d2_AIC_unscl_mod_ls.R and the memo cross-tab read.
#
# Replaces the S_conf_ls part of Preliminaries/HOR_CHP_conf_comp_old_bad.R, which was
# the only writer of this file.
#
# Input:  out_dir/HOR_CHP_TMB_all_mods.rds   (written by HOR_CHP_cand_mods_run.R)
# Output: out_dir/HOR_CHP_S_conf_ls.rds
# Setting (define before sourcing to override): folder results are written to / read from.
# Run from the TMB_SD_sthd root.
if (!exists("out_dir")) out_dir <- "output"

source("dAIC_threshold.R")
library(TMBhelp)

HOR_CHP_TMB_all_mods <- readRDS(file.path(out_dir, "HOR_CHP_TMB_all_mods.rds"))

AIC_DF_d2 <- HOR_CHP_TMB_all_mods$"AIC_DF_full" |> dplyr::filter(dAIC<dAIC_threshold)
AIC_DF_d2 <- AIC_DF_d2 |> dplyr::mutate(candmodID=match(dm,HOR_CHP_TMB_all_mods$allint_DMs))
# Akaike weights renormalised over the confidence set (the stored file has this column)
AIC_DF_d2$AICwt <- exp(-0.5*AIC_DF_d2$dAIC)/sum(exp(-0.5*AIC_DF_d2$dAIC))

HOR_CHP_mod_fits_d2_ls <- TMBhelp::fit_nested_tmb_mods(
  S_design_matrix_in = AIC_DF_d2$dm,
  dir_in = file.path(getwd(), "TMB"),
  TMB_data_baseline_in=HOR_CHP_TMB_all_mods$"TMB_data_baseline",
  calc_SEs=TRUE,
  cpp_name_in = "HOR_CHP_global")

S_coef_ls_confset <- lapply(1:length(HOR_CHP_mod_fits_d2_ls),function(ii){
  aa=HOR_CHP_mod_fits_d2_ls[[ii]]$"est_tab"
  S_coef=aa$"Estimate_SE"[aa$"Parameter"=="S_pars"]
  est_tmp=aa$"Estimate"[aa$"Parameter"=="S_pars"]
  se_tmp=aa$"Std..Error"[aa$"Parameter"=="S_pars"]
  c=data.frame(AICrank=ii,
               candmodID=AIC_DF_d2$candmodID[ii],
               dm=AIC_DF_d2$dm[ii],
               par_nm=names(HOR_CHP_TMB_all_mods$"par_nm_ind_ls"[[AIC_DF_d2$candmodID[ii]]]),
               S_parID=HOR_CHP_TMB_all_mods$"par_nm_ind_ls"[[AIC_DF_d2$candmodID[ii]]],
               S_coef,
               estimate=est_tmp,
               SE=se_tmp)
  rownames(c)=NULL
  return(c)})

S_coef_confset_DF <- do.call(rbind,S_coef_ls_confset)

S_conf_ls <- list(
  "AIC_DF_d2"=AIC_DF_d2,
  "S_coef_confset_DF"=S_coef_confset_DF)

saveRDS(S_conf_ls, file.path(out_dir, "HOR_CHP_S_conf_ls.rds"))
