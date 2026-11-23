# model_fitting/TCJ_CHP_conf_comp.R  (copy of Preliminaries/TCJ_CHP_conf_comp.R, renamed for the model_fitting pipeline)
# Setting (define before sourcing to override): folder results are written to / read from.
if (!exists("out_dir")) out_dir <- "output"
source("dAIC_threshold.R")

library(TMBhelp)
# devtools::load_all("../CVhelp")
TCJ_CHP_TMB_all_mods <- readRDS(file.path(out_dir, "TCJ_CHP_TMB_all_mods.rds"))

# loading functions that def work
# comp function is defined in the data_prep script (run_data_prep <- FALSE loads the function only)
run_data_prep <- FALSE
source("model_fitting/TCJ_CHP_surv_mods_data_prep.R")

# rlang::env_get_list()
# rlang::env_get(nm = ".dll")
# TMB:::getUserDLL()
# dyn.load(dyn.unload("TMB/TCJ_CHP_global"))
# 20 to 30
# TCJ_CHP_TMB_all_mods$AIC_DF_full[which(TCJ_CHP_TMB_all_mods$AIC_DF_full$dm=="111111111110000000000"),]
AIC_DF_full <- TCJ_CHP_TMB_all_mods$"AIC_DF_full"

AIC_DF_d2 <- TCJ_CHP_TMB_all_mods$"AIC_DF_full" |> dplyr::filter(dAIC<dAIC_threshold)
AIC_DF_d2 <- AIC_DF_d2 |> dplyr::mutate(candmodID=match(dm,TCJ_CHP_TMB_all_mods$allint_DMs))
# consider get_aic_wts()
AIC_DF_full$AICwt <- exp(-0.5*AIC_DF_full$dAIC)/sum(exp(-0.5*AIC_DF_full$dAIC))

# devtools::load_all("../TMBhelp")
TCJ_CHP_mod_fits_d2_ls <- TMBhelp::fit_nested_tmb_mods(
  S_design_matrix_in = AIC_DF_d2$dm, 
  dir_in =file.path(getwd(), "TMB"),
  TMB_data_baseline_in=TCJ_CHP_TMB_all_mods$"TMB_data_baseline",
  calc_SEs=TRUE,
  # TCJ_CHP fits deliberately use the HOR_CHP template: the two survival models share
  # one likelihood structure and no TCJ_CHP_global.cpp exists (see README.md, 'TMB templates').
  cpp_name_in = "HOR_CHP_global"
  )


TCJ_CHP_mod_fits_d2_ls$`1111110010101`$COV_mat
SD_obj <- TCJ_CHP_mod_fits_d2_ls$`1111110010101`$COV_mat
cov_matrix <- solve(SD_obj$jointPrecision)


TCJ_CHP_mod_fits_d2_ls[[1]]$est_tab$Estimate
S_coef_ls_confset <- lapply(1:length(TCJ_CHP_mod_fits_d2_ls),function(ii){
  aa=TCJ_CHP_mod_fits_d2_ls[[ii]]$"est_tab"
  S_coef=aa$"Estimate_SE"[aa$"Parameter"=="S_pars"]
  est_tmp=aa$"Estimate"[aa$"Parameter"=="S_pars"]
  se_tmp=aa$"Std..Error"[aa$"Parameter"=="S_pars"]
  c=data.frame(AICrank=ii,
               candmodID=AIC_DF_d2$candmodID[ii],
               dm=AIC_DF_d2$dm[ii],
               par_nm=names(TCJ_CHP_TMB_all_mods$"par_nm_ind_ls"[[AIC_DF_d2$candmodID[ii]]]),
               S_parID=TCJ_CHP_TMB_all_mods$"par_nm_ind_ls"[[AIC_DF_d2$candmodID[ii]]],
               S_coef,
               estimate=est_tmp,
               SE=se_tmp)
  rownames(c)=NULL
  return(c)})

S_coef_confset_DF <- do.call(rbind,S_coef_ls_confset)
S_conf_ls <- list(
  "AIC_DF_d2"=AIC_DF_d2,
  "S_coef_confset_DF"=S_coef_confset_DF)

saveRDS(S_conf_ls,file.path(out_dir, "TCJ_CHP_S_conf_ls.rds"))

head(TCJ_CHP_TMB_all_mods$"TMB_data_baseline")

# TCJ_CHP_mod_fits_d2_ls <- fit_nested_tmb_mods(
#   S_design_matrix_in = AIC_DF_d2$dm, 
#   dir_in =file.path(getwd(), "TMB"),
#   TMB_data_baseline_in=TCJ_CHP_TMB_all_mods$"TMB_data_baseline",
#   calc_SEs=TRUE,cpp_name_in = "TCJ_CHP_global")

stopifnot("run from the repository root (TMB_SD_sthd)" = file.exists("TMB_SD_sthd.Rproj"))

# setwd("..")
library(TMBhelp)
TCJ_CHP_mod_fits_d2_alt_sets_ls <- TMBhelp::fit_nested_tmb_mods(
  S_design_matrix_in = AIC_DF_d2$dm[1:2], 
  dir_in =file.path(getwd(), "TMB"),
  TMB_data_baseline_in=TCJ_CHP_TMB_all_mods$"TMB_data_baseline",
  calc_SEs=TRUE,
  cpp_name_in = "HOR_CHP_global"#,
  )

load("data/TCJ_CHP_mod_dat_ls.RData")


TCJ_CHP_comp_ls_unscl <- TCJ_CHP_comp(z_scale_vars = F,lvec.tcj_ls_in = lvec.tcj_ls,sdat.det.in=sdat.det.common2,x.df = x.df)

TCJ_CHP_comp_ls_unscl$XX_in_w_int_WYT
# tmp_scl <-attributes(scale(TCJ_CHP_comp_ls_unscl$XX_in_w_int_WYT[,grep(names(TCJ_CHP_comp_ls_unscl$XX_in_w_int_WYT),pattern="_x_",invert = T)]))
tmp_scl <-attributes(scale(TCJ_CHP_comp_ls_unscl$XX_in_w_int_WYT))


TCJ_CHP_comp_ls_scl <- TCJ_CHP_comp(z_scale_vars = T,lvec.tcj_ls_in = lvec.tcj_ls,sdat.det.in=sdat.det.common2,x.df = x.df)

str(TCJ_CHP_comp_ls_scl)
names(TCJ_CHP_comp_ls_scl)

# TCJ_CHP_comp_ls_scl

attributes(TCJ_CHP_comp_ls_scl$XX_in)




# attributes(TCJ_CHP_comp_ls_unscl$XX_in)
# attributes(TCJ_CHP_comp_ls_unscl$XX_in_w_int_WYT)

head(TCJ_CHP_comp_ls_unscl$XX_in_w_int_WYT)
head(TCJ_CHP_TMB_all_mods$XX_in_w_int_WYT)

# back-transformation of values
head(TCJ_CHP_TMB_all_mods$XX_in_w_int_WYT$flength*TCJ_CHP_TMB_all_mods$scale["flength"] + TCJ_CHP_TMB_all_mods$center["flength"])

head(TCJ_CHP_comp_ls_unscl$XX_in_w_int_WYT)
head(TCJ_CHP_comp_ls_unscl$XX_in)


################################################## #
# top-models with unscaled versions of the model
################################################## #

get_summ_tab_from_mod_ls(TCJ_CHP_mod_fits_d2_ls)

get_SEcoef_tab <- function(mod_sub_ls_in,mod_full_ls_in){
  out_DF=lapply(1:length(mod_sub_ls_in),function(ii){
  aa=mod_sub_ls_in[[ii]]$"est_tab"
  S_coef=aa$"Estimate_SE"[aa$"Parameter"=="S_pars"]
  est_tmp=aa$"Estimate"[aa$"Parameter"=="S_pars"]
  se_tmp=aa$"Std..Error"[aa$"Parameter"=="S_pars"]
  c=data.frame(AICrank=ii,
               candmodID=AIC_DF_d2$candmodID[ii],
               dm=AIC_DF_d2$dm[ii],
               par_nm=names(mod_full_ls_in$"par_nm_ind_ls"[[AIC_DF_d2$candmodID[ii]]]),
               S_parID=mod_full_ls_in$"par_nm_ind_ls"[[AIC_DF_d2$candmodID[ii]]],
               S_coef,
               estimate=est_tmp,
               SE=se_tmp)
  rownames(c)=NULL
  return(c)})
  return(out_DF)
}

S_coef_confset_ls <- get_SEcoef_tab(mod_full_ls_in = TCJ_CHP_TMB_all_mods,mod_sub_ls_in = TCJ_CHP_mod_fits_d2_ls)
S_coef_confset_DF <- do.call(rbind,S_coef_confset_ls)

# subset(S_coef_confset_DF,par_nm )
sdRELgrp_DF<- do.call(rbind,lapply(1:length(TCJ_CHP_mod_fits_d2_ls),function(x) subset(TCJ_CHP_mod_fits_d2_ls[[x]]$est_tab,Parameter=="logSD_RELGRP")))
rownames(sdRELgrp_DF) <- NULL

AIC_DF_d2$AICwt <- exp(-0.5*AIC_DF_d2$dAIC)/sum(exp(-0.5*AIC_DF_d2$dAIC))

S_coef_confset_DF$par_nm 
S_coef_confset_DF$AICwt <- AIC_DF_d2$AICwt[match(S_coef_confset_DF$dm,AIC_DF_d2$dm)]
library(ggplot2)
ggplot(data=S_coef_confset_DF,
         aes(x=factor(S_parID),
             y=estimate,ymin=estimate-SE*1.96,ymax=estimate+SE*1.96,
             size=AICwt*2.5,
             fill=factor(dm)
             )) + 
    geom_hline(yintercept = 0,linetype="solid",colour = "gray40") +
  geom_linerange(width=0,position = position_dodge(width=1),size=0.5) +
  # geom_errorbar(width=0,position = position_dodge(width=1)) + 
    geom_point(shape=21,position=position_dodge(width=1)) + 
  scale_x_discrete(limits=rev,labels=rev(unique(S_coef_confset_DF$par_nm)))  + coord_flip() 



# +
# tmp_in_ls <- list(TCJ_CHP_comp_ls$"TMB_data_baseline",
#                   TCJ_CHP_comp_ls_unscl$"TMB_data_baseline")
# # tmp_in_ls[[1]]$
# 
# tmp_out_ls <- list()
# for(ii in 1:2){
#   tmp_out_ls[[ii]] <- fit_nested_tmb_mods(
#     S_design_matrix_in = "111110000000000000000", 
#     dir_in =file.path(getwd(), "TMB"),
#     TMB_data_baseline_in=tmp_in_ls[[ii]],
#     # TMB_data_baseline_in=TCJ_CHP_comp_ls_unscl$"TMB_data_baseline",
#     calc_SEs=TRUE,cpp_name_in = "TCJ_CHP_global")}
# 
# tmp_out_coef_ls <- lapply(1:length(tmp_out_ls),function(ii){
#   data.frame(data_type=c("scl","unscl")[ii],tmp_out_ls[[ii]][[1]]$est_tab)})
# 
# tmp_outDF <- do.call(rbind,tmp_out_coef_ls)
# 
# tmp_outDF
# 
# tmp_out_coef_ls
# 
# 
# # tmp_out_ls[[1]]
# # get_SEcoef_tab(mod_sub_ls_in = tmp_out_ls)
# 
# 
# 
# 
# 
# TCJ_CHP_TMB_all_mods$XX_in_w_int_WYT  <- TCJ_CHP_TMB_all_mods$TCJ_CHP_comp_ls$"XX_in_w_int_WYT"
# 
# TCJ_CHP_TMB_all_mods$TCJ_CHP_comp_ls$XX_in
# TCJ_CHP_TMB_all_mods$"XX_in_w_int_WYT"
# # TCJ_CHP_TMB_all_mods$"XX_in_w_int_WYT"
# 
# TCJ_CHP_TMB_all_mods$TCJ_CHP_comp_ls$XX_in
# 
# subset(S_conf_ls$S_coef_confset_DF,AICrank==1)
# 
# TCJ_CHP_TMB_all_mods$center
# 
# ggplot(data=S_coef_confset_DF,
#        aes(x=factor(S_parID),
#            y=estimate,ymin=estimate-SE*1.96,ymax=estimate+SE*1.96,
#            fill=factor(AICrank))) + 
#   geom_errorbar(width=0,position = position_dodge(width=1)) + 
#   geom_point(size=2,shape=21,position=position_dodge(width=1)) + 
#   geom_hline(yintercept = 0,linetype="dotted")
# 
# S_coef_confset_DF[S_coef_confset_DF$S_parID %in% c(18,20,21),] |> arrange(par_nm)
# 
# table(S_coef_confset_DF$S_parID,S_coef_confset_DF$par_nm)
# 
# # S_coef_DF_d2 <- S_coef_DF |> filter(dm %in% AIC_DF_d2$dm)
# 
# coef_comp_d2 <- data.frame(table(S_coef_DF_d2$par_nm)) |> 
#   mutate(param=factor(Var1,levels=colnames(par_mat1))) |> 
#   arrange(param) |> select(-Var1) |> relocate(param) |>
#   rename(n_mods=Freq)
# coef_comp_d2
# 
# 
# 
# attributes(TCJ_CHP_TMB_all_mods.rds$TCJ_CHP_comp_ls$XX_in)$"scaled_vars"
# attributes(TCJ_CHP_TMB_all_mods.rds$TCJ_CHP_comp_ls$XX_in)$"center"
# attributes(TCJ_CHP_TMB_all_mods.rds$TCJ_CHP_comp_ls$XX_in)$"scale"
# 
# 
# 
# AIC_DF_tmp <- get_summ_tab_from_mod_ls(TCJ_CHP_TMB_mod_fits_ls)
# 
# # TCJ_CHP_mod_fits_1_1024_ls <- readRDS(file.path(out_dir, "TCJ_CHP_mod_fits_1_1024_ls.rds"))
# 
# # assignig parameter names associated w each model parameter names
# par_init_v_ls <- lapply(1:ncol(TCJ_CHP_TMB_mod_fits_ls),function(ii) {
#   TCJ_CHP_TMB_mod_fits_ls[,ii]$OPT$par})
# 
# XX_in_w_int_WYT <- TCJ_CHP_TMB_mod_fits_ls$XX_in_w_int_WYT
# par_mat1=do.call(rbind,strsplit(allint_DMs,split=""))=="1"
# colnames(par_mat1)=colnames(XX_in_w_int_WYT)
# par_nm_ind_ls <- apply(par_mat1,1,which,simplify = F)
# 
# 
# 
# 
# 
# get_DM_col_inds(allint_DMs)
# 
# # assignig parameter names associated w each model parameter names
# par_init_v_ls <- lapply(1:ncol(TCJ_CHP_TMB_mod_fits_ls),function(ii) {
#   TCJ_CHP_TMB_mod_fits_ls[,ii]$OPT$par})
# XX_in_w_int_WYT <- TCJ_CHP_TMB_mod_fits_ls$XX_in_w_int_WYT
# par_mat1=do.call(rbind,strsplit(allint_DMs,split=""))=="1"
# colnames(par_mat1)=colnames(XX_in_w_int_WYT)
# par_nm_ind_ls <- apply(par_mat1,1,which,simplify = F)
# 
# 
# AIC_DF_tmp <- get_summ_tab_from_mod_ls(TCJ_CHP_TMB_mod_fits_ls)
# 
# # old column method
# AIC_DF_tmp <- data.frame(
#   dm=allint_DMs,
#   k=sapply(1:ncol(TCJ_CHP_mod_fits_1_1024_ls),function(ii){TCJ_CHP_mod_fits_1_1024_ls[,ii]$"k"}),
#   AIC=sapply(1:ncol(TCJ_CHP_mod_fits_1_1024_ls),function(ii){TCJ_CHP_mod_fits_1_1024_ls[,ii]$"AIC"}))
# AIC_DF_full <- AIC_DF_tmp |> mutate(dAIC=AIC-min(AIC)) |> arrange(AIC)
# 
# S_coef_ls <- lapply(1:ncol(TCJ_CHP_mod_fits_1_1024_ls),function(ii){
#   aa=TCJ_CHP_mod_fits_1_1024_ls[,ii]$"est_tab"
#   S_coef=aa$"Estimate_SE"[aa$"Parameter"=="S_pars"]
#   c=data.frame(dmID=ii,dm=allint_DMs[ii],par_nm=names(par_nm_ind_ls[[ii]]),S_parID=par_nm_ind_ls[[ii]],S_coef)
#   rownames(c)=NULL
#   return(c)})
# 
# S_coef_DF <- do.call(rbind,S_coef_ls)
# 
# # looking at models where delta AIC <=2
# AIC_DF_d2 <- AIC_DF_full |> filter(dAIC<=2)
# AIC_DF_d2 <- AIC_DF_d2 |> mutate(candmodID=match(dm,allint_DMs))
# 
# S_coef_DF_d2 <- S_coef_DF |> filter(dm %in% AIC_DF_d2$dm)
# 
# coef_comp_d2 <- data.frame(table(S_coef_DF_d2$par_nm)) |> 
#   mutate(param=factor(Var1,levels=colnames(par_mat1))) |> 
#   arrange(param) |> select(-Var1) |> relocate(param) |>
#   rename(n_mods=Freq)
# coef_comp_d2
# 
# ################################################## #
# # re-run models and compute SEs for coefficients 
# ################################################## #
# 
# # bt <- proc.time()
# # TCJ_CHP_mod_fits_d2_ls=TMBhelp::mult_tmb_fits_ls_out(
# #   S_design_matrix_in = AIC_DF_d2$dm, # substitute the top 13 modls
# #   TMB_data_baseline_in=TMB_data_baseline,
# #   dir_in = file.path(getwd(), "TMB"),
# #   calc_SEs=TRUE,cpp_name_in = "TCJ_CHP_global")
# # proc.time()-bt
# 
# devtools::load_all("../TMBhelp")
# TCJ_CHP_mod_fits_d2_ls <- fit_nested_tmb_mods(
#   # S_design_matrix_in = allint_DMs[c(1:length(allint_DMs))],
#   # S_design_matrix_in = allint_DMs[c(1:3)],
#   S_design_matrix_in = AIC_DF_d2$dm, # substitute the top 13 modls
#   dir_in =file.path(getwd(), "TMB"),
#   TMB_data_baseline_in=TMB_data_baseline,
#   calc_SEs=TRUE,cpp_name_in = "TCJ_CHP_global")
# 
# TCJ_CHP_mod_fits_d2_ls
# 
# 
# 
# 
# # fit_nested_tmb_mods
# # list of estimates and SEs
# # sapply(1:ncol(TCJ_CHP_mod_fits_d2_ls),function(ii){TCJ_CHP_mod_fits_d2_ls[,ii]$"est_tab"$"Estimate_SE"})
# # TCJ_CHP_mod_fits_d2_ls
# 
# S_coef_ls_confset <- lapply(1:ncol(TCJ_CHP_mod_fits_d2_ls),function(ii){
#   aa=TCJ_CHP_mod_fits_d2_ls[,ii]$"est_tab"
#   S_coef=aa$"Estimate_SE"[aa$"Parameter"=="S_pars"]
#   c=data.frame(AICrank=ii,
#                candmodID=AIC_DF_d2$candmodID[ii],
#                dm=AIC_DF_d2$dm[ii],
#                par_nm=names(par_nm_ind_ls[[AIC_DF_d2$candmodID[ii]]]),
#                S_parID=par_nm_ind_ls[[AIC_DF_d2$candmodID[ii]]],
#                S_coef)
#   rownames(c)=NULL
#   return(c)})
# 
# 
# S_coef_ls_confset <- lapply(1:length(TCJ_CHP_mod_fits_d2_ls),function(ii){
#   aa=TCJ_CHP_mod_fits_d2_ls[[ii]]$"est_tab"
#   S_coef=aa$"Estimate_SE"[aa$"Parameter"=="S_pars"]
#   c=data.frame(AICrank=ii,
#                candmodID=AIC_DF_d2$candmodID[ii],
#                dm=AIC_DF_d2$dm[ii],
#                par_nm=names(par_nm_ind_ls[[AIC_DF_d2$candmodID[ii]]]),
#                S_parID=par_nm_ind_ls[[AIC_DF_d2$candmodID[ii]]],
#                S_coef)
#   rownames(c)=NULL
#   return(c)})
# 
# S_coef_confsetDF <-  do.call(rbind,S_coef_ls_confset)
# 
# 
# ################################################## #
# # Release-level random effects for models
# ################################################## #
# 
# TCJ_CHP_RE_d2_ls <- lapply(1:ncol(TCJ_CHP_mod_fits_d2_ls),function(ii){
#   aa=data.frame(ii,candmodID=AIC_DF_d2$candmodID[ii],dm=AIC_DF_d2$dm[ii],TCJ_CHP_mod_fits_d2_ls[,ii]$"RE_tab")
#   bb <- data.frame(param=rownames(aa),aa,rgID=c(1:19,NA))
#   bb$param[grep(x=bb$param,pattern="RELGRP_err")] <- "RELGRP_err"
#   bb
# })
# 
# TCJ_CHP_RE_d2_ls <- lapply(1:length(TCJ_CHP_mod_fits_d2_ls),function(ii){
#   aa=data.frame(ii,candmodID=AIC_DF_d2$candmodID[ii],
#                 dm=AIC_DF_d2$dm[ii],
#                 TCJ_CHP_mod_fits_d2_ls[[ii]]$"RE_tab")
#   bb <- data.frame(param=rownames(aa),aa)
#   bb$param[grep(x=bb$param,pattern="RELGRP_err")] <- "RELGRP_err"
#   bb
# })
# 
# 
# RE_tmp <-  do.call(rbind,TCJ_CHP_RE_d2_ls)
# RE_sigma_DF_d2 <- RE_tmp |> filter(param=="Sigma_RELGRP")
# RE_comb_DF_d2 <- RE_tmp |> filter(param=="RELGRP_err")
# rownames(RE_comb_DF_d2) =NULL
# rownames(RE_sigma_DF_d2) =NULL
# 
# RE_sigma_DF_d2
# RE_comb_DF_d2
# 
# # RE_sigma_DF_d2
# ggplot2::ggplot(data=RE_sigma_DF_d2,
#                 ggplot2::aes(y=Estimate,x=factor(ii),color=factor(ii),ymin=Estimate-Std..Error,ymax=Estimate+Std..Error)) + 
#   ggplot2::geom_errorbar(position = ggplot2::position_dodge(width=0.5),width=0) + 
#   ggplot2::geom_point(position = ggplot2::position_dodge(width=0.5))
# 
# ggplot2::ggplot(data=RE_comb_DF_d2,
#                 ggplot2::aes(y=Estimate,x=factor(rgID),color=factor(ii),ymin=Estimate-Std..Error,ymax=Estimate+Std..Error)) + 
#   ggplot2::geom_errorbar(position = ggplot2::position_dodge(width=1),width=0)  + 
#   ggplot2::geom_point(position = ggplot2::position_dodge(width=1)) 
# 
# # model="S_HOR-CHP"
# RE_comb_DF_d2
# 
# TCJ_CHP_mod_fits_d2_ls
# S_coef_confsetDF
# AIC_DF_d2
# 
# 
# 
# 
# # head(TCJ_CHP_mod_fits_d2_ls)
# # TCJ_CHP_mod_fits_d2_ls[,1]
# 
# 
# TCJ_CHP_mod_fits_1_1024_ls[,1]$est_tab
# sapply(1:ncol(TCJ_CHP_mod_fits_1_1024_ls),function(ii){TCJ_CHP_mod_fits_1_1024_ls[,ii]$"est_tab"$"Estimate_SE"})
# 
# 
# 
# devtools::load_all("../TMBhelp")
# TMBhelp::mult_tmb_fits_ls_out(
#   S_design_matrix_in = AIC_DF_d2$dm[1:2], # substitute the top 13 modls
#   TMB_data_baseline_in=TMB_data_baseline,
#   dir_in = file.path(getwd(), "TMB"),
#   calc_SEs=TRUE,cpp_name_in = "TCJ_CHP_global",simplify=F)
# 
# # dim(TCJ_CHP_mod_fits_1_1024_ls)
# sapply(1:ncol(TCJ_CHP_mod_fits_1_1024_ls),function(ii){TCJ_CHP_mod_fits_1_1024_ls[,ii]$"AIC"})
# sapply(1:length(TCJ_CHP_mod_fits_1_1024_ls),function(ii){TCJ_CHP_mod_fits_1_1024_ls[,ii]})
# 
# # 2800 seconds
# # saveRDS(TCJ_CHP_mod_fits_1_1024,"../output/TCJ_CHP_mod_fits_1_1024.rds")
# 
# # AIC table summarizing all model fits
# TCJ_CHP_mod_fits_1_1024 <- readRDS(file.path(out_dir, "TCJ_CHP_mod_fits_1_1024.rds"))
# str(TCJ_CHP_mod_fits_1_1024)
# 
# dir("output")
# 
# 
# 
# # TCJ_CHP_mod_fits_1_1024_ls[,1]$est_tab
# # code for extracting initialization format of inputs
# # par_init_ls <- lapply(1:ncol(TCJ_CHP_mod_fits_1_1024_ls),function(ii) {
# #   par_v=TCJ_CHP_mod_fits_1_1024_ls[,ii]$OPT$par
# #   list(
# #     "S_pars"=par_v[names(par_v)=="S_pars"],
# #     "P_mat"=matrix(par_v[names(par_v)=="P_pars"],
# #                    nrow=2))})
# 
# 
