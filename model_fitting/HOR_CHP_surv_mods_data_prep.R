############################################################################### #
# HOR-CHP survival model: data preparation and design matrices
#
# Two parts:
#   1. Function definitions (always loaded): scale_data_cols() and the legacy
#      HOR_CHP_comp(). The pipeline no longer calls HOR_CHP_comp(): the design
#      matrix, TMB data list and candidate model strings (allint_DMs) are built by
#      HOR_CHP_precomp() and HOR_CHP_modelprep() in R/HOR_CHP_comp_ALT.R, with the
#      formula-based design in R/HOR_CHP_designmat.R and R/legacy_adapters.R
#      (see model_fitting/README.md and the reference/ folders).
#   2. Data preparation (runs only if run_data_prep is TRUE, the default):
#      reduces the original Buchanan et al. (2021) workspaces to the objects
#      HOR_CHP_precomp() reads.
#
# Current design (settings in config.R): the 11 always-in terms are (Intercept),
#   WYT_wet, WYT_drought, route.facB, barrier.facTRUE, flength, Tmsd.hor.7dadm,
#   log.VNS.hor.5, SWP.hor.5, CVP.hor.5 and Qomt.hor.1net; the 10 toggled terms are
#   barrier x {Temp, VNS, SWP, CVP, OMT} and route x {Temp, VNS, SWP, CVP, OMT}
#   (2^10 = 1,024 models). WYT_ follows each row's water.yr (CHP_wyt = "water.yr")
#   and R_x_Temp uses Tmsd.hor.7dadm (HOR_CHP_r_x_temp). The legacy design used
#   Tclc.hor.3 in R_x_Temp and a release-group lookup for WYT_drought.
#
# Part 2 inputs (searched via input_dirs in config.R): tmp_img4FULL.Rdata
#                 Multiyear v2a AllVar v3.RData
# Part 2 outputs: data/HOR_CHP_mod_dat_ls.RData (read by HOR_CHP_precomp()),
#                 data/HOR_CHP_mod_dat_sep_ls.rds
# Run from the TMB_SD_sthd root.
#
# ----------------------------------------------------------------------------- #
# LEGACY: VARIABLE SUBSTITUTIONS AND COLUMNS USED BY HOR_CHP_comp()
# (kept for reference; the description below is of the legacy function)
# ----------------------------------------------------------------------------- #
# Rows:       x.df rows with a non-missing det.pred.qomt.rms (x.df_rem1)
# Kept as columns of x.df_SUB (HOR_CHP_VARS):
#             flength, log.QOUT, Tmsd.hor.7dadm, log.VNS.hor.5, SWP.hor.5,
#             CVP.hor.5, Qomt.hor.1net, Tclc.hor.3, VNS.rel, VNS.hor
# Design matrix columns (from Xc.3v8): (Intercept), year.fac2011 (built as the
#             complement of the 2012-2016 year indicators), year.fac2012-2016,
#             route.facB, barrier.facTRUE, route.facB:barrier.facTRUE
# Substitutions:
#   - year effects are dropped from the survival design matrix (columns 2-7) and
#     replaced by WYT_wet and WYT_drought indicators, derived from each release
#     group's water-year type (so predictions can be made outside study years)
#   - "log(VNS.hor.5)" from Xc.3v8 is NOT used; the log.VNS.hor.5 column of
#     x.df is used instead
#   - main effects (MAIN_EFFECT_IDS): flength, Tmsd.hor.7dadm, log.VNS.hor.5,
#     SWP.hor.5, CVP.hor.5, Qomt.hor.1net
#   - Tclc.hor.3 is NOT a main effect: it enters only through R_x_Temp, so the
#     route x temperature interaction uses CLC temperature (Tclc.hor.3) while the
#     barrier x temperature interaction uses MSD temperature (Tmsd.hor.7dadm)
#   - log.QOUT is used only for detection (P_col_ids = columns 1 and 20)
#   - z_scale_vars = TRUE centres and scales flength, Tmsd.hor.7dadm,
#     log.VNS.hor.5, SWP.hor.5, CVP.hor.5, Qomt.hor.1net and Tclc.hor.3
#     (scale_data_cols()); centre/scale are stored as attributes of XX_in
# Interactions: B_x_{Temp,VNS,SWP,CVP,OMT} and R_x_{Temp,VNS,SWP,CVP,OMT}
# Candidate set (allint_DMs): the first 11 terms are always in; the 10
#             interactions are toggled, giving 2^10 = 1,024 models. Expanded
#             candidate sets are built in dev/ (mk_model_strings()).
############################################################################### #

library(dplyr)

if (!exists("run_data_prep")) run_data_prep <- TRUE

# ---- Part 1: functions ------------------------------------------------------ #
#
source("R/scale_data_cols.R")

#' Title
#'
#' @param RData_pth_in  path of model data
#' @param z_scale_vars logical, scale and center variables or not
#'
#' @returns named list with vital 
#' 
#'
HOR_CHP_comp <- function(RData_pth_in="data/HOR_CHP_mod_dat_ls.RData",z_scale_vars=TRUE){
  load(RData_pth_in)
  
  # ommitting NAs
  x.df_rem1 <- subset(x.df,#tag!="1232531" &
                      !is.na(det.pred.qomt.rms))
  
  # for matching up inputs
  tmpDF <- x.df_rem1 %>% 
    select(tag,year,rel.grp,HOR.det, ,flength,barrier.fac) %>% 
    mutate(rel_grp=paste(year,rel.grp),
           model="survival",
           foc_area="HOR-CHP")
  
  WYT_ref_tab <- table(paste(x.df_rem1$year,x.df_rem1$rel.grp),x.df_rem1$water.yr)
  
  BASE_COLS <- c("form","rel.grp","year","water.yr","route.fac","barrier")
  HOR_CHP_VARS  <- c(
    "flength",         # Fork length
    "log.QOUT",        # OUTFLOW
    "Tmsd.hor.7dadm",  # MSD temp
    "log.VNS.hor.5",   # Vernalis flow
    "SWP.hor.5",       # Exports
    "CVP.hor.5",       # Exports
    "Qomt.hor.1net",   # Flow
    "Tclc.hor.3",
    "VNS.rel",
    "VNS.hor"
  )
  
  Xc.3v8 <- cbind(Xc.3v8,"year.fac2011"=(apply(Xc.3v8[,c("year.fac2012","year.fac2013","year.fac2014","year.fac2015","year.fac2016")],1,sum)==0)+0)
  x.df_SUB <- cbind(Xc.3v8[,c("(Intercept)","year.fac2011",
                              "year.fac2012","year.fac2013","year.fac2014","year.fac2015","year.fac2016",
                              "route.facB","barrier.facTRUE","route.facB:barrier.facTRUE"#,
                              # "log(VNS.hor.5)"
  )],
  x.df_rem1[,c(BASE_COLS,HOR_CHP_VARS)])
  # excluded row
  rw_ind <- which(x.df$tag!="1232531")
  XX_in <- x.df_SUB
  
  
  repl_tab <- data.frame(
    pattrn=c("X.Intercept.","route.facB.barrier.facTRUE"),
    rplace=c("(Intercept)","route.facB:barrier.facTRUE"))
  
  scl_var_nms <- c("flength","Tmsd.hor.7dadm","log.VNS.hor.5","SWP.hor.5" ,"CVP.hor.5","Qomt.hor.1net","Tclc.hor.3")
  
  if(z_scale_vars){
    XX_in <- scale_data_cols(
      repl_tab_in=repl_tab,
      x.df_SUB_in=x.df_SUB)
  } else{
    XX_in <- x.df_SUB
    var_scl_info_tmp <- list("scaled_vars" = scl_var_nms,
                         "center" = rep(0,length(scl_var_nms)),
                         "scale" =  rep(1,length(scl_var_nms)))
    attributes(XX_in) <- c(attributes(XX_in), var_scl_info_tmp)
  }
  
  
  # XX_in <- XX_in_alt
  # XX_in <- x.df_SUB # uncomment and run here to get unscaled estimates
  # XX_in <- XX_in_alt

  # vector of  4 linear combinations used in only for loop in CPP template
  P_lc_n <- length(c("P1_f1","Lam_f1",
                     "P11_f2f3","P12_f2f3"))
  # row_inds for different likelihood forms
  f1 <- which(XX_in[,"form"]==1)
  f2 <- which(XX_in[,"form"]==2)
  f3 <- which(XX_in[,"form"]==3)
  
  BASIC_IND_COLS <- c(1:9)
  MAIN_EFFECT_IDS <- match(c("flength","Tmsd.hor.7dadm","log.VNS.hor.5","SWP.hor.5","CVP.hor.5","Qomt.hor.1net"),names(XX_in))
  
  # only in the intercept
  # INT_ONLY_IDS <- match(c("Tclc.hor.3"),names(XX_in))
  
  S_col_add_ids <- c(BASIC_IND_COLS,MAIN_EFFECT_IDS)
  Bar_col_id <- 9 # column with BAR
  Rte_col_id <- 8 # column with ROUTE
  
  # additive design matrix
  XX_in_add <- XX_in[,S_col_add_ids]
  
  # adding one-way interactions design matrix
  XX_in_w_int <- cbind(XX_in_add,
                       data.frame(
                         # Barrier interactions
                         B_x_Temp=XX_in_add[,Bar_col_id]*XX_in_add[,c("Tmsd.hor.7dadm")],
                         B_x_VNS=XX_in_add[,Bar_col_id]*XX_in_add[,c("log.VNS.hor.5")],
                         B_x_SWP=XX_in_add[,Bar_col_id]*XX_in_add[,c("SWP.hor.5")],
                         B_x_CVP=XX_in_add[,Bar_col_id]*XX_in_add[,c("CVP.hor.5")],
                         B_x_OMT=XX_in_add[,Bar_col_id]*XX_in_add[,c("Qomt.hor.1net")],
                         # Route interactions
                         R_x_Temp=XX_in[,Rte_col_id]*XX_in[,c("Tclc.hor.3")], # does not take column from main effect
                         R_x_VNS=XX_in_add[,Rte_col_id]*XX_in_add[,c("log.VNS.hor.5")],
                         R_x_SWP=XX_in_add[,Rte_col_id]*XX_in_add[,c("SWP.hor.5")],
                         R_x_CVP=XX_in_add[,Rte_col_id]*XX_in_add[,c("CVP.hor.5")],
                         R_x_OMT=XX_in_add[,Rte_col_id]*XX_in_add[,c("Qomt.hor.1net")]
                       ))
  
  # columns used to model survival (intercept and outflow)
  P_col_ids <- c(1,20)
  
  rel_grps_df_tmp <- data.frame(
    rel_grps=unique(paste(XX_in$year,XX_in$rel.grp)),
    t(sapply(unique(paste(XX_in$year,XX_in$rel.grp)),function(x) strsplit(x,split = " ")[[1]])))
  rel_grps_df_tmp[,2]=as.numeric(rel_grps_df_tmp[,2])
  names(rel_grps_df_tmp)[c(2,3)]=c("year","rel_nm")
  
  rel_grps_df_tmp <- rel_grps_df_tmp[order(rel_grps_df_tmp$year,rel_grps_df_tmp$rel_nm),]
  rel_grps_df_tmp$relgrpF <- factor(rel_grps_df_tmp$rel_grps,levels = rel_grps_df_tmp$rel_grps)
  
  tmp_tab <- which(WYT_ref_tab!=0,arr.ind = T)
  rel_grps_df <- data.frame(rel_grps_df_tmp,
                            WYT=c("wet","dry","drought")[tmp_tab[tmp_tab[,"row"],"col"]],
                            R_ind=as.numeric(rel_grps_df_tmp$relgrpF),
                            TMB_ind=as.numeric(rel_grps_df_tmp$relgrpF)-1)
  
  relgrps_ind_IN <- rel_grps_df$TMB_ind[match(paste(XX_in$year,XX_in$rel.grp),rel_grps_df$rel_grps)]
  
  WYT_wet=(relgrps_ind_IN %in% rel_grps_df$TMB_ind[rel_grps_df$WYT=="wet"])+0
  WYT_drought=(relgrps_ind_IN %in% rel_grps_df$TMB_ind[rel_grps_df$WYT=="drought"])+0
  
  XX_in_w_int_noYRS <- XX_in_w_int[-c(2:7)]
  XX_in_w_int_WYT <- cbind(XX_in_w_int_noYRS[,"(Intercept)"],
                           WYT_wet,WYT_drought,
                           XX_in_w_int_noYRS[,which(names(XX_in_w_int_noYRS)!="(Intercept)")])
  names(XX_in_w_int_WYT)[1] <- "(Intercept)"
  
  
  # XX_in
  
  
  TMB_data_baseline <- list(
    "n_relgrps"=19,
    "relgrps_ind"=relgrps_ind_IN,
    "XX"=as.matrix(XX_in[,P_col_ids]),
    "XX_s"=as.matrix(XX_in_w_int_WYT),
    "Dat"=sdat.det.common2,
    "col_Dat1"=match(paste(c("CHP.BBR"),c("11","10","01","00"),sep="."),colnames(sdat.det)) - 1,
    "col_Dat2"=match(paste(c("CHP.2x"),c("11","10","01","00"),sep="."),colnames(sdat.det)) - 1,
    "col_Dat3"=match(paste(c("CHP.3x"),c("111","110","101","100","011","010","001","000"),sep="."),colnames(sdat.det))-1,
    "ind_all"=1:nrow(sdat.det.common2)-1,
    "f1_ind"=f1-1,
    "f2_ind"=f2-1,
    "f3_ind"=f3-1,
    "lvec1"=lvec1.common2,
    "lvec2"=lvec2.common2,
    "col_S"=(1:ncol(XX_in_w_int_WYT)),
    "col_P11"=(1:2)-1,
    "col_P12"=(1:2)-1,
    "y1to7"=(1:7)-1,
    "y1to4"=(1:4)-1,
    "y1256"=c(1,2,5,6)-1,
    "y1357"=c(1,3,5,7)-1,
    "y5to7"=(5:7)-1,
    "y347"=c(3,4,7)-1,
    "y246"=c(2,4,6)-1,
    "y123"=(1:3)-1,
    "y12"=(1:2)-1,
    "y13"=c(1,3)-1,
    "P_lc_n"=P_lc_n)
  
  # initial values of detection probability
  P_mat_base <- matrix(0,nrow=length(TMB_data_baseline$col_P11),
                       ncol=TMB_data_baseline$P_lc_n)
  
  S_pars_input <- rep(0,ncol(TMB_data_baseline$XX_s))
  
  # specifying combinations of parameters to add or remove
  allint_DMs <- apply(cbind(rep(c("11111111111"),1024),
                            expand.grid(
                              c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1)))
                      ,1,paste0,collapse="")
  
  # extracting column ids based on vectors that represent design matrices
  # get_DM_col_inds(allint_DMs[1])
  # get_DM_col_inds(allint_DMs[1024])
  
  if(!z_scale_vars){attributes(XX_in) <- c(attributes(XX_in), var_scl_info_tmp)}
  
  
  out_ls=list("allint_DMs"=allint_DMs,
              "TMB_data_baseline"=TMB_data_baseline,
              "XX_in"=XX_in,
              "XX_in_w_int_WYT"=XX_in_w_int_WYT)
  
  return(out_ls)
}


# ---- Part 2: data preparation ----------------------------------------------- #

# local(): the loaded workspaces do not stay in the global environment
# (the *_cand_mods_run.R scripts save the global environment to out_dir)
if (run_data_prep) local({

  if (!exists("input_path")) source("model_fitting/config.R")
  load(input_path("tmp_img4FULL.Rdata"))
  load(input_path("Multiyear v2a AllVar v3.RData"))

  # objs_for_surv_HOR_CHP_mods <- c("x.df","Xc.3v8","sdat.det.common2","sdat.det","lvec1.common2","lvec2.common2")

  objs_for_surv_HOR_CHP_mods <- c("x.df","Xc.3v8","sdat.det.common2","sdat.det","lvec1.common2","lvec2.common2")
  
  save(list=objs_for_surv_HOR_CHP_mods,file="data/HOR_CHP_mod_dat_ls.RData")
  
  
    HOR_CHP_mod_dat_sep_ls <- list(
    "x.df"=x.df,
    "Xc.3v8"=Xc.3v8,
    "sdat.det.common2"=sdat.det.common2,
    "sdat.det"=sdat.det,
    "lvec1.common2"=lvec1.common2,
    "lvec2.common2"=lvec2.common2)

  saveRDS(HOR_CHP_mod_dat_sep_ls,"data/HOR_CHP_mod_dat_sep_ls.rds")
})
