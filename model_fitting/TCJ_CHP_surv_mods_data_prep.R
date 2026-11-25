############################################################################### #
# TCJ-CHP survival model: data preparation and design matrices
#
# Two parts:
#   1. Function definition (always loaded): the legacy TCJ_CHP_comp(). The
#      pipeline no longer calls it: the design matrix, TMB data list and candidate
#      model strings (allint_DMs) are built by TCJ_CHP_precomp() and
#      TCJ_CHP_modelprep() in R/TCJ_CHP_comp_ALT.R, with the formula-based design
#      in R/TCJ_CHP_designmat.R, R/HOR_CHP_designmat.R and R/legacy_adapters.R
#      (see model_fitting/README.md and the reference/ folders).
#   2. Data preparation (runs only if run_data_prep is TRUE, the default):
#      reduces the original Buchanan et al. (2021) workspace to the objects that
#      TCJ_CHP_precomp() reads.
#
# Current design (setting CHP_wyt in config.R): the 4 always-in terms are
#   (Intercept), WYT_wet, WYT_drought and TCJroute.facF; the 9 toggled terms are
#   flength, Qomt.tcj.1net, log.VNS.tcj, CVP.tcj.4, SWP.tcj.4, Tmsd.tcj.7dadm and
#   route x {SWP, CVP, OMT}, minus models with a route interaction but not its
#   main effect (216 of 512 remain). WYT_ follows each row's water.yr (the legacy
#   design used a release-group lookup for WYT_drought).
#
# Part 2 inputs (searched via input_dirs in config.R): tmp_img4FULL.Rdata, lvec.tcj_ls.rds
# Part 2 outputs: data/TCJ_CHP_mod_dat_ls.RData
# Run from the TMB_SD_sthd root.
#
# ----------------------------------------------------------------------------- #
# LEGACY: VARIABLE SUBSTITUTIONS AND COLUMNS USED BY TCJ_CHP_comp()
# (kept for reference; the description below is of the legacy function)
# ----------------------------------------------------------------------------- #
# Rows:       x.df rows with a non-missing TCJ.det
# Route:      TCJroute.fac is built from the MAC / TRN indicators ("A" = MAC,
#             "F" = TRN); the design-matrix column is TCJroute.facF
# Kept as columns of x.df_SUB (TCJ_CHP_VARS):
#             flength, log.QOUT, Qomt.tcj.1net, VNS.tcj, CVP.tcj.4, SWP.tcj.4,
#             Tmsd.tcj.7dadm, Tclc.tcj.3
# Design matrix columns: (Intercept), year.fac2012-2016 (2011 is the baseline),
#             TCJroute.facF (7 columns, BASIC_IND_COLS = 1:7)
# Substitutions:
#   - year effects are dropped (columns 2-6) and replaced by WYT_wet and
#     WYT_drought indicators derived from each release group's water-year type
#   - VNS.tcj is renamed log.VNS.tcj and then log-transformed (the rename must
#     happen first; the log is taken of the renamed column)
#   - missing Qomt.tcj.1net values are replaced by 0 before scaling (z-scaled
#     runs only)
#   - main effects (MAIN_EFFECT_IDS): flength, Qomt.tcj.1net, VNS.tcj, CVP.tcj.4,
#     SWP.tcj.4, Tmsd.tcj.7dadm
#   - Tclc.tcj.3 and log.QOUT are kept in the data but are not main effects;
#     log.QOUT is used only for detection (P_col_ids = columns 1 and 15)
#   - z_scale_vars = TRUE centres and scales flength, log.VNS.tcj,
#     Qomt.tcj.1net, CVP.tcj.4, SWP.tcj.4 and Tmsd.tcj.7dadm
# Interactions: route x {SWP, CVP, OMT} only (R_x_SWP, R_x_CVP, R_x_OMT); the
#             route x Temp and route x VNS terms are commented out in the code
# Other:      n_relgrps = 18 is hard-coded in TMB_data_baseline
# Candidate set (allint_DMs): 4 fixed terms plus 9 toggled terms (2^9 = 512),
#             minus models with an interaction but not its main effect
#             (216 models remain)
############################################################################### #

if (!exists("run_data_prep")) run_data_prep <- TRUE

# ---- Part 1: function ------------------------------------------------------- #

#' Title
#'
#' @param RData_pth_in  path of model data
#' @param z_scale_vars logical, scale and center variables or not
#'
#' @returns named list with vital 
#' 
#'
TCJ_CHP_comp <- function(
  sdat.det.in,
  lvec.tcj_ls_in,
    x.df
    # RData_pth_in="data/HOR_CHP_mod_dat_ls.RData",
    ,
    z_scale_vars=TRUE
    ){
  # load(RData_pth_in)

  lvec1.tcj <- lvec.tcj_ls_in$lvec1.tcj
  lvec2.tcj <- lvec.tcj_ls_in$lvec2.tcj

  # ommitting NAs not actually used
  x.df_rem1 <- subset(x.df,#tag!="1232531" &
                      !is.na(det.pred.qomt.rms))

  # for matching up inputs
  # tmpDF <- x.df_rem1 %>%
  #   select(tag,year,rel.grp,HOR.det, ,flength,barrier.fac) %>%
  #   mutate(rel_grp=paste(year,rel.grp),
  #          model="survival",
  #          foc_area="HOR-CHP")

  x.df_rem1 <- subset(x.df,!is.na(x.df$TCJ.det) )#& !is.na(det.pred.qomt.rms))
  
  WYT_ref_tab <- table(paste(x.df_rem1$year,x.df_rem1$rel.grp),x.df_rem1$water.yr)
  nrow(x.df_rem1)
  sdat.det.tcj <- as.matrix(x.df_rem1[,match(colnames(sdat.det.in),names(x.df_rem1))])
  
  
  # 
  BASE_COLS <- c("form","rel.grp","year","water.yr","route.fac","barrier")
  # HOR_CHP_VARS  <- c(
  #   "flength",         # Fork length
  #   "log.QOUT",        # OUTFLOW
  #   "Tmsd.hor.7dadm",  # MSD temp
  #   "log.VNS.hor.5",   # Vernalis flow
  #   "CVP.hor.5",       # Exports
  #   "SWP.hor.5",       # Exports
  #   "Qomt.hor.1net",   # Flow
  #   "Tclc.hor.3",
  #   "VNS.rel",
  #   "VNS.hor"
  # )      # CLC temp
  TCJ_CHP_VARS  <- c(
    "flength",         # Fork length
    "log.QOUT",        # OUTFLOW
    # SW_picks
    "Qomt.tcj.1net", #don't love the "net" part
    "VNS.tcj",
    "CVP.tcj.4",
    "SWP.tcj.4",
    "Tmsd.tcj.7dadm",
    "Tclc.tcj.3"
  )      # CLC temp
  
  x.df_rem1$TCJroute=factor(apply(x.df_rem1[,c("MAC","TRN")],1,function(x) c("A","F")[which(x==1)]))
  x.df_rem1$TCJroute.fac=factor(x.df_rem1$TCJroute)
  x.df_rem1$day.TCJ.fac=factor(x.df_rem1$day.TCJ)
  
  
  # x.df_rem1
  X.tcjchp.1=model.matrix(rel.grp~log.QOUT + year.fac + TCJroute.fac + barrier.fac + flength + drought.fac, #day.TCJ.fac + night.TCJ.fac + tod.TCJ.fac + 
                          data=x.df_rem1)
  Xtc.1=X.tcjchp.1
  
  x.df_SUB <- cbind(Xtc.1[,
                          c("(Intercept)","year.fac2012","year.fac2013","year.fac2014","year.fac2015","year.fac2016","TCJroute.facF"
  )],
  x.df_rem1[,c(BASE_COLS,TCJ_CHP_VARS)])

  
  rw_ind <- 1:nrow(x.df_SUB)
  
  # 
  XX_in <- x.df_SUB

  
  # repl_tab <- data.frame(
  #   pattrn=c("X.Intercept.","route.facB.barrier.facTRUE"),
  #   rplace=c("(Intercept)","route.facB:barrier.facTRUE"))
  # 
  # scl_var_nms <- c("flength","Tmsd.hor.7dadm","log.VNS.hor.5","SWP.hor.5" ,"CVP.hor.5","Qomt.hor.1net","Tclc.hor.3")
  scl_var_nms <- c(
    "flength",
    "Qomt.tcj.1net", #don't love the "net" part
    "VNS.tcj",
    "CVP.tcj.4",
    "SWP.tcj.4",
    "Tmsd.tcj.7dadm",
    "Tclc.tcj.3"
  )




  # vector of  4 linear combinations used in only for loop in CPP template
  P_lc_n <- length(c("P1_f1","Lam_f1",
                     "P11_f2f3","P12_f2f3"))
  # row_inds for different likelihood forms
  f1 <- which(XX_in[,"form"]==1)
  f2 <- which(XX_in[,"form"]==2)
  f3 <- which(XX_in[,"form"]==3)

  BASIC_IND_COLS <- c(1:7)
  # MAIN_EFFECT_IDS <- match(c("flength","Tmsd.hor.7dadm","log.VNS.hor.5","SWP.hor.5","CVP.hor.5","Qomt.hor.1net"),names(XX_in))
  MAIN_EFFECT_IDS <- match(c(
    "flength",
    "Qomt.tcj.1net", #don't love the "net" part
    "VNS.tcj",
    "CVP.tcj.4",
    "SWP.tcj.4",
    "Tmsd.tcj.7dadm"
    ),names(XX_in))
  S_col_add_ids <- c(BASIC_IND_COLS,MAIN_EFFECT_IDS)

  # # additive design matrix
  
  # log transform vns a head of everything
  colnames(XX_in)[colnames(XX_in)=="VNS.tcj"]="log.VNS.tcj"
  XX_in[,"log.VNS.tcj"] <- log(XX_in[,"log.VNS.tcj"])
  
  if(z_scale_vars){
    
    # sapply()
    
    XX_in[,"flength"] <- scale(XX_in[,"flength"])
    XX_in[,"log.VNS.tcj"] <- scale(XX_in[,"log.VNS.tcj"])
    XX_in[,"Qomt.tcj.1net"] <- ifelse(is.na(XX_in[,"Qomt.tcj.1net"]),0,XX_in[,"Qomt.tcj.1net"]) # substitution of NAs
    XX_in[,"Qomt.tcj.1net"] <- scale(XX_in[,"Qomt.tcj.1net"])
    XX_in[,"CVP.tcj.4"] <- scale(XX_in[,"CVP.tcj.4"])
    XX_in[,"SWP.tcj.4"] <- scale(XX_in[,"SWP.tcj.4"])
    XX_in[,"Tmsd.tcj.7dadm"] <- scale(XX_in[,"Tmsd.tcj.7dadm"])
    
    # attributes()

    }
  
  # if(z_scale_vars){
  #   XX_in <- CVhelp::scale_data_cols(
  #     scale_cols = scl_var_nms,
  #     repl_tab_in=repl_tab,
  #     x.df_SUB_in=x.df_SUB)
  # } else{
  #   XX_in <- x.df_SUB
  #   var_scl_info_tmp <- list("scaled_vars" = scl_var_nms,
  #                            "center" = rep(0,length(scl_var_nms)),
  #                            "scale" =  rep(1,length(scl_var_nms)))
  #   attributes(XX_in) <- c(attributes(XX_in), var_scl_info_tmp)
  # }
  XX_in_add <- XX_in[,S_col_add_ids]
  
  
  
  # XX_in_add[,"Tclc.tcj.3"] <- scale(XX_in_add[,"Tclc.tcj.3"])
  
  Rte_col_id <- 7
  # 
  # # adding one-way interactions design matrix
  XX_in_w_int <- cbind(XX_in_add,
                       data.frame(
                        # Route interactions
                         # R_x_Temp=XX_in_add[,Rte_col_id]*XX_in_add[,c("Tmsd.tcj.7dadm")], # does not take column from main effect
                         # R_x_VNS=XX_in_add[,Rte_col_id]*XX_in_add[,c("log.VNS.tcj")],
                         R_x_SWP=XX_in_add[,Rte_col_id]*XX_in_add[,c("SWP.tcj.4")],
                         R_x_CVP=XX_in_add[,Rte_col_id]*XX_in_add[,c("CVP.tcj.4")],
                         R_x_OMT=XX_in_add[,Rte_col_id]*XX_in_add[,c("Qomt.tcj.1net")]
                       )
  )
  # 
  # # columns used to model survival (intercept and outflow)
  P_col_ids <- c(1,15)
  # 
  rel_grps_df_tmp <- data.frame(
    rel_grps=unique(paste(XX_in$year,XX_in$rel.grp)),
    t(sapply(unique(paste(XX_in$year,XX_in$rel.grp)),function(x) strsplit(x,split = " ")[[1]])))
  rel_grps_df_tmp[,2]=as.numeric(rel_grps_df_tmp[,2])
  names(rel_grps_df_tmp)[c(2,3)]=c("year","rel_nm")

  rel_grps_df_tmp <- rel_grps_df_tmp[order(rel_grps_df_tmp$year,rel_grps_df_tmp$rel_nm),]
  rel_grps_df_tmp$relgrpF <- factor(rel_grps_df_tmp$rel_grps,levels = rel_grps_df_tmp$rel_grps)
  
  rel_grps_df_tmp
  

  tmp_tab <- which(WYT_ref_tab!=0,arr.ind = T)
  rel_grps_df <- data.frame(rel_grps_df_tmp,
                            WYT=c("wet","dry","drought")[tmp_tab[tmp_tab[,"row"],"col"]],
                            R_ind=as.numeric(rel_grps_df_tmp$relgrpF),
                            TMB_ind=as.numeric(rel_grps_df_tmp$relgrpF)-1)

  relgrps_ind_IN <- rel_grps_df$TMB_ind[match(paste(XX_in$year,XX_in$rel.grp),rel_grps_df$rel_grps)]

  WYT_wet=(relgrps_ind_IN %in% rel_grps_df$TMB_ind[rel_grps_df$WYT=="wet"])+0
  WYT_drought=(relgrps_ind_IN %in% rel_grps_df$TMB_ind[rel_grps_df$WYT=="drought"])+0

  XX_in_w_int_noYRS <- XX_in_w_int[-c(2:6)]
  XX_in_w_int_WYT <- cbind(XX_in_w_int_noYRS[,"(Intercept)"],
                           WYT_wet,WYT_drought,
                           XX_in_w_int_noYRS[,which(names(XX_in_w_int_noYRS)!="(Intercept)")])
  names(XX_in_w_int_WYT)[1] <- "(Intercept)"

  # list("WYT_ref_tab"=WYT_ref_tab,
  #      "rel_grps_df_tmp"=rel_grps_df_tmp,
  #      "XX_in_w_int_WYT"=XX_in_w_int_WYT)
  # all_potential_cols2include <- c("(Intercept)","year.fac2012","year.fac2013","year.fac2014","year.fac2015","year.fac2016","TCJroute.facF","flength","Qomt.tcj.1net","VNS.tcj","CVP.tcj.4","SWP.tcj.4","Tmsd.tcj.7dadm","Tclc.tcj.3")
  # S_col_ids <- match(all_potential_cols2include,colnames(XX_in))
  
  # S_col_ids <- 1:9
  S_col_ids <- nrow(XX_in_w_int_WYT)

  # XX_in
  # 
  # 
  TMB_data_baseline <- list(
    "n_relgrps"=18,
    "relgrps_ind"=relgrps_ind_IN,
    "XX"=as.matrix(XX_in[,P_col_ids]),
    "XX_s"=as.matrix(XX_in_w_int_WYT),
    "Dat"=sdat.det.tcj,
    "col_Dat1"=match(paste(c("CHP.BBR"),c("11","10","01","00"),sep="."),colnames(sdat.det.tcj)) - 1,
    "col_Dat2"=match(paste(c("CHP.2x"),c("11","10","01","00"),sep="."),colnames(sdat.det.tcj)) - 1,
    "col_Dat3"=match(paste(c("CHP.3x"),c("111","110","101","100","011","010","001","000"),sep="."),colnames(sdat.det.tcj))-1,
    "ind_all"=1:length(rw_ind)-1,
    "f1_ind"=f1-1,
    "f2_ind"=f2-1,
    "f3_ind"=f3-1,
    "lvec1"=lvec1.tcj,
    "lvec2"=lvec2.tcj,
    # "col_S"=S_col_ids-1,
    "col_S"=(1:length(S_col_ids))-1,
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

  allint_DMs_raw <- apply(cbind(rep(c("1111"), 2^9), expand.grid(c(0, 1), c(0, 1), c(0, 1),c(0, 1), c(0, 1), c(0, 1),c(0, 1), c(0, 1), c(0, 1)))
                          , 1, paste0, collapse = "")
  DM_mat_raw<- do.call(rbind,lapply(1:length(allint_DMs_raw),function(x) as.numeric(strsplit(allint_DMs_raw[x],split="")[[1]]) ))
  colnames(DM_mat_raw) <- colnames(TMB_data_baseline$XX_s)
  
  bad_ind <- which(DM_mat_raw[,"Qomt.tcj.1net"]==0 & DM_mat_raw[,"R_x_OMT"]==1 | DM_mat_raw[,"CVP.tcj.4"]==0 & DM_mat_raw[,"R_x_CVP"]==1 | DM_mat_raw[,"SWP.tcj.4"]==0 & DM_mat_raw[,"R_x_SWP"]==1)
  allint_DMs_raw[bad_ind]
  DM_mat <- DM_mat_raw[-c(bad_ind),]
  allint_DMs <- allint_DMs_raw[-c(bad_ind)]
  
  list("allint_DMs"=allint_DMs,
       "WYT_ref_tab"=WYT_ref_tab,
       "rel_grps_df_tmp"=rel_grps_df_tmp,
       "XX_in_add"=XX_in_add,
       "XX_in_w_int_WYT"=XX_in_w_int_WYT,
       "TMB_data_baseline"=TMB_data_baseline,
       "XX_in"=XX_in)
  
  
}


# ---- Part 2: data preparation ----------------------------------------------- #

# local(): the loaded workspace does not stay in the global environment
# (the *_cand_mods_run.R scripts save the global environment to out_dir)
if (run_data_prep) local({

  if (!exists("input_path")) {
    source("model_fitting/config.R")
    source("R/input_path.R")
  }
  load(input_path("tmp_img4FULL.Rdata"))
  lvec.tcj_ls <- readRDS(input_path("lvec.tcj_ls.rds"))

  objs_for_surv_TCJ_CHP_mods <- c("x.df","lvec.tcj_ls","sdat.det.common2")
  save(list=objs_for_surv_TCJ_CHP_mods,file="data/TCJ_CHP_mod_dat_ls.RData")
})
