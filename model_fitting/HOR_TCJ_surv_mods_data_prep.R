############################################################################### #
# HOR-TCJ survival models: reduce original model output to a data frame (kf.df)
# Input:  (searched via input_dirs in config.R) Multiyear v2 extra SJL_TCJ.RData
# Output: data/HOR_TCJ_survmod_data.RData (kf.df), read by HOR_TCJ_cand_mods_run.R
# Run from the TMB_SD_sthd root.
############################################################################### #

if (!exists("input_path")) source("model_fitting/config.R")

# local(): the large input workspace does not stay in the global environment
# (the *_cand_mods_run.R scripts save the global environment to out_dir)
local({

load(input_path("Multiyear v2 extra SJL_TCJ.RData"))

kf.df=ibm.df[!is.na(ibm.df$HOR.site) & ibm.df$HOR.site=="SJL" & ibm.df$YrRel!="2014-S1" & !ibm.df$censor,]
kf.df$drought=FALSE
kf.df$drought[kf.df$year %in% c(2013:2015)]=TRUE
kf.df$drought.fac=factor(kf.df$drought)
kf.df$Year=factor(kf.df$year)
kf.df$YrRel=gsub("S","",kf.df$YrRel)
kf.df$YrRel=factor(as.character(kf.df$YrRel))
kf.df$detd=(kf.df$TCJ==1)
kf.df$pMAC=NA
kf.df$pTRN=NA
for(i in kf.df$YrRel) {kf.df$pMAC[kf.df$YrRel==i]=par.df$Estimate[par.df$YrRel==i & par.df$parameter=="pMAC"]; kf.df$pTRN[kf.df$YrRel==i]=par.df$Estimate[par.df$YrRel==i & par.df$parameter=="pTRN"]}

kf.df$log.VNS.hor.5=log(kf.df$VNS.hor.5)
kf.df$log.IE.hor.5=log(kf.df$IE.hor.5)

WYT_year_ref <- data.frame(
  year=c(2011,2012,2013,2014,2015,2016),
  WYT=c("wet","dry","drought","drought","drought","dry"),
  SJ_CDEC_WYT=c("Wet","Dry","Critical","Critical","Critical","Dry"))

kf.df$WYT <- factor(WYT_year_ref$WYT[match(kf.df$year,WYT_year_ref$year)],levels = c("dry","drought","wet"))

save(list = c("kf.df"),file = "data/HOR_TCJ_survmod_data.RData")

})
