############################################################################### #
# HOR and TCJ routing models: reduce the original route-selection workspace to
# the data frames read by HOR_cand_mods_run.R and TCJ_cand_mods_run.R
# Input:  (searched via input_dirs in config.R) Sth route selection revision_1_v4.RData
# Output: data/route_dfs_HOR_TCJ.rds, data/HOR_TCJ_route_mod_data.RData
# Run from the TMB_SD_sthd root.
############################################################################### #

if (!exists("input_path")) {
  source("model_fitting/config.R")
  source("R/input_path.R")
}

# local(): the large input workspace does not stay in the global environment
# (HOR_cand_mods_run.R and TCJ_cand_mods_run.R save the global environment to out_dir)
local({

load(input_path("Sth route selection revision_1_v4.RData"))

route_dfs_HOR_TCJ <- list("hor.df"=hor.df,"hor.stn.df"=hor.stn.df,
                          "tcj.df"=tcj.df,"tcj.stn.df"=tcj.stn.df,
                          "d.df"=d.df)

saveRDS(route_dfs_HOR_TCJ,"data/route_dfs_HOR_TCJ.rds")
save(list = c("route_dfs_HOR_TCJ"),file = "data/HOR_TCJ_route_mod_data.RData")

})


