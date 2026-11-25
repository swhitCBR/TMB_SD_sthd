TCJ_CHP_global_formula <- ~ WYT_ + TCJroute.fac + flength + Qomt.tcj.1net + log.VNS.tcj + CVP.tcj.4 + SWP.tcj.4 + Tmsd.tcj.7dadm +
  TCJroute.fac:SWP.tcj.4 + TCJroute.fac:CVP.tcj.4 + TCJroute.fac:Qomt.tcj.1net

TCJ_CHP_fixed_terms <- c("WYT_", "TCJroute.fac")
