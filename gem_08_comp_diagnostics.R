# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_comp_diagnostics
# Description:  Function to take all sim restuls and calculate diagnostics: 
# RRMSE - for comparing across parameters - which parameters are estimated well
# MAPE - mean absolute percent error - how much bias is there in estimates
# Coverage - power to estimate parameters (does credible interval cover truth)

################################# Arguments ####################################

# p_results:
#       Detection probability predictions produced with gem_run_model function.
# pg_results:
#       Genetic sign detection probability predictions produced with gem_run_model function.
# nf_results:
#       Adult female predictions produced with gem_run_model function.
# nm_results:
#       Adult male predictions produced with gem_run_model function.
# be_results:
#       Birth event predictions produced with gem_run_model function.
# wf_results:
#       New female predictions produced with gem_run_model function.
# wm_results:
#       New male predictions produced with gem_run_model function.
# wnf_results:
#       Total females (adult and new) alive at the end of a time step predictions 
#       produced with gem_run_model function.
# wnm_results:
#       Total males (adult and new) alive at the end of a time step predictions 
#       produced with gem_run_model function.
# wtot_results:
#       Total individuals (adult and new) alive at the end of a time step predictions 
#       produced with gem_run_model function.
# z_results:
#       Population occupancy at the end of each time step predictions produced 
#       with gem_run_model function.
# surv_results:
#       Survival predictions produced with gem_run_model function.

################################## Output ######################################

# diagdata: list saved to global environment. Contains model diagnostics for 
# prediction values from the formatted GEM model output for 
# detection probability (p), detection probability of genetic sign (pg), 
# adult females (nf), adult males (nm), birth events (be), 
# new females (wf), new males (wm), total adult and new females (wnf), total 
# adult and new males (wnm), total adult and new individuals (wtot), population 
# occupancy (z), and survival (surv).

################################# Function #####################################

gem_comp_diagnostics <- function(p_results, pg_results, nf_results, nm_results, 
                                 be_results, wf_results, wm_results, wnf_results, 
                                 wnm_results, wtot_results, z_results, 
                                 surv_results){
  
  # Calculate RRMSE, MAPE, COVERAGE - by parameter (across all sims and timesteps)
  # Note that 999 is used below to denote time steps where the population was 0
  
  # NF results
  nfsmall <-nf_results[complete.cases(nf_results$rhat),]
  nfsmall$convergance <-ifelse(sum(nfsmall$rhat)/(nrow(nfsmall))<=(nrow(nfsmall)*1.06),1,0)
  nfsmall$rrmse <-NA
  nfsmall$mape <-NA
  nfsmall$maie <-NA
  nfsmall$coverage <-NA
  m <-mean(nfsmall$true)
  for (i in 1:nrow(nfsmall)){
    nfsmall[i,11] <-sqrt((1/nrow(nfsmall))*sum((nfsmall[i,7]-nfsmall[i,4])^2))/m
    nfsmall[i,12] <-ifelse(nfsmall[i,7]==0, 999,((abs(nfsmall[i,7] - nfsmall[i,4]))/(nfsmall[i,7]))*100)
    nfsmall[i,13] <-ifelse(nfsmall[i,7]==0, 999, abs(nfsmall[i,7] - nfsmall[i,4]))
    nfsmall[i,14] <-ifelse(nfsmall[i,7]>=nfsmall[i,5] & nfsmall[i,7]<=nfsmall[i,6],1,0)
  }
  nfsmall$cpcnt <-(sum(nfsmall$coverage)/nrow(nfsmall))*100

  # NM results
  nmsmall <-nm_results[complete.cases(nm_results$rhat),]
  nmsmall$convergance <-ifelse(sum(nmsmall$rhat)/(nrow(nmsmall))<=(nrow(nmsmall)*1.06),1,0)
  nmsmall$rrmse <-NA
  nmsmall$mape <-NA
  nmsmall$maie <-NA
  nmsmall$coverage <-NA
  m <-mean(nmsmall$true)
  for (i in 1:nrow(nmsmall)){
    nmsmall[i,11] <-sqrt((1/nrow(nmsmall))*sum((nmsmall[i,7]-nmsmall[i,4])^2))/m
    nmsmall[i,12] <-ifelse(nmsmall[i,7]==0, 999,((abs(nmsmall[i,7] - nmsmall[i,4]))/(nmsmall[i,7]))*100)
    nmsmall[i,13] <-ifelse(nmsmall[i,7]==0, 999, abs(nmsmall[i,7] - nmsmall[i,4]))
    nmsmall[i,14] <-ifelse(nmsmall[i,7]>=nmsmall[i,5] & nmsmall[i,7]<=nmsmall[i,6],1,0)
  }
  nmsmall$cpcnt <-(sum(nmsmall$coverage)/nrow(nmsmall))*100
  
  # BE results
  besmall <-be_results[complete.cases(be_results$rhat),]
  besmall$convergance <-ifelse(sum(besmall$rhat)/(nrow(besmall))<=(nrow(besmall)*1.06),1,0)
  besmall$rrmse <-NA
  besmall$mape <-NA
  besmall$maie <-NA
  besmall$coverage <-NA
  m <-mean(besmall$true)
  for (i in 1:nrow(besmall)){
    besmall[i,11] <-sqrt((1/nrow(besmall))*sum((besmall[i,7]-besmall[i,4])^2))/m
    besmall[i,12] <-ifelse(besmall[i,7]==0, 999,((abs(besmall[i,7] - besmall[i,4]))/(besmall[i,7]))*100)
    besmall[i,13] <-ifelse(besmall[i,7]==0, 999, abs(besmall[i,7] - besmall[i,4]))
    besmall[i,14] <-ifelse(besmall[i,7]>=besmall[i,5] & besmall[i,7]<=besmall[i,6],1,0)
  }
  besmall$cpcnt <-(sum(besmall$coverage)/nrow(besmall))*100
  
  # WF results
  wfsmall <-wf_results[complete.cases(wf_results$rhat),]
  wfsmall$convergance <-ifelse(sum(wfsmall$rhat)/(nrow(wfsmall))<=(nrow(wfsmall)*1.06),1,0)
  wfsmall$rrmse <-NA
  wfsmall$mape <-NA
  wfsmall$maie <-NA
  wfsmall$coverage <-NA
  m <-mean(wfsmall$true)
  for (i in 1:nrow(wfsmall)){
    wfsmall[i,11] <-sqrt((1/nrow(wfsmall))*sum((wfsmall[i,7]-wfsmall[i,4])^2))/m
    wfsmall[i,12] <-ifelse(wfsmall[i,7]==0, 999,((abs(wfsmall[i,7] - wfsmall[i,4]))/(wfsmall[i,7]))*100)
    wfsmall[i,13] <-ifelse(wfsmall[i,7]==0, 999, abs(wfsmall[i,7] - wfsmall[i,4]))
    wfsmall[i,14] <-ifelse(wfsmall[i,7]>=wfsmall[i,5] & wfsmall[i,7]<=wfsmall[i,6],1,0)
  }
  wfsmall$cpcnt <-(sum(wfsmall$coverage)/nrow(wfsmall))*100
  
  # WM results
  wmsmall <-wm_results[complete.cases(wm_results$rhat),]
  wmsmall$convergance <-ifelse(sum(wmsmall$rhat)/(nrow(wmsmall))<=(nrow(wmsmall)*1.06),1,0)
  wmsmall$rrmse <-NA
  wmsmall$mape <-NA
  wmsmall$maie <-NA
  wmsmall$coverage <-NA
  m <-mean(wmsmall$true)
  for (i in 1:nrow(wmsmall)){
    wmsmall[i,11] <-sqrt((1/nrow(wmsmall))*sum((wmsmall[i,7]-wmsmall[i,4])^2))/m
    wmsmall[i,12] <-ifelse(wmsmall[i,7]==0, 999,((abs(wmsmall[i,7] - wmsmall[i,4]))/(wmsmall[i,7]))*100)
    wmsmall[i,13] <-ifelse(wmsmall[i,7]==0, 999, abs(wmsmall[i,7] - wmsmall[i,4]))
    wmsmall[i,14] <-ifelse(wmsmall[i,7]>=wmsmall[i,5] & wmsmall[i,7]<=wmsmall[i,6],1,0)
  }
  wmsmall$cpcnt <-(sum(wmsmall$coverage)/nrow(wmsmall))*100
  
  # WNF results
  wnfsmall <-wnf_results[complete.cases(wnf_results$rhat),]
  wnfsmall$convergance <-ifelse(sum(wnfsmall$rhat)/(nrow(wnfsmall))<=(nrow(wnfsmall)*1.06),1,0)
  wnfsmall$rrmse <-NA
  wnfsmall$mape <-NA
  wnfsmall$maie <-NA
  wnfsmall$coverage <-NA
  m <-mean(wnfsmall$true)
  for (i in 1:nrow(wnfsmall)){
    wnfsmall[i,11] <-sqrt((1/nrow(wnfsmall))*sum((wnfsmall[i,7]-wnfsmall[i,4])^2))/m
    wnfsmall[i,12] <-ifelse(wnfsmall[i,7]==0, 999,((abs(wnfsmall[i,7] - wnfsmall[i,4]))/(wnfsmall[i,7]))*100)
    wnfsmall[i,13] <-ifelse(wnfsmall[i,7]==0, 999, abs(wnfsmall[i,7] - wnfsmall[i,4]))
    wnfsmall[i,14] <-ifelse(wnfsmall[i,7]>=wnfsmall[i,5] & wnfsmall[i,7]<=wnfsmall[i,6],1,0)
  }
  wnfsmall$cpcnt <-(sum(wnfsmall$coverage)/nrow(wnfsmall))*100
  
  # wnm results
  wnmsmall <-wnm_results[complete.cases(wnm_results$rhat),]
  wnmsmall$convergance <-ifelse(sum(wnmsmall$rhat)/(nrow(wnmsmall))<=(nrow(wnmsmall)*1.06),1,0)
  wnmsmall$rrmse <-NA
  wnmsmall$mape <-NA
  wnmsmall$maie <-NA
  wnmsmall$coverage <-NA
  m <-mean(wnmsmall$true)
  for (i in 1:nrow(wnmsmall)){
    wnmsmall[i,11] <-sqrt((1/nrow(wnmsmall))*sum((wnmsmall[i,7]-wnmsmall[i,4])^2))/m
    wnmsmall[i,12] <-ifelse(wnmsmall[i,7]==0, 999,((abs(wnmsmall[i,7] - wnmsmall[i,4]))/(wnmsmall[i,7]))*100)
    wnmsmall[i,13] <-ifelse(wnmsmall[i,7]==0, 999, abs(wnmsmall[i,7] - wnmsmall[i,4]))
    wnmsmall[i,14] <-ifelse(wnmsmall[i,7]>=wnmsmall[i,5] & wnmsmall[i,7]<=wnmsmall[i,6],1,0)
  }
  wnmsmall$cpcnt <-(sum(wnmsmall$coverage)/nrow(wnmsmall))*100
  
  # WTOT results
  wtotsmall <-wtot_results[complete.cases(wtot_results$rhat),]
  wtotsmall$convergance <-ifelse(sum(wtotsmall$rhat)/(nrow(wtotsmall))<=(nrow(wtotsmall)*1.06),1,0)
  wtotsmall$rrmse <-NA
  wtotsmall$mape <-NA
  wtotsmall$maie <-NA
  wtotsmall$coverage <-NA
  m <-mean(wtotsmall$true)
  for (i in 1:nrow(wtotsmall)){
    wtotsmall[i,11] <-sqrt((1/nrow(wtotsmall))*sum((wtotsmall[i,7]-wtotsmall[i,4])^2))/m
    wtotsmall[i,12] <-ifelse(wtotsmall[i,7]==0, 999,((abs(wtotsmall[i,7] - wtotsmall[i,4]))/(wtotsmall[i,7]))*100)
    wtotsmall[i,13] <-ifelse(wtotsmall[i,7]==0, 999, abs(wtotsmall[i,7] - wtotsmall[i,4]))
    wtotsmall[i,14] <-ifelse(wtotsmall[i,7]>=wtotsmall[i,5] & wtotsmall[i,7]<=wtotsmall[i,6],1,0)
  }
  wtotsmall$cpcnt <-(sum(wtotsmall$coverage)/nrow(wtotsmall))*100

  # Z results
  zsmall <-z_results
  zsmall$coverage <-NA
  for (i in 1:nrow(zsmall)){
   zsmall[i,10] <-ifelse(zsmall[i,7]>=zsmall[i,5] & zsmall[i,7]<=zsmall[i,6],1,0)
  }
  zsmall$cpcnt <-(sum(zsmall$coverage)/nrow(zsmall))*100

  # Detection results
  psmall <-p_results
  psmall$convergance <-ifelse(sum(psmall$rhat)/(nrow(psmall))<=(nrow(psmall)*1.06),1,0)
  psmall$rrmse <-NA
  psmall$mape <-NA
  psmall$coverage <-NA
  for (i in 1:nrow(psmall)){
    psmall[i,8] <-sqrt((1/nrow(psmall))*sum((psmall[i,4]- psmall[i,1])^2))/(mean(psmall$pred))
    psmall[i,9] <-(abs(psmall[i,4] - psmall[i,1])/(psmall[i,4]))*100
    psmall[i,10] <-ifelse(psmall[i,4]>=psmall[i,2] & psmall[i,4]<=psmall[i,3],1,0)
  }
  psmall$cpcnt <-(sum(psmall$coverage)/nrow(psmall))*100
  
  # Genetic detection results
  pgsmall <-pg_results
  pgsmall$convergance <-ifelse(sum(pgsmall$rhat)/(nrow(pgsmall))<=(nrow(pgsmall)*1.06),1,0)
  pgsmall$rrmse <-NA
  pgsmall$mape <-NA
  pgsmall$coverage <-NA
  for (i in 1:nrow(pgsmall)){
    pgsmall[i,8] <-sqrt((1/nrow(pgsmall))*sum((pgsmall[i,4]- pgsmall[i,1])^2))/(mean(pgsmall$pred))
    pgsmall[i,9] <-(abs(pgsmall[i,4] - pgsmall[i,1])/(pgsmall[i,4]))*100
    pgsmall[i,10] <-ifelse(pgsmall[i,4]>=pgsmall[i,2] & pgsmall[i,4]<=pgsmall[i,3],1,0)
  }
  pgsmall$cpcnt <-(sum(pgsmall$coverage)/nrow(pgsmall))*100
  
  # Survival results
  ssmall <-surv_results[complete.cases(surv_results$rhat),]
  ssmall$convergance <-ifelse(sum(ssmall$rhat)/(nrow(ssmall))<=(nrow(ssmall)*1.06),1,0)
  ssmall$rrmse <-NA
  ssmall$mape <-NA
  ssmall$coverage <-NA
  for (i in 1:nrow(ssmall)){
    ssmall[i,9] <-sqrt((1/nrow(pgsmall))*sum((ssmall[i,4]- ssmall[i,1])^2))/(mean(ssmall$pred))
    ssmall[i,10] <-(abs(ssmall[i,4] - ssmall[i,1])/(ssmall[i,4]))*100
    ssmall[i,11] <-ifelse(ssmall[i,4]>=ssmall[i,2] & ssmall[i,4]<=ssmall[i,3],1,0)
  }
  ssmall$cpcnt <-(sum(ssmall$coverage)/nrow(ssmall))*100
  
  diagdata <- list("p_results"= psmall, "pg_results"= pgsmall, "nm_results" = nmsmall, "nf_results" = nfsmall, 
                   "be_results" = besmall, "wf_results" = wfsmall, "wm_results" = wmsmall,
                   "wnf_results" = wnfsmall, "wnm_results" = wnmsmall,
                   "wtot_results" = wtotsmall, "z_results" = zsmall, 
                   "s_results" = ssmall)
  list2env(diagdata ,.GlobalEnv)
  
}