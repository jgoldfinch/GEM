# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_tbl_diagnostics
# Description:  Function to take all sim restuls diagnostics and put them in 
# tables
# RRMSE - relative root mean square error
# MAPE - mean absolute percent error
# Coverage - percentage of simulations the credible interval covers truth
# Fp - false positive rate

################################# Arguments ####################################

# p_results:
#       Detection probability diagnostics produced with gem_comp_diagnostics 
#       function.
# pg_results:
#       Genetic sign detection probability diagnostics produced with 
#       gem_comp_diagnostics function.
# nf_results:
#       Adult female predictions diagnostics produced with gem_comp_diagnostics 
#       function.
# nm_results:
#       Adult male predictions diagnostics produced with gem_comp_diagnostics 
#       function.
# be_results:
#       Birth event predictions diagnostics produced with gem_comp_diagnostics 
#       function.
# wf_results:
#       New female predictions diagnostics produced with gem_comp_diagnostics 
#       function.
# wm_results:
#       New male predictions diagnostics produced with gem_comp_diagnostics 
#       function.
# wnf_results:
#       Total females (adult and new) alive at the end of a time step diagnostics 
#       produced with gem_comp_diagnostics function.
# wnm_results:
#       Total males (adult and new) alive at the end of a time step diagnostics 
#       produced with gem_comp_diagnostics function.
# wtot_results:
#       Total individuals (adult and new) alive at the end of a time step 
#       diagnostics produced with gem_comp_diagnostics function.
# z_results:
#       Population occupancy at the end of each time step diagnostics produced 
#       with gem_comp_diagnostics function.
# surv_results:
#       Survival diagnostics produced with gem_comp_diagnostics function.

################################## Output ######################################

# Full tables and summary tables for adult females (nf), adult males (nm), birth events (be), 
# new females (wf), new males (wm), total adult and new females (wnf), total 
# adult and new males (wnm), total adult and new individuals (wtot), population 
# occupancy (z), survival (surv), detection probability (p), and
# detection probability of genetic sign (pg).

################################# Function #####################################

gem_tbl_ds <- function(p_results, pg_results, nf_results, nm_results, be_results, wf_results,
                       wm_results, wnf_results, wnm_results, wtot_results, 
                       z_results, s_results){

  filename = paste0 ("./")
  
  # Note that 999 was used in the previous function (gem_comp_diagnostics) to 
  # denote time steps where the population was 0, so it is used to filter out rows 
  # in this function for calculating diagnostics
  
  # NF results
  write.csv(nf_results, paste0(filename, "nf_results_all.csv"))
  nf_results_summary <- as.data.frame(cbind(sum(nf_results$convergance)/nrow(nf_results), 
                              sum(nf_results[nf_results$mape !=999,11])/nrow(nf_results[nf_results$mape !=999,]),
                              sum(nf_results[nf_results$mape !=999,12])/nrow(nf_results[nf_results$mape !=999,]),
                              (sum(nf_results$coverage)/nrow(nf_results))*100,
                              sum(nf_results[nf_results$mape !=999,13])/nrow(nf_results[nf_results$mape !=999,]),
                              (sum(nrow(nf_results[nf_results$mape ==999,])/nrow(nf_results)))*100))
  colnames(nf_results_summary) <-c("convergence","rrmse","mape","coverage","maie", "fp")
  nf_results_summary$var <-"nf"
  write.csv(nf_results_summary, paste0(filename, "nf_results_summary.csv"))
  
  # NM results
  write.csv(nm_results, paste0(filename, "nm_results_all.csv"))
  nm_results_summary <- as.data.frame(cbind(sum(nm_results$convergance)/nrow(nm_results), 
                                            sum(nm_results[nm_results$mape !=999,11])/nrow(nm_results[nm_results$mape !=999,]),
                                            sum(nm_results[nm_results$mape !=999,12])/nrow(nm_results[nm_results$mape !=999,]),
                                            (sum(nm_results$coverage)/nrow(nm_results))*100,
                                            sum(nm_results[nm_results$mape !=999,13])/nrow(nm_results[nm_results$mape !=999,]),
                                            (sum(nrow(nm_results[nm_results$mape ==999,])/nrow(nm_results)))*100))
  colnames(nm_results_summary) <-c("convergence","rrmse","mape","coverage","maie", "fp")
  nm_results_summary$var <-"nm"
  write.csv(nm_results_summary, paste0(filename, "nm_results_summary.csv"))
  
  # BE results
  write.csv(be_results, paste0(filename, "be_results_all.csv"))
  be_results_summary <- as.data.frame(cbind(sum(be_results$convergance)/nrow(be_results), 
                                            sum(be_results[be_results$mape !=999,11])/nrow(be_results[be_results$mape !=999,]),
                                            sum(be_results[be_results$mape !=999,12])/nrow(be_results[be_results$mape !=999,]),
                                            (sum(be_results$coverage)/nrow(be_results))*100,
                                            sum(be_results[be_results$mape !=999,13])/nrow(be_results[be_results$mape !=999,]),
                                            (sum(nrow(be_results[be_results$mape ==999,])/nrow(be_results)))*100))
  colnames(be_results_summary) <-c("convergence","rrmse","mape","coverage","maie", "fp")
  be_results_summary$var <-"be"
  write.csv(be_results_summary, paste0(filename, "be_results_summary.csv"))
  
  # WF results
  write.csv(wf_results, paste0(filename, "wf_results_all.csv"))
  wf_results_summary <- as.data.frame(cbind(sum(wf_results$convergance)/nrow(wf_results), 
                                            sum(wf_results[wf_results$mape !=999,11])/nrow(wf_results[wf_results$mape !=999,]),
                                            sum(wf_results[wf_results$mape !=999,12])/nrow(wf_results[wf_results$mape !=999,]),
                                            (sum(wf_results$coverage)/nrow(wf_results))*100,
                                            sum(wf_results[wf_results$mape !=999,13])/nrow(wf_results[wf_results$mape !=999,]),
                                            (sum(nrow(wf_results[wf_results$mape ==999,])/nrow(wf_results)))*100))
  colnames(wf_results_summary) <-c("convergence","rrmse","mape","coverage","maie", "fp")
  wf_results_summary$var <-"wf"
  write.csv(wf_results_summary, paste0(filename, "wf_results_summary.csv"))

  # WM results
  write.csv(wm_results, paste0(filename, "wm_results_all.csv"))
  wm_results_summary <- as.data.frame(cbind(sum(wm_results$convergance)/nrow(wm_results), 
                                              sum(wm_results[wm_results$mape !=999,11])/nrow(wm_results[wm_results$mape !=999,]),
                                              sum(wm_results[wm_results$mape !=999,12])/nrow(wm_results[wm_results$mape !=999,]),
                                              (sum(wm_results$coverage)/nrow(wm_results))*100,
                                              sum(wm_results[wm_results$mape !=999,13])/nrow(wm_results[wm_results$mape !=999,]),
                                              (sum(nrow(wm_results[wm_results$mape ==999,])/nrow(wm_results)))*100))
  colnames(wm_results_summary) <-c("convergence","rrmse","mape","coverage","maie", "fp")
  wm_results_summary$var <-"wm"
  write.csv(wm_results_summary, paste0(filename, "wm_results_summary.csv"))
  
  # WNF results
  write.csv(wnf_results, paste0(filename, "wnf_results_all.csv"))
  wnf_results_summary <- as.data.frame(cbind(sum(wnf_results$convergance)/nrow(wnf_results), 
                                            sum(wnf_results[wnf_results$mape !=999,11])/nrow(wnf_results[wnf_results$mape !=999,]),
                                            sum(wnf_results[wnf_results$mape !=999,12])/nrow(wnf_results[wnf_results$mape !=999,]),
                                            (sum(wnf_results$coverage)/nrow(wnf_results))*100,
                                            sum(wnf_results[wnf_results$mape !=999,13])/nrow(wnf_results[wnf_results$mape !=999,]),
                                            (sum(nrow(wnf_results[wnf_results$mape ==999,])/nrow(wnf_results)))*100))
  colnames(wnf_results_summary) <-c("convergence","rrmse","mape","coverage","maie", "fp")
  wnf_results_summary$var <-"wnf"
  write.csv(wnf_results_summary, paste0(filename, "wnf_results_summary.csv"))
  
  # WNM results
  write.csv(wnm_results, paste0(filename, "wnm_results_all.csv"))
  wnm_results_summary <- as.data.frame(cbind(sum(wnm_results$convergance)/nrow(wnm_results), 
                                             sum(wnm_results[wnm_results$mape !=999,11])/nrow(wnm_results[wnm_results$mape !=999,]),
                                             sum(wnm_results[wnm_results$mape !=999,12])/nrow(wnm_results[wnm_results$mape !=999,]),
                                             (sum(wnm_results$coverage)/nrow(wnm_results))*100,
                                             sum(wnm_results[wnm_results$mape !=999,13])/nrow(wnm_results[wnm_results$mape !=999,]),
                                             (sum(nrow(wnm_results[wnm_results$mape ==999,])/nrow(wnm_results)))*100))
  colnames(wnm_results_summary) <-c("convergence","rrmse","mape","coverage","maie", "fp")
  wnm_results_summary$var <-"wnm"
  write.csv(wnm_results_summary, paste0(filename, "wnm_results_summary.csv"))

  # WTOT results
  write.csv(wtot_results, paste0(filename, "wtot_results_all.csv"))
  wtot_results_summary <- as.data.frame(cbind(sum(wtot_results$convergance)/nrow(wtot_results), 
                                            sum(wtot_results[wtot_results$mape !=999,11])/nrow(wtot_results[wtot_results$mape !=999,]),
                                            sum(wtot_results[wtot_results$mape !=999,12])/nrow(wtot_results[wtot_results$mape !=999,]),
                                            (sum(wtot_results$coverage)/nrow(wtot_results))*100,
                                            sum(wtot_results[wtot_results$mape !=999,13])/nrow(wtot_results[wtot_results$mape !=999,]),
                                            (sum(nrow(wtot_results[wtot_results$mape ==999,])/nrow(wtot_results)))*100))
  colnames(wtot_results_summary) <-c("convergence","rrmse","mape","coverage","maie", "fp")
  wtot_results_summary$var <-"wtot"
  write.csv(wtot_results_summary, paste0(filename, "wtot_results_summary.csv"))

  # Detection results
  write.csv(p_results, paste0(filename, "p_results_all.csv"))
  p_results_summary <- as.data.frame(cbind(sum(p_results$convergance)/nrow(p_results),
                                              sum(p_results$rrmse)/nrow(p_results),
                                              sum(p_results$mape)/nrow(p_results),
                                              (sum(p_results$coverage)/nrow(p_results))*100))
  colnames(p_results_summary) <-c("convergence","rrmse","mape","coverage")
  p_results_summary$var <-"p"
  write.csv(p_results_summary, paste0(filename, "p_results_summary.csv"))
  
  # Genetic detection results
  write.csv(pg_results, paste0(filename, "pg_results_all.csv"))
  pg_results_summary <- as.data.frame(cbind(sum(pg_results$convergance)/nrow(pg_results),
                                           sum(pg_results$rrmse)/nrow(pg_results),
                                           sum(pg_results$mape)/nrow(pg_results),
                                           (sum(pg_results$coverage)/nrow(pg_results))*100))
  colnames(pg_results_summary) <-c("convergence","rrmse","mape","coverage")
  pg_results_summary$var <-"pg"
  write.csv(pg_results_summary, paste0(filename, "pg_results_summary.csv"))

  # Survival results
  write.csv(s_results, paste0(filename, "s_results_all.csv"))
  s_results_summary <- as.data.frame(cbind(sum(s_results$convergance)/nrow(s_results),
                                           sum(s_results$rrmse)/nrow(s_results),
                                           sum(s_results$mape)/nrow(s_results),
                                           (sum(s_results$coverage)/nrow(s_results))*100))
  colnames(s_results_summary) <-c("convergence","rrmse","mape","coverage")
  s_results_summary$var <-"s"
  write.csv(s_results_summary, paste0(filename, "s_results_summary.csv"))

  # Z results
  write.csv(z_results, paste0(filename, "z_results_all.csv"))
  z_results_summary <- as.data.frame((sum(z_results$coverage)/nrow(z_results))*100)
  colnames(z_results_summary) <-c("coverage")
  z_results_summary$var <-"z"
  write.csv(z_results_summary, paste0(filename, "z_results_summary.csv"))

  # All biological results
  all_results_summary <-rbind(nf_results_summary, nm_results_summary, be_results_summary,
                              wf_results_summary, wm_results_summary, wnf_results_summary,
                              wnm_results_summary, wtot_results_summary)
  write.csv(all_results_summary, paste0(filename, "all_results_summary.csv"))
  
}
