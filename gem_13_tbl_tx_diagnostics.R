# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_tbl_diagnostics
# Description:  Function to take all GEM transition probability diagnostics 
# and save them in tables 

################################# Arguments ####################################

# psi44_results:
#       Diagnostics for GEM transition probability between state 4 and state 4
#       calculated from gem_calc_tx_diagnostics function.
# psi43_results:
#       Diagnostics for GEM transition probability between state 4 and state 3
#       calculated from gem_calc_tx_diagnostics function.
# psi42_results:
#       Diagnostics for GEM transition probability between state 4 and state 2
#       calculated from gem_calc_tx_diagnostics function.
# psi41_results:
#       Diagnostics for GEM transition probability between state 4 and state 1
#       calculated from gem_calc_tx_diagnostics function.
# psi33_results:
#       Diagnostics for GEM transition probability between state 3 and state 3
#       calculated from gem_calc_tx_diagnostics function.
# psi32_results:
#       Diagnostics for GEM transition probability between state 3 and state 2
#       calculated from gem_calc_tx_diagnostics function. 
# psi31_results:
#       Diagnostics for GEM transition probability between state 3 and state 1
#       calculated from gem_calc_tx_diagnostics function. 
# psi22_results:
#       Diagnostics for GEM transition probability between state 2 and state 2
#       calculated from gem_calc_tx_diagnostics function. 
# psi21_results:
#       Diagnostics for GEM transition probability between state 2 and state 1
#       calculated from gem_calc_tx_diagnostics function. 

################################## Output ######################################

# Full tables and summary tables for GEM transition probabilities and diagnostics.

################################# Function #####################################

gem_tbl_txd <- function(psi44_results, psi43_results, psi42_results, psi41_results,
                        psi33_results, psi32_results, psi31_results, psi22_results,
                        psi21_results){

  # Psi results summary table
  summarypsi <-as.data.frame(matrix(data=NA,nrow=9,ncol=5))
  colnames(summarypsi) <-c("variable","mape","coverage","rrmse","fp")
  summarypsi$variable <-c("psi44","psi43","psi42","psi41","psi33","psi32","psi31","psi22","psi21")
  
  summarypsi[1,2] <-mean(psi44_results[psi44_results$mape !=999,14])
  summarypsi[1,3] <-(sum(psi44_results$coverage)/nrow(psi44_results))*100
  summarypsi[1,4] <-mean(psi44_results[psi44_results$mape !=999,15])
  summarypsi[1,5] <- (nrow(psi44_results[psi44_results$mape ==999,])/(nrow(psi44_results)))*100
  
  summarypsi[2,2] <-mean(psi43_results[psi43_results$mape !=999,14])
  summarypsi[2,3] <-(sum(psi43_results$coverage)/nrow(psi43_results))*100
  summarypsi[2,4] <-mean(psi43_results[psi43_results$mape !=999,15])
  summarypsi[2,5] <- (nrow(psi43_results[psi43_results$mape ==999,])/(nrow(psi43_results)))*100
  
  summarypsi[3,2] <-mean(psi42_results[psi42_results$mape !=999,14])
  summarypsi[3,3] <-(sum(psi42_results$coverage)/nrow(psi42_results))*100
  summarypsi[3,4] <-mean(psi42_results[psi42_results$mape !=999,15])
  summarypsi[3,5] <- (nrow(psi42_results[psi42_results$mape ==999,])/(nrow(psi42_results)))*100
  
  summarypsi[4,2] <-mean(psi44_results[psi44_results$mape !=999,14])
  summarypsi[4,3] <-(sum(psi44_results$coverage)/nrow(psi44_results))*100
  summarypsi[4,4] <-mean(psi44_results[psi44_results$mape !=999,15])
  summarypsi[4,5] <- (nrow(psi44_results[psi44_results$mape ==999,])/(nrow(psi44_results)))*100
  
  # summarypsi[5,2] <-mean(psi33_results[psi33_results$mape !=999,13])
  # summarypsi[5,3] <-(sum(psi33_results$coverage)/nrow(psi33_results))*100
  # summarypsi[5,4] <-mean(psi33_results[psi33_results$mape !=999,14])
  # summarypsi[5,5] <- (nrow(psi33_results[psi33_results$mape ==999,])/(nrow(psi33_results)))*100
  # 
  # summarypsi[6,2] <-mean(psi32_results[psi32_results$mape !=999,13])
  # summarypsi[6,3] <-(sum(psi32_results$coverage)/nrow(psi32_results))*100
  # summarypsi[6,4] <-mean(psi32_results[psi32_results$mape !=999,14])
  # summarypsi[6,5] <- (nrow(psi32_results[psi32_results$mape ==999,])/(nrow(psi32_results)))*100
  # 
  # summarypsi[7,2] <-mean(psi31_results[psi31_results$mape !=999,13])
  # summarypsi[7,3] <-(sum(psi31_results$coverage)/nrow(psi31_results))*100
  # summarypsi[7,4] <-mean(psi31_results[psi31_results$mape !=999,14])
  # summarypsi[7,5] <- (nrow(psi31_results[psi31_results$mape ==999,])/(nrow(psi31_results)))*100
  # 
  # summarypsi[8,2] <-mean(psi22_results[psi22_results$mape !=999,12])
  # summarypsi[8,3] <-(sum(psi22_results$coverage)/nrow(psi22_results))*100
  # summarypsi[8,4] <-mean(psi22_results[psi22_results$mape !=999,13])
  # summarypsi[8,5] <-(nrow(psi22_results[psi22_results$mape ==999,])/(nrow(psi22_results)))*100
  # 
  # summarypsi[9,2] <-mean(psi21_results[psi21_results$mape !=999,12])
  # summarypsi[9,3] <-(sum(psi21_results$coverage)/nrow(psi21_results))*100
  # summarypsi[9,4] <-mean(psi21_results[psi21_results$mape !=999,13])
  # summarypsi[9,5] <-(nrow(psi21_results[psi21_results$mape ==999,])/(nrow(psi21_results)))*100
  # 
  write.csv(summarypsi, file = "./psi_results_summary.csv")
  
  # Psi results full
  write.csv(psi44_results, file = "./psi44_results_all.csv")
  write.csv(psi43_results, file = "./psi43_results_all.csv")
  write.csv(psi42_results, file = "./psi42_results_all.csv")
  write.csv(psi41_results, file = "./psi41_results_all.csv")
  # write.csv(psi33_results, paste0(filename, "psi33_results_all.csv"))
  # write.csv(psi32_results, paste0(filename, "psi32_results_all.csv"))
  # write.csv(psi31_results, paste0(filename, "psi31_results_all.csv"))
  # write.csv(psi22_results, paste0(filename, "psi22_results.csv"))
  # write.csv(psi21_results, paste0(filename, "psi21_results.csv"))

}
