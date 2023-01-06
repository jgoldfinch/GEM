# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_diag_tx
# Description:  function to calculate GEM transition probability diagnostics
# RRMSE - relative root mean square error
# MAPE - mean absolute percent error
# Coverage - percentage of simulations the credible interval covers truth

################################# Arguments ####################################

# tprob:
#       true occupancy of the population calculated from gem_format_tx function.
# twtot:
#       true total individuals at the end of each time step calculated from 
#       gem_format_tx function.
# twnf:
#       true total females (new and existing adult) at the end of each time step 
#       calculated from gem_format_tx function.
# twnm:
#       true total males (new and existing adult) at the end of each time step 
#       calculated from gem_format_tx function.
# psi21p:
#       predicted GEM transition probability between state 2 and state 1
# psi22p:
#       predicted GEM transition probability for staying in state 2
# psi31p:
#       predicted GEM transition probability between state 3 and state 1
# psi32p:
#       predicted GEM transition probability between state 3 and state 2
# psi33p:
#       predicted GEM transition probability between state 3 and state 3
# psi41p:
#       predicted GEM transition probability between state 4 and state 1
# psi42p:
#       predicted GEM transition probability between state 4 and state 2
# psi43p:
#       predicted GEM transition probability between state 4 and state 3
# psi44p:
#       predicted GEM transition probability between state 4 and state 4

################################## Output ######################################

# txcalcdata: list saved to global environment. Contains transition probability
# diagnostics including RRMSE, MAPE, and coverage.

################################# Function #####################################

gem_calc_txd <- function(tprob, twtot, twnf, twnm, psi21p, psi22p, psi31p, 
                         psi32p, psi33p, psi41p, psi42p, psi43p, psi44p){
   
# Calculate true probabilities
   
   true_state <-melt(tprob)
   colnames(true_state) <-c("tstate","t","rep")
   true_state$prob2 <-ifelse(true_state$tstate ==2,1,0)
   true_state$prob3 <-ifelse(true_state$tstate ==3,1,0)
   true_state$prob4 <-ifelse(true_state$tstate ==4,1,0)
   
   true_wtot <-melt(twtot)
   colnames(true_wtot) <-c("twtot","t","rep")
   true_wf <-melt(twnf)
   colnames(true_wf) <-c("twnf","t","rep")
   true_wm <-melt(twnm)
   colnames(true_wm) <-c("twnm","t","rep")
   true_wtot <-cbind(true_wtot[,2:3], true_wtot[,1], true_wf$twnf, true_wm$twnm)
   true_wtot$s <-0.7
   colnames(true_wtot) <-c("t","rep","twtot","twnf", "twnm","s")
   
   true_wtot$psi44 <-NA
   true_wtot$psi43 <-NA
   true_wtot$psi42 <-NA
   true_wtot$psi41 <-NA
   true_wtot$psi33 <-NA
   true_wtot$psi32 <-NA
   true_wtot$psi31 <-NA
   true_wtot$psi22 <-NA
   true_wtot$psi21 <-NA
   true_wtot$prob4l <-ifelse(true_wtot$twnf==1 & true_wtot$twnm==1,0,1)
   
   for (i in 1:nrow(true_wtot)){
      
      # Psi41
      true_wtot[i,10] <-(dbinom(0,true_wtot[i,3],true_wtot[i,6]))*true_state[i,6]
      
      # Psi42
      true_wtot[i,9] <-(dbinom(1,true_wtot[i,3],true_wtot[i,6]))*true_state[i,6]
      
      # Psi43
      true_wtot[i,8] <- ((pbinom(1, true_wtot[i,5], true_wtot[i,6], lower.tail = FALSE)*
                          dbinom(0, true_wtot[i,4], true_wtot[i,6])) +
                         (pbinom(1, true_wtot[i,4], true_wtot[i,6], lower.tail = FALSE)*
                          dbinom(0, true_wtot[i,5], true_wtot[i,6])))*true_state[i,6]*true_wtot[i,16]
      
      # Psi44
      true_wtot[i,7] <-(1-(true_wtot[i,8]+true_wtot[i,9]+true_wtot[i,10]))*true_state[i,6]
      
      
      # Psi31
      true_wtot[i,13] <-(dbinom(0,true_wtot[i,3],true_wtot[i,6]))*true_state[i,5]
     
      # Psi32
      true_wtot[i,12] <-(dbinom(1,true_wtot[i,3],true_wtot[i,6]))*true_state[i,5]
      
      # Psi33
      true_wtot[i,11] <-(1-(true_wtot[i,12]+true_wtot[i,13]))*true_state[i,5]
      
      
      # Psi21
      true_wtot[i,15] <-(dbinom(0,true_wtot[i,3],true_wtot[i,6]))*true_state[i,4]
      
      # Psi22
      true_wtot[i,14] <-(1-true_wtot[i,15])*true_state[i,4]
      
   }
   
   # Merge true with predicted and calculate values
   
   # PSI 44
   psi44results <-merge(true_wtot,psi44p, by = c("t","rep"), all = TRUE)
   psi44results2 <-as.data.frame(psi44results[,c(1:5,7:10,17:19)])
   colnames(psi44results2) <-c("t","rep","true_total","true_wnf","true_wnm","true_psi44","true_psi43","true_psi42", "true_psi41", "predicted_psi44", "predicted_psi44_lci","predicted_psi44_uci")
   psi44results2$coverage <-ifelse(psi44results2$true_psi44 >= psi44results2$predicted_psi44_lci & psi44results2$true_psi44 <=psi44results2$predicted_psi44_uci,1,0)
   psi44results2$mape <- ifelse(psi44results2$true_psi44==0,999,(abs(psi44results2$true_psi44 - psi44results2$predicted_psi44)/psi44results2$true_psi44)*100)
   psi44small <-psi44results2[complete.cases(psi44results2),]
   m <-mean(psi44small$true_psi44)
   psi44small$rrmse <-NA
   for (i in 1:nrow(psi44small)){
      psi44small[i,15] <-sqrt((1/nrow(psi44small))*sum((psi44small[i,6] - psi44small[i,10])^2))/m
   }
   
   # PSI 43
   psi43results <-merge(true_wtot,psi43p, by = c("t","rep"), all = TRUE)
   psi43results2 <-as.data.frame(psi43results[,c(1:5,7:10,17:19)])
   colnames(psi43results2) <-c("t","rep","true_total","true_wnf","true_wnm","true_psi44","true_psi43","true_psi42", "true_psi41", "predicted_psi43", "predicted_psi43_lci","predicted_psi43_uci")
   psi43results2$coverage <-ifelse(psi43results2$true_psi43 >= psi43results2$predicted_psi43_lci & psi43results2$true_psi43 <=psi43results2$predicted_psi43_uci,1,0)
   psi43results2$mape <- ifelse(psi43results2$true_psi43==0,999,(abs(psi43results2$true_psi43 - psi43results2$predicted_psi43)/psi43results2$true_psi43)*100)
   psi43small <-psi43results2[complete.cases(psi43results2),]
   m <-mean(psi43small$true_psi43)
   psi43small$rrmse <-NA
   for (i in 1:nrow(psi43small)){
      psi43small[i,15] <-sqrt((1/nrow(psi43small))*sum((psi43small[i,7] - psi43small[i,10])^2))/m
   }
   
   # PSI 42
   psi42results <-merge(true_wtot,psi42p, by = c("t","rep"), all = TRUE)
   psi42results2 <-as.data.frame(psi42results[,c(1:5,7:10,17:19)])
   colnames(psi42results2) <-c("t","rep","true_total","true_wnf","true_wnm","true_psi44","true_psi43","true_psi42", "true_psi41", "predicted_psi42", "predicted_psi42_lci","predicted_psi42_uci")
   psi42results2$coverage <-ifelse(psi42results2$true_psi42 >= psi42results2$predicted_psi42_lci & psi42results2$true_psi42 <=psi42results2$predicted_psi42_uci,1,0)
   psi42results2$mape <- ifelse(psi42results2$true_psi42==0,999,(abs(psi42results2$true_psi42 - psi42results2$predicted_psi42)/psi42results2$true_psi42)*100)
   psi42small <-psi42results2[complete.cases(psi42results2),]
   m <-mean(psi42small$true_psi42)
   psi42small$rrmse <-NA
   for (i in 1:nrow(psi42small)){
      psi42small[i,15] <-sqrt((1/nrow(psi42small))*sum((psi42small[i,8] - psi42small[i,10])^2))/m
   }
   
   # PSI 41
   psi41results <-merge(true_wtot,psi41p, by = c("t","rep"), all = TRUE)
   psi41results2 <-as.data.frame(psi41results[,c(1:5,7:10,17:19)])
   colnames(psi41results2) <-c("t","rep","true_total","true_wnf","true_wnm","true_psi44","true_psi43","true_psi42", "true_psi41", "predicted_psi41", "predicted_psi41_lci","predicted_psi41_uci")
   psi41results2$coverage <-ifelse(psi41results2$true_psi41 >= psi41results2$predicted_psi41_lci & psi41results2$true_psi41 <=psi41results2$predicted_psi41_uci,1,0)
   psi41results2$mape <- ifelse(psi41results2$true_psi41==0,999,(abs(psi41results2$true_psi41 - psi41results2$predicted_psi41)/psi41results2$true_psi41)*100)
   psi41small <-psi41results2[complete.cases(psi41results2),]
   m <-mean(psi41small$true_psi41)
   psi41small$rrmse <-NA
   for (i in 1:nrow(psi41small)){
      psi41small[i,15] <-sqrt((1/nrow(psi41small))*sum((psi41small[i,9] - psi41small[i,10])^2))/m
   }

 # # PSI 33
 # psi33results <-merge(true_wtot,psi33p, by = c("t","rep"), all = TRUE)
 # psi33results2 <-as.data.frame(psi33results[,c(1:5,11:13,17:19)])
 # colnames(psi33results2) <-c("t","rep","true_total","true_wnf","true_wnm","true_psi33","true_psi32", "true_psi31", "predicted_psi33","predicted_psi33_lci","predicted_psi33_uci")
 # psi33results2$coverage <-ifelse(psi33results2$true_psi33 >= psi33results2$predicted_psi33_lci & psi33results2$true_psi33 <=psi33results2$predicted_psi33_uci,1,0)
 # psi33results2$mape <- ifelse(psi33results2$true_psi33==0,999,(abs(psi33results2$true_psi33 - psi33results2$predicted_psi33)/psi33results2$true_psi33)*100)
 # psi33small <-psi33results2[complete.cases(psi33results2),]
 # m <-mean(psi33small$true_psi33)
 # psi33small$rrmse <-NA
 # for (i in 1:nrow(psi33small)){
 #    psi33small[i,14] <-sqrt((1/nrow(psi33small))*sum((psi33small[i,6] - psi33small[i,9])^2))/m
 # }
 # 
 # # PSI 32
 # psi32results <-merge(true_wtot,psi32p, by = c("t","rep"), all = TRUE)
 # psi32results2 <-as.data.frame(psi32results[,c(1:5,11:13,17:19)])
 # colnames(psi32results2) <-c("t","rep","true_total","true_wnf","true_wnm","true_psi33","true_psi32", "true_psi31", "predicted_psi32","predicted_psi32_lci","predicted_psi32_uci")
 # psi32results2$coverage <-ifelse(psi32results2$true_psi32 >= psi32results2$predicted_psi32_lci & psi32results2$true_psi32 <=psi32results2$predicted_psi32_uci,1,0)
 # psi32results2$mape <- ifelse(psi32results2$true_psi32==0,999,(abs(psi32results2$true_psi32 - psi32results2$predicted_psi32)/psi32results2$true_psi32)*100)
 # psi32small <-psi32results2[complete.cases(psi32results2),]
 # m <-mean(psi32small$true_psi32)
 # psi32small$rrmse <-NA
 # for (i in 1:nrow(psi32small)){
 #    psi32small[i,14] <-sqrt((1/nrow(psi32small))*sum((psi32small[i,7] - psi32small[i,9])^2))/m
 # }
 # 
 # # PSI 31
 # psi31results <-merge(true_wtot,psi31p, by = c("t","rep"), all = TRUE)
 # psi31results2 <-as.data.frame(psi31results[,c(1:5,11:13,17:19)])
 # colnames(psi31results2) <-c("t","rep","true_total","true_wnf","true_wnm","true_psi33","true_psi32", "true_psi31", "predicted_psi31","predicted_psi31_lci","predicted_psi31_uci")
 # psi31results2$coverage <-ifelse(psi31results2$true_psi31 >= psi31results2$predicted_psi31_lci & psi31results2$true_psi31 <=psi31results2$predicted_psi31_uci,1,0)
 # psi31results2$mape <- ifelse(psi31results2$true_psi31==0,999,(abs(psi31results2$true_psi31 - psi31results2$predicted_psi31)/psi31results2$true_psi31)*100)
 # psi31small <-psi31results2[complete.cases(psi31results2),]
 # m <-mean(psi31small$true_psi31)
 # psi31small$rrmse <-NA
 # for (i in 1:nrow(psi32small)){
 #    psi31small[i,14] <-sqrt((1/nrow(psi31small))*sum((psi31small[i,8] - psi31small[i,9])^2))/m
 # }
 # 
 # # PSI 22
 # psi22results <-merge(true_wtot,psi22p, by = c("t","rep"), all = TRUE)
 # psi22results2 <-as.data.frame(psi22results[,c(1:5,14:15,17:19)])
 # colnames(psi22results2) <-c("t","rep","true_total","true_wnf","true_wnm", "true_psi22", "true_psi21","predicted_psi22","predicted_psi22_lci","predicted_psi22_uci")
 # psi22results2$coverage <-ifelse(psi22results2$true_psi22 >= psi22results2$predicted_psi22_lci & psi22results2$true_psi22 <=psi22results2$predicted_psi22_uci,1,0)
 # psi22results2$mape <- ifelse(psi22results2$true_psi22==0,999,(abs(psi22results2$true_psi22 - psi22results2$predicted_psi22)/psi22results2$true_psi22)*100)
 # psi22small <-psi22results2[complete.cases(psi22results2),]
 # m <-mean(psi22small$true_psi22)
 # psi22small$rrmse <-NA
 # for (i in 1:nrow(psi22small)){
 #    psi22small[i,13] <-sqrt((1/nrow(psi22small))*sum((psi22small[i,6] - psi22small[i,8])^2))/m
 # }
 # 
 # # PSI 21
 # psi21results <-merge(true_wtot,psi21p, by = c("t","rep"), all = TRUE)
 # psi21results2 <-as.data.frame(psi21results[,c(1:5,14:15,17:19)])
 # colnames(psi21results2) <-c("t","rep","true_total","true_wnf","true_wnm", "true_psi22", "true_psi21","predicted_psi21","predicted_psi21_lci","predicted_psi21_uci")
 # psi21results2$coverage <-ifelse(psi21results2$true_psi21 >= psi21results2$predicted_psi21_lci & psi21results2$true_psi21 <=psi21results2$predicted_psi21_uci,1,0)
 # psi21results2$mape <- ifelse(psi21results2$true_psi21==0,999,(abs(psi21results2$true_psi21 - psi21results2$predicted_psi21)/psi21results2$true_psi21)*100)
 # psi21small <-psi21results2[complete.cases(psi21results2),]
 # m <-mean(psi21small$true_psi21)
 # psi21small$rrmse <-NA
 # for (i in 1:nrow(psi21small)){
 #    psi21small[i,13] <-sqrt((1/nrow(psi21small))*sum((psi21small[i,7] - psi21small[i,8])^2))/m
 # }
      
####################### Write data to global environment #######################
txcalcdata <- list("psi44_results" = psi44small, "psi43_results" = psi43small, 
                   "psi42_results" = psi42small, "psi41_results" = psi41small) 
                   # "psi33_results" = psi33small, "psi32_results" = psi32small,
                   # "psi31_results" = psi31small, "psi22_results" = psi22small,
                   # "psi21_results" = psi21small)

      list2env(txcalcdata,.GlobalEnv)
}
