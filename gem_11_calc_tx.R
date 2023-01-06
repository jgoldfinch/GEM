# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Name: gem_calc_tx
# Description:  function to calculate predicted GEM transition probabilities 
#               from model, which are changes in populations with four states 
#               (not present, single individual present, multiple individuals 
#               of a single sex present, and multiple individuals with both 
#               sexes present)

################################# Arguments ####################################

# prob2:
#       List of whether the population was in GEM state 2 (single individual) during 
#       each time step (prob2 = 1 if population was in state 2, if not = 0) 
#       calculated from gem_format_tx function.
# prob3:
#       List of whether the population was in GEM state 3 (multiple individuals, 
#       single sex) during each time step (prob3 = 1 if population was in state 
#       3, if not = 0) calculated from gem_format_tx function.
# prob4:
#       List of whether the population was in GEM state 3 (multiple individuals, 
#       both sexes) during each time step (prob4 = 1 if population was in state 
#       4, if not = 0) calculated from gem_format_tx function.
# wtot: Formatted total individuals (adult and new) alive at the end of a time 
#       step predictions.
# wnf:  Formatted total females (adult and new) alive at the end of a time step 
#       predictions.
# wnm:  Formatted total males (adult and new) alive at the end of a time step 
#       predictions.
# s_results: 
#       Survival predictions produced with gem_run_model function.
# n.timestep: 
#       Number of time steps you want to observe. Whole number.
# n.reps: 
#       Number of simulation replicates.

################################## Output ######################################

# txcalcdata: list saved to global environment. Contains predicted transition 
# probabilities for state 4 (psi44p, psi43p, psi42p, psi41p), state 3 (psi33p,
# psi32p, psi31p), and state 2 (psi22p, psi21p).

################################# Function #####################################

gem_calc_tx <- function(prob2, prob3, prob4, wtot, wnf, wnm, s_results, 
                        n.timestep, n.reps){
   
   # Create empty lists
   psi44 <-lapply(psi44<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi44_mean <-lapply(psi44_mean<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi44_lowci <-lapply(psi44_lowci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi44_highci <-lapply(psi44_highci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
   psi43 <-lapply(psi43<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi43_mean <-lapply(psi43_mean<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi43_lowci <-lapply(psi43_lowci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi43_highci <-lapply(psi43_highci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
   psi42 <-lapply(psi42<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi42_mean <-lapply(psi42_mean<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi42_lowci <-lapply(psi42_lowci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi42_highci <-lapply(psi42_highci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
   psi41 <-lapply(psi41<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi41_mean <-lapply(psi41_mean<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi41_lowci <-lapply(psi41_lowci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi41_highci <-lapply(psi41_highci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
   psi33 <-lapply(psi33<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi33_mean <-lapply(psi33_mean<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi33_lowci <-lapply(psi33_lowci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi33_highci <-lapply(psi33_highci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
   psi32 <-lapply(psi32<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi32_mean <-lapply(psi32_mean<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi32_lowci <-lapply(psi32_lowci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi32_highci <-lapply(psi32_highci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
   psi31 <-lapply(psi31<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi31_mean <-lapply(psi31_mean<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi31_lowci <-lapply(psi31_lowci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi31_highci <-lapply(psi31_highci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
   psi22 <-lapply(psi22<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi22_mean <-lapply(psi22_mean<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi22_lowci <-lapply(psi22_lowci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi22_highci <-lapply(psi22_highci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
   psi21 <-lapply(psi21<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi21_mean <-lapply(psi21_mean<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi21_lowci <-lapply(psi21_lowci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   psi21_highci <-lapply(psi21_highci<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
for (k in 1:n.reps){
for (i in 1:n.timestep){  
   
      # State 4 transitions
      psi41[[k]][[i]] <-dbinom(0,wtot[[k]][[i]],s_results[[k]])
      psi42[[k]][[i]] <-dbinom(1,wtot[[k]][[i]],s_results[[k]])
      psi43[[k]][[i]] <-(pbinom(1, wnm[[k]][[i]], s_results[[k]], lower.tail = FALSE)*dbinom(0,wnf[[k]][[i]], s_results[[k]]))+
                        (pbinom(1, wnf[[k]][[i]], s_results[[k]], lower.tail = FALSE)*dbinom(0,wnm[[k]][[i]], s_results[[k]]))
      psi44[[k]][[i]] <-1-(psi41[[k]][[i]]+psi42[[k]][[i]]+psi43[[k]][[i]])
      
      
      # psi44 estimate for each time step
      psi44_lowci[[k]][[i]] <-(quantile(psi44[[k]][[i]], probs = c(0.025)))*prob4[[k]][[i]]
      psi44_mean[[k]][[i]] <-(mean(psi44[[k]][[i]]))*prob4[[k]][[i]]
      psi44_highci[[k]][[i]] <-(quantile(psi44[[k]][[i]], probs = c(0.975)))*prob4[[k]][[i]]
      
      # psi43 estimate for each time step
      psi43_lowci[[k]][[i]] <-(quantile(psi43[[k]][[i]], probs = c(0.025)))*prob4[[k]][[i]]
      psi43_mean[[k]][[i]] <-(mean(psi43[[k]][[i]]))*prob4[[k]][[i]]
      psi43_highci[[k]][[i]] <-(quantile(psi43[[k]][[i]], probs = c(0.975)))*prob4[[k]][[i]]
      
      # psi42 estimate for each time step
      psi42_lowci[[k]][[i]] <-(quantile(psi42[[k]][[i]], probs = c(0.025)))*prob4[[k]][[i]]
      psi42_mean[[k]][[i]] <- (mean(psi42[[k]][[i]]))*prob4[[k]][[i]]
      psi42_highci[[k]][[i]] <-(quantile(psi42[[k]][[i]], probs = c(0.975)))*prob4[[k]][[i]]
      
      # psi41 estimate for each time step
      psi41_lowci[[k]][[i]] <-(quantile(psi41[[k]][[i]], probs = c(0.025)))*prob4[[k]][[i]]
      psi41_mean[[k]][[i]] <- (mean(psi41[[k]][[i]]))*prob4[[k]][[i]]
      psi41_highci[[k]][[i]] <-(quantile(psi41[[k]][[i]], probs = c(0.975)))*prob4[[k]][[i]]
      
      # State 3 transitions
      psi31[[k]][[i]] <-dbinom(0,wtot[[k]][[i]],s_results[[k]])
      psi32[[k]][[i]] <-dbinom(1,wtot[[k]][[i]],s_results[[k]])
      psi33[[k]][[i]] <-1-(psi31[[k]][[i]]+psi32[[k]][[i]])

      # psi33 estimate for each time step
      psi33_lowci[[k]][[i]] <-(quantile(psi33[[k]][[i]], probs = c(0.025)))*prob3[[k]][[i]]
      psi33_mean[[k]][[i]] <-(mean(psi33[[k]][[i]]))*prob3[[k]][[i]]
      psi33_highci[[k]][[i]] <-(quantile(psi33[[k]][[i]], probs = c(0.975)))*prob3[[k]][[i]]
   
      # psi32 estimate for each time step
      psi32_lowci[[k]][[i]] <-(quantile(psi32[[k]][[i]], probs = c(0.025)))*prob3[[k]][[i]]
      psi32_mean[[k]][[i]] <- (mean(psi32[[k]][[i]]))*prob3[[k]][[i]]
      psi32_highci[[k]][[i]] <-(quantile(psi32[[k]][[i]], probs = c(0.975)))*prob3[[k]][[i]]
      
      # psi31 estimate for each time step
      psi31_lowci[[k]][[i]] <-(quantile(psi31[[k]][[i]], probs = c(0.025)))*prob3[[k]][[i]]
      psi31_mean[[k]][[i]] <- (mean(psi31[[k]][[i]]))*prob3[[k]][[i]]
      psi31_highci[[k]][[i]] <-(quantile(psi31[[k]][[i]], probs = c(0.975)))*prob3[[k]][[i]]
      
      # State 2 transitions
      psi21[[k]][[i]] <-dbinom(0,wtot[[k]][[i]],s_results[[k]])
      psi22[[k]][[i]] <-1-psi21[[k]][[i]]
      
      # psi22 estimate for each time step
      psi22_lowci[[k]][[i]] <-(quantile(psi22[[k]][[i]], probs = c(0.025)))*prob2[[k]][[i]]
      psi22_mean[[k]][[i]] <-(mean(psi22[[k]][[i]]))*prob2[[k]][[i]]
      psi22_highci[[k]][[i]] <-(quantile(psi22[[k]][[i]], probs = c(0.975)))*prob2[[k]][[i]]
      
      # psi21 estimate for each time step
      psi21_lowci[[k]][[i]] <-(quantile(psi21[[k]][[i]], probs = c(0.025)))*prob2[[k]][[i]]
      psi21_mean[[k]][[i]] <-(mean(psi21[[k]][[i]]))*prob2[[k]][[i]]
      psi21_highci[[k]][[i]] <-(quantile(psi21[[k]][[i]], probs = c(0.975)))*prob2[[k]][[i]]
      
}
}
   
# Summarize and simplify formats
   
   # State 4
   # Stay in 4
      psi44_means <-melt(psi44_mean)
      colnames(psi44_means) <-c("psi44_mean","t","rep")
      
      psi44_highcis <-melt(psi44_highci)
      colnames(psi44_highcis) <-c("psi44_highci","t","rep")
      
      psi44_lowcis <-melt(psi44_lowci)
      colnames(psi44_lowcis) <-c("psi44_lowci","t","rep")
      
      psi44_all <-cbind(psi44_means$t, psi44_means$rep, psi44_means$psi44_mean, psi44_lowcis$psi44_lowci, psi44_highcis$psi44_highci)
      colnames(psi44_all) <-c("t","rep","mean", "lowci", "uci")
      psi44_all <-psi44_all[psi44_all[,3] != 0, ]
      
      # Transition to single sex, multiple individuals
      psi43_means <-melt(psi43_mean)
      colnames(psi43_means) <-c("psi43_mean","t","rep")
      
      psi43_highcis <-melt(psi43_highci)
      colnames(psi43_highcis) <-c("psi43_highci","t","rep")
      
      psi43_lowcis <-melt(psi43_lowci)
      colnames(psi43_lowcis) <-c("psi43_lowci","t","rep")
      
      psi43_all <-cbind(psi43_means$t, psi43_means$rep, psi43_means$psi43_mean, psi43_lowcis$psi43_lowci, psi43_highcis$psi43_highci)
      colnames(psi43_all) <-c("t","rep","mean", "lowci", "uci")
      psi43_all <-psi43_all[psi43_all[,3] != 0, ]
      
      # Transition to single individual
      psi42_means <-melt(psi42_mean)
      colnames(psi42_means) <-c("psi42_mean","t","rep")
      
      psi42_highcis <-melt(psi42_highci)
      colnames(psi42_highcis) <-c("psi42_highci","t","rep")
      
      psi42_lowcis <-melt(psi42_lowci)
      colnames(psi42_lowcis) <-c("psi42_lowci","t","rep")
      
      psi42_all <-cbind(psi42_means$t, psi42_means$rep, psi42_means$psi42_mean, psi42_lowcis$psi42_lowci, psi42_highcis$psi42_highci)
      colnames(psi42_all) <-c("t","rep","mean", "lowci", "uci")
      psi42_all <-psi42_all[psi42_all[,3] != 0, ]
      
      # Transition to not present
      psi41_means <-melt(psi41_mean)
      colnames(psi41_means) <-c("psi41_mean","t","rep")
      
      psi41_highcis <-melt(psi41_highci)
      colnames(psi41_highcis) <-c("psi41_highci","t","rep")
      
      psi41_lowcis <-melt(psi41_lowci)
      colnames(psi41_lowcis) <-c("psi41_lowci","t","rep")
      
      psi41_all <-cbind(psi41_means$t, psi41_means$rep, psi41_means$psi41_mean, psi41_lowcis$psi41_lowci, psi41_highcis$psi41_highci)
      colnames(psi41_all) <-c("t","rep","mean", "lowci", "uci")
      psi41_all <-psi41_all[psi41_all[,3] != 0, ]
      
   # State 3
      # Stay in 3
      psi33_means <-melt(psi33_mean)
      colnames(psi33_means) <-c("psi33_mean","t","rep")
     
      psi33_highcis <-melt(psi33_highci)
      colnames(psi33_highcis) <-c("psi33_highci","t","rep")
      
      psi33_lowcis <-melt(psi33_lowci)
      colnames(psi33_lowcis) <-c("psi33_lowci","t","rep")
      
      psi33_all <-cbind(psi33_means$t, psi33_means$rep, psi33_means$psi33_mean, psi33_lowcis$psi33_lowci, psi33_highcis$psi33_highci)
      colnames(psi33_all) <-c("t","rep","mean", "lowci", "uci")
      psi33_all <-psi33_all[psi33_all[,3] != 0, ]
      
      # Transition to single individual
      psi32_means <-melt(psi32_mean)
      colnames(psi32_means) <-c("psi32_mean","t","rep")
      
      psi32_highcis <-melt(psi32_highci)
      colnames(psi32_highcis) <-c("psi32_highci","t","rep")
      
      psi32_lowcis <-melt(psi32_lowci)
      colnames(psi32_lowcis) <-c("psi32_lowci","t","rep")
      
      psi32_all <-cbind(psi32_means$t, psi32_means$rep, psi32_means$psi32_mean, psi32_lowcis$psi32_lowci, psi32_highcis$psi32_highci)
      colnames(psi32_all) <-c("t","rep","mean", "lowci", "uci")
      psi32_all <-psi32_all[psi32_all[,3] != 0, ]
      
      # Transition to not present
      psi31_means <-melt(psi31_mean)
      colnames(psi31_means) <-c("psi31_mean","t","rep")
      
      psi31_highcis <-melt(psi31_highci)
      colnames(psi31_highcis) <-c("psi31_highci","t","rep")
      
      psi31_lowcis <-melt(psi31_lowci)
      colnames(psi31_lowcis) <-c("psi31_lowci","t","rep")
      
      psi31_all <-cbind(psi31_means$t, psi31_means$rep, psi31_means$psi31_mean, psi31_lowcis$psi31_lowci, psi31_highcis$psi31_highci)
      colnames(psi31_all) <-c("t","rep","mean", "lowci", "uci")
      psi31_all <-psi31_all[psi31_all[,3] != 0, ]
      
      
      # State 2
      # Stay in 2
      psi22_means <-melt(psi22_mean)
      colnames(psi22_means) <-c("psi22_mean","t","rep")
      
      psi22_highcis <-melt(psi22_highci)
      colnames(psi22_highcis) <-c("psi22_highci","t","rep")
      
      psi22_lowcis <-melt(psi22_lowci)
      colnames(psi22_lowcis) <-c("psi22_lowci","t","rep")
      
      psi22_all <-cbind(psi22_means$t, psi22_means$rep, psi22_means$psi22_mean, psi22_lowcis$psi22_lowci, psi22_highcis$psi22_highci)
      colnames(psi22_all) <-c("t","rep","mean", "lowci", "uci")
      psi22_all <-psi22_all[psi22_all[,3] != 0, ]
      
      # Transition to not present
      psi21_means <-melt(psi21_mean)
      colnames(psi21_means) <-c("psi21_mean","t","rep")
      
      psi21_highcis <-melt(psi21_highci)
      colnames(psi21_highcis) <-c("psi21_highci","t","rep")
      
      psi21_lowcis <-melt(psi21_lowci)
      colnames(psi21_lowcis) <-c("psi21_lowci","t","rep")
      
      psi21_all <-cbind(psi21_means$t, psi21_means$rep, psi21_means$psi21_mean, psi21_lowcis$psi21_lowci, psi21_highcis$psi21_highci)
      colnames(psi21_all) <-c("t","rep","mean", "lowci", "uci")
      psi21_all <-psi21_all[psi21_all[,3] != 0, ]
      
####################### Write data to global environment #######################
txcalcdata <- list("psi44p" = psi44_all, "psi43p" = psi43_all, 
                   "psi42p" = psi42_all, "psi41p" = psi41_all, 
                   "psi33p" = psi33_all, "psi32p" = psi32_all,
                   "psi31p" = psi31_all, "psi22p" = psi22_all,
                   "psi21p" = psi21_all)

      list2env(txcalcdata,.GlobalEnv)
}
