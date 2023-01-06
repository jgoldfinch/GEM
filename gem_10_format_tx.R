# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_format_tx
# Description:  function to format data needed to calculate predicted GEM state 
#               transition probabilities from GEM model output. GEM state 
#               transitions are changes in populations with four states 
#               (not present, single individual present, multiple individuals of
#               a single sex present, and multiple individuals with both sexes 
#               present)

################################# Arguments ####################################

# n.group:
#       Number of populations (groups). Whole number.
# s.group: 
#       Starting number of females and males in each population. Whole number.
# n.timestep: 
#       Number of time steps you want to observe. Whole number.
# n.states: 
#       Number of GEM population states. Whole number.
# s.surv: 
#       Adult survival. Number between 0 and 1. 
# p.litter: 
#       Probability of having a litter. Number between 0 and 1.
# n.litter: 
#       Number of individuals per litter. Number between 0 and 1.
# sr.litter: 
#       Sex ratio of males to females per litter. Number between 0 and 1.
# p:
#       Detection probability. Number between 0 and 1.
# n.visits: 
#       Number of visits to each group to survey.
# n.reps: 
#       Number of simulation replicates.

################################## Output ######################################

# txpreddata: list saved to global environment. Contains formatted model output 
# data necessary to calculate transitions. This includes: 

# lists of the entire model output distributions for total predicted 
# individuals at the end of each time step (wtot), total predicted females 
# (new and existing adult) at the end of each time step (wnm), total predicted 
# females (new and existing adult) at the end of each time step (wnf)

# summarized record of when each simulation was in each population state 
# (prob4 = 1 if population was in state 4, if not = 0,
# prob3 = 1 if population was in state 3, if not = 0, prob2 = 1 if population was
# in state 2, if not = 0)

# lists of the true values for total individuals at the end of each 
# time step (twtot), total females (new and existing adult) at the end of each
# time step (twnm), total females (new and existing adult) at the end of each 
# time step (twnf), occupancy of the population (tprob)

# lists of the entire model output distributions for predicted survival (s_results)

################################# Function #####################################

gem_format_tx <- function(n.group, s.group, n.timestep, n.states, s.surv, p.litter, 
                          n.litter, sr.litter, p, n.visits, n.reps){
   
   # Create empty lists
   wtot <-lapply(wtot<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   wnf <-lapply(wnf<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   wnm <-lapply(wnm<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   prob1<-lapply(prob1<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   prob2<-lapply(prob2<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   prob3<-lapply(prob3<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   prob4<-lapply(prob4<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   
   twtot <-lapply(twtot<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   twnf <-lapply(twnf<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   twnm <-lapply(twnm<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   tprob<-lapply(tprob1<-vector(mode = 'list',n.reps),function(x) x<-vector(mode='list',n.timestep))
   s_results <- vector("list", n.reps)
   
for (k in 1:n.reps){
      name <-paste("rep",k,"_", n.timestep, "t", "_", n.group,"g", "_", s.group,"st", "_", 
                   n.states, "states", "_", s.surv, "surv", "_", p.litter, "plitter",
                   "_", n.litter, "perlitter","_", sr.litter, "sr","_", p,"p", "_", 
                   n.visits, "visits", sep="")
      
      filename = paste0 ("./", name)
      load(filename)
      
      s_results[[k]] <- finalout[[1]]$sims.list$s.surv
      
for(i in 1:n.timestep){  
   
      wtot[[k]][[i]] <- finalout[[1]]$sims.list$wtot[,1,i,1]
      wnf[[k]][[i]] <- finalout[[1]]$sims.list$wnf[,1,i,1]
      wnm[[k]][[i]] <- finalout[[1]]$sims.list$wnm[,1,i,1]
      prob1[[k]][[i]] <-ifelse(round(finalout[[1]]$mean$z[1,i,1])==1,1,0)
      prob2[[k]][[i]] <-ifelse(round(finalout[[1]]$mean$z[1,i,1])==2,1,0)
      prob3[[k]][[i]] <-ifelse(round(finalout[[1]]$mean$z[1,i,1])==3,1,0)
      prob4[[k]][[i]] <-ifelse(round(finalout[[1]]$mean$z[1,i,1])==4,1,0)
      
      twtot[[k]][[i]] <- finalout[[2]][finalout[[2]]$variable =="wtot" & finalout[[2]]$time ==i,4]
      twnf[[k]][[i]] <- finalout[[2]][finalout[[2]]$variable =="wnf" & finalout[[2]]$time ==i,4]
      twnm[[k]][[i]] <- finalout[[2]][finalout[[2]]$variable =="wnm" & finalout[[2]]$time ==i,4]
      tprob[[k]][[i]] <-finalout[[2]][finalout[[2]]$variable =="z" & finalout[[2]]$time ==i,4]
}
}
      
####################### Write data to global environment #######################
txpreddata <- list("wtot" = wtot, "wnm" = wnm, "wnf" = wnf, "prob1" = prob1,
                   "prob2" = prob2, "prob3" = prob3, "prob4" = prob4, 
                   "twnf" = twnf, "twnm" = twnm,"twtot" = twtot, "tprob" = tprob, 
                   "s_results" = s_results)

      list2env(txpreddata,.GlobalEnv)
}
