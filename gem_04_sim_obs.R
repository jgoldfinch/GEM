# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_sim_obs
# Description:  function to observe biological conditions of one or multiple
#               populations with four states (not present, single 
#               individual present, multiple individuals of a single sex present,
#               and multiple individuals with both sexes present)

################################# Arguments ####################################

# p:
#       Detection probability. Number between 0 and 1.
# nf: 
#       Adult females from population which you observe. Array generated with
#       gem_sim_bio function.
# nm: 
#       Adult males from population which you observe. Array generated with
#       gem_sim_bio function.
# n.visits: 
#       Number of visits to each group to survey.
# n.group:
#       Number of populations (groups). Whole number.
# n.timestep: 
#       Number of time steps you want to observe. Whole number.
# n.states: 
#       Number of GEM population states. Whole number.
# z: 
#       Occupancy of populations. Array generated with gem_sim_bio function.
# pgenetic: 
#       Detection probability of genetic sign. Number between 0 and 1.

################################## Output ######################################

# obsdata list saved to global environment. Contains p (defined above), y (observed
# count dataframe), n.obs (number of observations), group (group index for all
# observations), timestep (time index for all observations)

################################# Function #####################################

gem_sim_obs <- function(p, nf, nm, n.visits, n.group, n.timestep, n.states, z, pgenetic){
  
### Create matrices/arrays to hold data
  # Present, counts, genetic sign, females genetic sign, males genetic sign
  yp <- yc <- yg <- ygf <- ygm <-array(data = NA, dim = c(n.visits, n.timestep, n.group))

for (k in 1:n.visits){ 
  for (i in 1:n.group){         
    for(j in 2:n.timestep){ 
    
    ## Time 1
    ## Is the species present?
      yp[k,1,i] <- rbinom(1,ifelse(z[1,1,i]>1,1,0),1-((1-p)^(nf[1,1,i]+nm[1,1,i])))
      
    ## Time 2 and beyond
    ## Is the species present?
      yp[k,j,i] <- rbinom(1,ifelse(z[1,j,i]>1,1,0),1-((1-p)^(nf[1,j,i]+nm[1,j,i])))
      
    ## Are multiple individuals present?
      yc[k,j,i] <- ifelse(any(yp[,j-1,i]==1),rbinom(1,nf[1,j,i]+nm[1,j,i],p),NA)
    
    ## Are females and males present?
      yg[k,j,i] <- ifelse(any(yc[,j-1,i]>1),rbinom(1,yc[1,j,i],pgenetic),NA)
      ygf[k,j,i] <-rhyper(1,nf[1,j,i],nm[1,j,i],yg[k,j,i])
      ygm[k,j,i] <-yg[k,j,i]-ygf[k,j,i]
    }   
  }
}

################################################################################
  obsdata <- list("p"= p, "yp"= yp,"yc"= yc, "yg" = yg, "ygf"=ygf, "ygm" = ygm,
                  "pgenetic" = pgenetic, "n.visits" = n.visits)
  list2env(obsdata ,.GlobalEnv)
}

