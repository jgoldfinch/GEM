# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_sim_bio
# Description:  function to simulate biological conditions of one or multiple
#               small populations with four states (not present, single 
#               individual present, multiple individuals of a single sex present,
#               and multiple individuals with both sexes present)

################################# Arguments ####################################

# n.group:
#       Number of populations (groups). Whole number.
# s.group: 
#       Starting number of females and males in each population. Whole number.
# n.timestep: 
#       Number of timesteps. Whole number.
# n.states: 
#       Number of states. Whole number.
# s.surv: 
#       Adult survival. Number between 0 and 1.
# p.litter: 
#       Probability of having a litter. Number between 0 and 1.
# n.litter: 
#       Number of individuals per litter. Number between 0 and 1.
# sr.litter: 
#       Sex ratio of males to females per litter. Number between 0 and 1.

################################## Output ######################################

# biodata list saved to global environment. Contains n.group, s.group, n.timestep,
# n.states, s.surv, p.litter, n.litter, sr.litter (definitions above), biological groups,
# N_long (for population graphing)

################################# Function #####################################

gem_sim_bio <- function(n.group, s.group, n.timestep,n.states, s.surv, p.litter, 
                        n.litter, sr.litter){
### Create matrices/arrays to hold data
# Start abundance
nm <-array(data = NA, dim = c(1,n.timestep,n.group),
          dimnames = list(c("nm"),
                          c(1:n.timestep),
                          c(1:n.group)))
nf <-array(data = NA, dim = c(1,n.timestep,n.group),
           dimnames = list(c("nf"),
                           c(1:n.timestep),
                           c(1:n.group)))

# Birth events
be <-array(data = NA, dim = c(1,n.timestep,n.group),
           dimnames = list(c("be"),
                           c(1:n.timestep),
                           c(1:n.group)))

# New individuals
ni <-array(data = NA, dim = c(1,n.timestep,n.group),
           dimnames = list(c("ni"),
                           c(1:n.timestep),
                           c(1:n.group)))

# New males
wm <-array(data = NA, dim = c(1,n.timestep,n.group),
           dimnames = list(c("wm"),
                           c(1:n.timestep),
                           c(1:n.group)))

# New females
wf <-array(data = NA, dim = c(1,n.timestep,n.group),
           dimnames = list(c("wf"),
                           c(1:n.timestep),
                           c(1:n.group)))

# End abundance
wnm <-array(data = NA, dim = c(1,n.timestep,n.group),
           dimnames = list(c("wnm"),
                           c(1:n.timestep),
                           c(1:n.group)))
wnf <-array(data = NA, dim = c(1,n.timestep,n.group),
           dimnames = list(c("wnf"),
                           c(1:n.timestep),
                           c(1:n.group)))

wtot <-array(data = NA, dim = c(1,n.timestep,n.group),
            dimnames = list(c("wtot"),
                            c(1:n.timestep),
                            c(1:n.group)))

# Lambda
lam <-array(data = NA, dim = c(1,n.timestep,n.group),
          dimnames = list(c("l"),
                          c(1:n.timestep),
                          c(1:n.group)))

# Z - states
z <-array(data = NA, dim = c(1,n.timestep,n.group),
            dimnames = list(c("z"),
                            c(1:n.timestep),
                            c(1:n.group)))

nfa1 <-NULL
nma1 <-NULL
  
# Reproduction
litter_prob<-matrix(NA,n.timestep,n.group)
state_breed <-matrix(NA,n.timestep,n.group)

# Survival
S <-matrix(NA,n.group,n.timestep)

# Probability of litter
L <-matrix(NA,n.group,n.timestep) 
  
# Fill in survival and probability of litter
for (l in 1:n.group){ 
for (k in 1:n.timestep){ 
  S[l,k] <-s.surv 
  L[l,k] <-p.litter 
}
}

############################ Initial time 1 values ############################# 
## Initial states for time 1
# Initial population size at time 1
for (i in 1:n.group){ 
  nfa1[i] <- rpois(1,s.group) # females
  nma1[i] <- rpois(1,s.group) # males
}
######################### Loop to generate population ########################## 
# Population
for (i in 1:n.group){         
  for(j in 2:n.timestep){     
## Time 1
# Entering individuals/state
  nm[1,1,i] <- as.numeric(nma1[i]) 
  nf[1,1,i] <- as.numeric(nfa1[i]) 
  z[1,1,i] <-ifelse(nm[1,1,i] >= 1 & nf[1,1,i] >=1, 4,
             ifelse(nm[1,1,i] >  1 & nf[1,1,i] == 0|nf[1,1,i] > 1 & nm[1,1,i] == 0, 3,
             ifelse(nm[1,1,i] == 1 & nf[1,1,i] == 0|nf[1,1,i] ==1 & nm[1,1,i] == 0, 2,
             ifelse(nm[1,1,i] == 0 & nf[2,1,i] == 0,1,NA))))
  
# Breeding possible
  state_breed[1,i] <-as.numeric(ifelse(z[1,1,i]==4,1,0))
  litter_prob[1,i] <-state_breed[1,i]*L[i,1] 
  
# Births
  be[1,1,i] <- rbinom(1,nf[1,1,i],litter_prob[1,i]) 

# New individuals 
  ni[1,1,i] <- be[1,1,i]*n.litter
  
  wm[1,1,i] <- rbinom(1,ni[1,1,i], sr.litter) # new males
  wf[1,1,i] <- ni[1,1,i] - wm[1,1,i] # new females
  
# Totals
  wnm[1,1,i] <- wm[1,1,i] + nm[1,1,i] # total males at time 1
  wnf[1,1,i] <- wf[1,1,i] + nf[1,1,i] # total females at time 1
  wtot[1,1,i] <-wnm[1,1,i] + wnf[1,1,i] # total at time 1
  
# Lambda
  lam[1,1,i] <- NA

## Time 2 and beyond
## Entering individuals/state
  nm[1,j,i] <- rbinom(1,wnm[1,j-1,i],S[i,j])
  nf[1,j,i] <- rbinom(1,wnf[1,j-1,i],S[i,j])
  z[1,j,i] <-ifelse(nm[1,j,i] >= 1 & nf[1,j,i] >=1, 4,
             ifelse(nm[1,j,i] > 1  & nf[1,j,i] == 0|nf[1,j,i] > 1 & nm[1,j,i] == 0, 3,
             ifelse(nm[1,j,i] == 1 & nf[1,j,i] == 0|nf[1,j,i] ==1 & nm[1,j,i] == 0, 2,
             ifelse(nm[1,j,i] == 0 & nf[1,j,i] == 0,1,NA))))

# Breeding possible
  state_breed[j,i] <-as.numeric(ifelse(z[1,j,i]==4,1,0))
  litter_prob[j,i] <-state_breed[j,i]*L[i,j] 
  
# Births
  be[1,j,i] <- rbinom(1, nf[1,j,i], litter_prob[j,i])
  
# New individuals 
  ni[1,j,i] <- be[1,j,i]*n.litter
  
  wm[1,j,i] <- rbinom(1,ni[1,j,i],sr.litter) # new males
  wf[1,j,i] <- ni[1,j,i] - wm[1,j,i] # new females
  
# Totals
  wnm[1,j,i] <- wm[1,j,i] + nm[1,j,i] # total males
  wnf[1,j,i] <- wf[1,j,i] + nf[1,j,i] # total females 
  wtot[1,j,i] <- wnm[1,j,i] + wnf[1,j,i] # total

  # Lambda
  lam[1,j,i] <- wtot[1,j,i]/wtot[1,j-1,i]

}
}

################################################################################
### Create long data frames for plotting and result/diagnostic comparisons

  nm2 <-as.data.frame(melt(nm))
  colnames(nm2) <-c("variable","time","group","value")
  
  nf2 <-as.data.frame(melt(nf))
  colnames(nf2) <-c("variable","time","group","value")
  
  z2 <-as.data.frame(melt(z))
  colnames(z2) <-c("variable","time","group","value")
  
  be2 <-as.data.frame(melt(be))
  colnames(be2) <-c("variable","time","group","value")
  
  ni2 <-as.data.frame(melt(ni))
  colnames(ni2) <-c("variable","time","group","value")
  
  wm2 <-as.data.frame(melt(wm))
  colnames(wm2) <-c("variable","time","group","value")
  
  wf2 <-as.data.frame(melt(wf))
  colnames(wf2) <-c("variable","time","group","value")
  
  wnm2 <-as.data.frame(melt(wnm))
  colnames(wnm2) <-c("variable","time","group","value")
  
  wnf2 <-as.data.frame(melt(wnf))
  colnames(wnf2) <-c("variable","time","group","value")
  
  wtot2 <-as.data.frame(melt(wtot))
  colnames(wtot2) <-c("variable","time","group","value")
  
  lam2 <-as.data.frame(melt(lam))
  colnames(lam2) <-c("variable","time","group","value")
  
  n2 <-rbind(nm2,nf2,z2,be2,ni2,wm2,wf2,wnm2,wnf2,wtot2,lam2)

####################### Write data to global environment #######################
  biodata <- list("n.group"= n.group, "s.group"= s.group,"n.timestep"=n.timestep, 
                  "n.states"= n.states, "s.surv" = s.surv, "S"=S, 
                  "p.litter"= p.litter, "L"=L, "n.litter"= n.litter, 
                  "sr.litter" = sr.litter, "N_long" = n2, "breed"= state_breed, 
                  "nf.init" = nfa1, "nm.init"= nma1, "nf" = nf, "nm" = nm,"z" = z,
                  "be" = be, "ni" = ni, "wm" = wm, "wf" = wf, "wnm" = wnm, "wnf" = wnf,
                  "wtot" = wtot, "lam" = lam)
  list2env(biodata,.GlobalEnv)
}


