# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_sim_tx
# Description:  function to calculate GEM transition probabilities, which are changes in
#               four population states (not present, single 
#               individual present, multiple individuals of a single sex present,
#               and multiple individuals with both sexes present)

################################# Arguments ####################################

# wtot:
#       Total individuals in the population at the end of each time step. 
#       Array generated with gem_sim_bio function.
# wnm:
#       Total males in the population at the end of each time step. 
#       Array generated with gem_sim_bio function.
# wnf:
#       Total females in the population at the end of each time step. 
#       Array generated with gem_sim_bio function.
# n.group:
#       Number of populations (groups). Whole number.
# n.timestep: 
#       Number of time steps. Whole number.
# n.states: 
#       Number of GEM population states. Whole number.

################################## Output ######################################

# txdata list saved to global environment. Contains transition probability matrix
# in multiple forms (array and long), and transition probability matricies

################################# Function #####################################

gem_sim_tx <- function(wtot, wnm, wnf, n.group, n.timestep, n.states){

### Time 1
### Create matrices/arrays to hold data

psi11 <-psi21 <-psi22 <-psi31 <-psi32 <-psi33 <-psi41 <- psi42 <-psi43 <-psi44 <-matrix(NA,n.timestep,n.group)

prob1 <-prob2 <-prob3 <-prob4 <-matrix(NA,n.timestep,n.group)

TPM <-array(data = NA, dim = c(n.states, n.states,n.timestep,n.group))


for (i in 1:n.group){         
  for(j in 2:n.timestep){  

# Time 1
prob1[1,i] <-ifelse(z[1,1,i]==1,1,0)
psi11[1,i] <-(dbinom (0, wtot[1,1,i], S[i,1])*rbinom(1,1,prob1[1,i]))/(rbinom(1,1,prob1[1,i]))
    
prob2[1,i] <-ifelse(z[1,1,i]==2,1,0)
psi21[1,i] <- (dbinom (0, wtot[1,1,i], S[i,1])*rbinom(1,1,prob2[1,i]))/(rbinom(1,1,prob2[1,i]))
psi22[1,i] <- 1-psi21[1,i]

prob3[1,i] <-ifelse(z[1,1,i]==3,1,0)
psi31[1,i] <- (dbinom (0, wtot[1,1,i], S[i,1])*rbinom(1,1,prob3[1,i]))/(rbinom(1,1,prob3[1,i]))
psi32[1,i] <- (((dbinom(1, wnm[1,1,i], S[i,1])*dbinom(0,wnf[1,1,i],S[i,1]))+
                (dbinom(1, wnf[1,1,i], S[i,1])*dbinom(0,wnm[1,1,i],S[i,1])))*rbinom(1,1,prob3[1,i]))/rbinom(1,1,prob3[1,i])
psi33[1,i] <- 1-psi32[1,i]-psi31[1,i]

prob4[1,i] <-ifelse(z[1,1,i]==4,1,0)
psi41[1,i] <- (dbinom (0, wtot[1,1,i], S[i,1])*rbinom(1,1,prob4[1,i]))/(rbinom(1,1,prob4[1,i]))
psi42[1,i] <- (dbinom (1, wtot[1,1,i], S[i,1])*rbinom(1,1,prob4[1,i]))/(rbinom(1,1,prob4[1,i]))
psi43[1,i] <- (((pbinom(1, wnm[1,1,i], S[i,1], lower.tail = FALSE)*dbinom(0,wnf[1,1,i],S[i,1]))+
                (pbinom(1, wnf[1,1,i], S[i,1], lower.tail = FALSE)*dbinom(0,wnm[1,1,i],S[i,1])))*rbinom(1,1,prob4[1,i]))/rbinom(1,1,prob4[1,i])
psi44[1,i] <- 1-psi43[1,i]-psi42[1,i]-psi41[1,i]

 TPM[,,1,i] <- matrix(c(

   psi11[1,i], 0, 0, 0,
   psi21[1,i], 1-psi21[1,i], 0, 0,
   psi31[1,i], psi32[1,i], psi33[1,i], 0,
   psi41[1,i], psi42[1,i], psi43[1,i], psi44[1,i]),

   nrow=n.states, ncol=n.states, byrow=TRUE)

## Time 2 and beyond
prob1[j,i] <-ifelse(wtot[1,j,i]==1,1,0) 
psi11[j,i] <-(dbinom (0, wtot[1,j,i], S[i,j])*rbinom(1,1,prob1[j,i]))/(rbinom(1,1,prob1[j,i]))
 
prob2[j,i] <-ifelse(wtot[1,j,i]==2,1,0)
psi21[j,i] <- (dbinom (0, wtot[1,j,i], S[i,j])*rbinom(1,1,prob2[j,i]))/(rbinom(1,1,prob2[j,i]))
psi22[j,i] <- 1-psi21[j,i]

prob3[j,i] <-ifelse(z[1,j,i]==3,1,0)
psi31[j,i] <- (dbinom (0, wtot[1,j,i], S[i,j])*rbinom(1,1,prob3[j,i]))/(rbinom(1,1,prob3[j,i]))
psi32[j,i] <- (((dbinom(1, wnm[1,j,i], S[i,j])*dbinom(0,wnf[1,j,i],S[i,j]))+
                (dbinom(1, wnf[1,j,i], S[i,j])*dbinom(0,wnm[1,j,i],S[i,j])))*rbinom(1,1,prob3[j,i]))/rbinom(1,1,prob3[j,i])
psi33[j,i] <- 1-psi32[j,i]-psi31[j,i]

prob4[j,i] <-ifelse(z[1,j,i]==4,1,0)
psi41[j,i] <- (dbinom (0, wtot[1,j,i], S[i,j])*rbinom(1,1,prob4[j,i]))/(rbinom(1,1,prob4[j,i]))
psi42[j,i] <- (dbinom (1, wtot[1,j,i], S[i,j])*rbinom(1,1,prob4[j,i]))/(rbinom(1,1,prob4[j,i]))
psi43[j,i] <- (((pbinom(1, wnm[1,j,i], S[i,j], lower.tail = FALSE)*dbinom(0,wnf[1,j,i],S[i,j]))+
                (pbinom(1, wnf[1,j,i], S[i,j], lower.tail = FALSE)*dbinom(0,wnm[1,j,i],S[i,j])))*rbinom(1,1,prob4[j,i]))/rbinom(1,1,prob4[j,i])
psi44[j,i] <- 1-psi43[j,i]-psi42[j,i]-psi41[j,i]

 TPM[,,j,i] <- matrix(c(
   
   psi11[j,i], 0, 0, 0,
   psi21[j,i], psi22[j,i], 0, 0,
   psi31[j,i], psi32[j,i], psi33[j,i], 0,
   psi41[j,i], psi42[j,i], psi43[j,i], psi44[j,i]),
   
   nrow=n.states, ncol=n.states, byrow=TRUE)
}
}

# Transition probabilities
TPM2 <- provideDimnames(TPM , sep = "_", 
                        base = list(as.character(1:n.states),as.character(1:n.states),
                                    as.character(1:(n.timestep)),as.character(1:n.group)))
t <-melt(TPM2)
colnames(t)<-c("startstate","endstate","time","group","transitionprob")


####################### Write data to global environment #######################
txdata <- list("TPM" = TPM, "TPM_long"= t,
               "psi41"= psi41, "psi42"= psi42, "psi43"= psi43,"psi44"= psi44,
               "psi31"= psi31, "psi32"= psi32, "psi33"= psi33,
               "psi21"= psi21, "psi22"= psi22, "psi11" = psi11, 
               "prob1"= prob1, "prob2"= prob2, "prob3" = prob3, "prob4" = prob4)
list2env(txdata,.GlobalEnv)
}
