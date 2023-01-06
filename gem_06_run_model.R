# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_run_model
# Description: function to run Bayesian GEM model

################################# Arguments ####################################

# yp: 
#       Observation of presence. Array generated with gem_sim_obs function. Whole
#       number.
# yc: 
#       Observation of counts. Array generated with gem_sim_obs function. Whole number.
# yg: 
#       Observation of genetic sign (subset of yc). Array generated with 
#       gem_sim_obs function. Whole number.
# ygf: 
#       Females observed from genetic sign (subset of yg). Array generated with 
#       gem_sim_obs function. Whole number.
# params: 
#       Parameters for the model to track. Vector of characters.
# n.group: 
#       Number of populations that were observed. Whole number.
# s.group: 
#       Starting number of females and males in each population. Whole number.
# n.timestep: 
#       Number of time steps you want to observe. Whole number.
# n.visits: 
#       Number of visits to each group to survey.
# n.litter: 
#       Number of individuals per litter. Number between 0 and 1.
# sr: 
#       Sex ratio of males to females per litter. Number between 0 and 1.
# n.iter: 
#       Number of iterations to run each MCMC chain. Whole number.
# n.burnin: 
#       Number of burn iterations discarded in MCMC chain runs. Whole number.
# s.surv.init: 
#       Initial value for survival. Number between 0 and 1.

################################## Output ######################################

# out: JAGS UI object with Bayesian model results that summarize info for the 
# parameters set with the params term

################################# Function #####################################

gem_run_model <- function(yp, yc, yg, ygf, params, n.group, s.group, n.timestep, 
                          n.visits, n.litter, sr, n.iter, n.burnin, s.surv.init){

sink("Model.txt")
cat("
    model{
    ## Priors 
    # Detection
    p ~ dunif(0, 1)
    s.surv ~ dunif(0.1,1)
    p.litter ~ dunif(0,1)
    pgenetic ~ dunif(0,1)
    lambda ~ dgamma(0.001,0.001)

	  for(i in 1:n.group){  
	  for(j in 1:n.timestep){
		  S[i,j] <- s.surv 
		  L[i,j] <- p.litter 
		  sexratio[i,j] <- sr
		  nl[i,j] <-n.litter
    }
	  }
    
    for(i in 1:n.group){
      nfi[i] ~ dpois(lambda) 
	    nmi[i] ~ dpois(lambda)
    }

    ## Biological model
    for(i in 1:n.group){ 
       
      
      ## Time 1
      
      # Entering individuals/state
	    nf[1,1,i] <-nfi[i] 
	    nm[1,1,i] <-nmi[i] 
	    z[1,1,i] <- ifelse(nm[1,1,i] >= 1 && nf[1,1,i] >=1,4,
                  ifelse(nm[1,1,i] > 1 && nf[1,1,i] == 0 || nf[1,1,i] > 1 && nm[1,1,i] == 0,3,
                  ifelse(nm[1,1,i] == 1 && nf[1,1,i] == 0 || nf[1,1,i] == 1 && nm[1,1,i] == 0,2,
                  ifelse(nm[1,1,i] == 0 && nf[1,1,i] == 0,1,99))))
                  
      z2[1,1,i] <- ifelse(z[1,1,i] > 1,1,0)


	    # Breeding possible
      state_breed[1,i] <-ifelse(z[1,1,i]==4,1,0)
      litter_prob[1,i] <-inprod(state_breed[1,i],L[i,1])

      # Births
      be[1,1,i] ~ dbin(litter_prob[1,i], nf[1,1,i]) 
      
      # New individuals
      ni[1,1,i] <- inprod(be[1,1,i],nl[i,1])

      wm[1,1,i] ~ dbin(sexratio[i,1], ni[1,1,i])
      wf[1,1,i] <- ni[1,1,i] - wm[1,1,i]

      # Totals
      wnf[1,1,i] <- wf[1,1,i] + nf[1,1,i]
      wnm[1,1,i] <- wm[1,1,i] + nm[1,1,i]
      wtot[1,1,i] <- wnm[1,1,i] + wnf[1,1,i]
      
    for(j in 2:n.timestep){ 
      ## Time 2+
      # Entering individuals/state
	    nf[1,j,i] ~ dbin(S[i,j], wnf[1,j-1,i])
	    nm[1,j,i] ~ dbin(S[i,j], wnm[1,j-1,i]) 
	    z[1,j,i] <- ifelse(nm[1,j,i] >= 1 && nf[1,j,i] >=1,4,
                  ifelse(nm[1,j,i] >  1 && nf[1,j,i] == 0 || nf[1,j,i] >  1 && nm[1,j,i] == 0,3,
                  ifelse(nm[1,j,i] == 1 && nf[1,j,i] == 0 || nf[1,j,i] == 1 && nm[1,j,i] == 0,2,
                  ifelse(nm[1,j,i] == 0 && nf[1,j,i] == 0,1,99))))
                  
      z2[1,j,i] <- ifelse(z[1,j,i] > 1,1,0)
    
	    # Breeding possible
      state_breed[j,i] <-ifelse(z[1,j,i]==4,1,0)
      litter_prob[j,i] <-inprod(state_breed[j,i],L[i,j])
      
      # Births
      be[1,j,i] ~ dbin(litter_prob[j,i], nf[1,j,i]) 
      
      # New individuals
      ni[1,j,i] <- inprod(be[1,j,i],nl[i,j])
      
      wm[1,j,i] ~ dbin(sexratio[i,j], ni[1,j,i])
      wf[1,j,i] <- ni[1,j,i] - wm[1,j,i]
      
      # Totals
      wnf[1,j,i] <- wf[1,j,i] + nf[1,j,i]
      wnm[1,j,i] <- wm[1,j,i] + nm[1,j,i]
      wtot[1,j,i] <- wnm[1,j,i] + wnf[1,j,i]
      
      # Lambda
      lam[j,i] <- (wtot[1,j,i])/(wtot[1,j-1,i])
      
    }
    }
 
    ## Observation model
    for (k in 1:n.visits){ 
    for(i in 1:n.group){ 
      
    ## Time 1
    # Is the species present?
      yp[k,1,i] ~ dbin((1-((1-p)^(nf[1,1,i]+nm[1,1,i]))), z2[1,1,i])
    
    for(j in 2:n.timestep){  
    ## Time 2 +
    # Is the species present?
      yp[k,j,i] ~ dbin((1-((1-p)^(nf[1,j,i]+nm[1,j,i]))), z2[1,j,i])
      
    # Are multiple individuals present?
      yc[k,j,i] ~ dbin(p, nf[1,j,i]+nm[1,j,i])
    
    # Are females and males present?
      yg[k,j,i] ~ dbin(pgenetic, yc[k,j,i])
      ygf[k,j,i] ~ dhyper(nf[1,j,i],nm[1,j,i],yg[k,j,i],1)

      }
     }
    }
    
}
    ", fill = TRUE)
sink()

# Main data
data <- list(n.group = n.group, s.group = s.group, n.timestep = n.timestep,
             n.litter = n.litter, sr = sr, n.visits = n.visits,
             yp = yp, yc = yc, yg = yg, ygf = ygf)

# Initial value data
inits <- function() {
  list(nfi = nf.init, nmi = nm.init, s.surv = rtruncnorm(1, a=0, b=1, mean = s.surv.init, sd = 0.15))
}

# Run model and save output as object out
out<- jags(data=data, inits=inits, parameters.to.save = params, "Model.txt",
           n.chains=3, n.thin=10, n.iter=n.iter, n.burnin=n.burnin, n.adapt=5000,
           parallel = TRUE)

modeldata <- list("out"= out)
list2env(modeldata ,.GlobalEnv)
}
