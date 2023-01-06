# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_format_results
# Description:  Function to format result (observations, predictions, truth) 
#               from Bayesian model of variables (not model diagnostics) for 
#               plotting and calculating diagnostics

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
# pgenetic: 
#       Detection probability of genetic sign. Number between 0 and 1.
# n.visits: 
#       Number of visits to each group to survey.
# n.reps: 
#       Number of simulation replicates.
# s.surv.init: 
#       Initial value for survival. Number between 0 and 1.

################################## Output ######################################

# modeldata: list saved to global environment. Contains formatted prediction values
# from the GEM model for detection probability (p), detection probability of 
# genetic sign (pg), adult females (nf), adult males (nm), birth events (be), 
# new females (wf), new males (wm), total adult and new females (wnf), total 
# adult and new males (wnm), total adult and new individuals (wtot), population 
# occupancy (z), and survival (surv).

################################# Function #####################################

gem_format_results <- function(n.group, s.group, n.timestep, n.states, s.surv, p.litter, n.litter, 
                               sr.litter, p, pgenetic, n.visits, n.reps, s.surv.init){

  nf_results <- nm_results <- be_results <- wf_results <- wm_results <-pg_results <-
  wnf_results <- wnm_results <-wtot_results <-p_results <- z_results <-surv_results <-list()

  for (s in 1:n.reps){
    name <-paste("rep",s,"_", n.timestep, "t", "_", n.group,"g", "_", s.group,"st", "_", 
                 n.states, "states", "_", s.surv, "surv", "_", p.litter, "plitter",
                 "_", n.litter, "perlitter","_", sr.litter, "sr","_", p,"p", "_", 
                 n.visits, "visits", sep="")
    
    filename = paste0 ("./", name)
    load(filename)

    # NF results
    nf_results[[s]] <-melt(finalout[[1]]$mean$nf)
    colnames(nf_results[[s]]) <-c("var","t","group","pred")
    nf_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$nf)[,4]
    nf_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$nf)[,4]
    nf_results[[s]]$true <-finalout[[2]][finalout[[2]]$variable=="nf",][,4]
    nf_results[[s]]$rep <-s
    colnames(nf_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep")
    nf_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$nf,3))[,4]
    colnames(nf_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep", "rhat")
    nf_results2 <-bind_rows(nf_results)
    
    # NM results
    nm_results[[s]] <-melt(finalout[[1]]$mean$nm)
    colnames(nm_results[[s]]) <-c("var","t","group","pred")
    nm_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$nm)[,4]
    nm_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$nm)[,4]
    nm_results[[s]]$true <-finalout[[2]][finalout[[2]]$variable=="nm",][,4]
    nm_results[[s]]$rep <-s
    colnames(nm_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep")
    nm_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$nm,3))[,4]
    colnames(nm_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep", "rhat")
    nm_results2 <-bind_rows(nm_results)
    
    # BE results
    be_results[[s]] <-melt(finalout[[1]]$mean$be)
    colnames(be_results[[s]]) <-c("var","t","group","pred")
    be_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$be)[,4]
    be_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$be)[,4]
    be_results[[s]]$true <-finalout[[2]][finalout[[2]]$variable=="be",][,4]
    be_results[[s]]$rep <-s
    colnames(be_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep")
    be_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$be,3))[,4]
    colnames(be_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep", "rhat")
    be_results2 <-bind_rows(be_results)
    
    # wF results
    wf_results[[s]] <-melt(finalout[[1]]$mean$wf)
    colnames(wf_results[[s]]) <-c("var","t","group","pred")
    wf_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$wf)[,4]
    wf_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$wf)[,4]
    wf_results[[s]]$true <-finalout[[2]][finalout[[2]]$variable=="wf",][,4]
    wf_results[[s]]$rep <-s
    colnames(wf_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep")
    wf_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$wf,3))[,4]
    colnames(wf_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep", "rhat")
    wf_results2 <-bind_rows(wf_results)
    
    # wM results
    wm_results[[s]] <-melt(finalout[[1]]$mean$wm)
    colnames(wm_results[[s]]) <-c("var","t","group","pred")
    wm_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$wm)[,4]
    wm_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$wm)[,4]
    wm_results[[s]]$true <-finalout[[2]][finalout[[2]]$variable=="wm",][,4]
    wm_results[[s]]$rep <-s
    colnames(wm_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep")
    wm_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$wm,3))[,4]
    colnames(wm_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep", "rhat")
    wm_results2 <-bind_rows(wm_results)
    
    # wNF results
    wnf_results[[s]] <-melt(finalout[[1]]$mean$wnf)
    colnames(wnf_results[[s]]) <-c("var","t","group","pred")
    wnf_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$wnf)[,4]
    wnf_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$wnf)[,4]
    wnf_results[[s]]$true <-finalout[[2]][finalout[[2]]$variable=="wnf",][,4]
    wnf_results[[s]]$rep <-s
    colnames(wnf_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep")
    wnf_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$wnf,3))[,4]
    colnames(wnf_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep", "rhat")
    wnf_results2 <-bind_rows(wnf_results)
    
    # wNM results
    wnm_results[[s]] <-melt(finalout[[1]]$mean$wnm)
    colnames(wnm_results[[s]]) <-c("var","t","group","pred")
    wnm_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$wnm)[,4]
    wnm_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$wnm)[,4]
    wnm_results[[s]]$true <-finalout[[2]][finalout[[2]]$variable=="wnm",][,4]
    wnm_results[[s]]$rep <-s
    colnames(wnm_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep")
    wnm_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$wnm,3))[,4]
    colnames(wnm_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep", "rhat")
    wnm_results2 <-bind_rows(wnm_results)
  
    # WTOT results
    wtot_results[[s]] <-melt(finalout[[1]]$mean$wtot)
    colnames(wtot_results[[s]]) <-c("var","t","group","pred")
    wtot_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$wtot)[,4]
    wtot_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$wtot)[,4]
    wtot_results[[s]]$true <-finalout[[2]][finalout[[2]]$variable=="wtot",][,4]
    wtot_results[[s]]$rep <-s
    colnames(wtot_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep")
    wtot_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$wtot,3))[,4]
    colnames(wtot_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep", "rhat")
    wtot_results2 <-bind_rows(wtot_results)
    
    # Survival results
    surv_results[[s]] <-melt(finalout[[1]]$mean$s)
    colnames(surv_results[[s]]) <-c("pred")
    surv_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$s)[,1]
    surv_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$s)[,1]
    surv_results[[s]]$true <-0.70
    surv_results[[s]]$rep <-s
    surv_results[[s]]$var <-"s"
    surv_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$s,3))[,1]
    surv_results2 <-bind_rows(surv_results)
    
    # Z results
    z_results[[s]] <-melt(finalout[[1]]$mean$z)
    colnames(z_results[[s]]) <-c("var","t","group","pred")
    z_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$z)[,4]
    z_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$z)[,4]
    z_results[[s]]$true <-finalout[[2]][finalout[[2]]$variable=="z",][,4]
    z_results[[s]]$rep <-s
    colnames(z_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep")
    z_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$z,3))[,4]
    colnames(z_results[[s]]) <-c("var","t","group","pred","predlci","preduci","true","rep", "rhat")
    z_results2 <-bind_rows(z_results)
    
    # Detection results
    p_results[[s]] <-melt(finalout[[1]]$mean$p)
    colnames(p_results[[s]]) <-c("pred")
    p_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$p)[,1]
    p_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$p)[,1]
    colnames(p_results[[s]]) <-c("pred","predlci","preduci")
    p_results[[s]]$true <-p
    p_results[[s]]$rep <-s
    p_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$p,3))[,1]
    p_results2 <-bind_rows(p_results)
    
    # Genetic sign detection results
    pg_results[[s]] <-melt(finalout[[1]]$mean$pgenetic)
    colnames(pg_results[[s]]) <-c("pred")
    pg_results[[s]]$predlci <-melt(finalout[[1]]$q2.5$p)[,1]
    pg_results[[s]]$preduci <-melt(finalout[[1]]$q97.5$p)[,1]
    colnames(pg_results[[s]]) <-c("pred","predlci","preduci")
    pg_results[[s]]$true <-pgenetic
    pg_results[[s]]$rep <-s
    pg_results[[s]]$rhat <-melt(round(finalout[[1]]$Rhat$p,3))[,1]
    pg_results2 <-bind_rows(pg_results)
    
  }
  
  modeldata <- list("p_results"= p_results2, "pg_results"= pg_results2,
                    "nf_results"= nf_results2, "nm_results"= nm_results2, 
                    "be_results" = be_results2, "wf_results" = wf_results2,
                    "wm_results" = wm_results2, "wnf_results" = wnf_results2, 
                    "wnm_results" = wnm_results2, "wtot_results"= wtot_results2,
                    "z_results" = z_results2, "surv_results" = surv_results2)
  
  list2env(modeldata ,.GlobalEnv)
  
  }
