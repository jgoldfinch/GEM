# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

# MASTER FILE

################################################################################
## 01. Setup

# Set working directory
setwd("C:/Users/jgolding/Documents/publications/gem_pub_1/GitHub")

# Load functions 
source("./gem_01_load_packages.R")
source("./gem_02_sim_bio.R")
source("./gem_03_sim_tx.R")
source("./gem_04_sim_obs.R")
source("./gem_05_format_init.R")
source("./gem_06_run_model.R")
source("./gem_07_format_results.R")
source("./gem_08_comp_diagnostics.R")
source("./gem_09_tbl_diagnostics.R")
source("./gem_10_format_tx.R")
source("./gem_11_calc_tx.R")
source("./gem_12_calc_tx_diagnostics.R")
source("./gem_13_tbl_tx_diagnostics.R")

################################################################################
n.reps <-1

for(q in 1:n.reps){

## 02. Biological process - simulate populations to observe
  gem_sim_bio(n.group = 1, s.group = 7, n.timestep = 11,
              n.states = 4, s.surv = 0.7, p.litter = 0.5,
              n.litter = 2, sr.litter=0.5)

## 03. Calculate the GEM state transition probabilities 
  gem_sim_tx(wtot=wtot, wnm=wnm, wnf=wnf, n.group, n.timestep, n.states)

## 04. Observe data (observation process) 
  gem_sim_obs(p=0.63, nf=nf, nm=nm, n.visits=3, n.group, n.timestep, n.states,
              z=z, pgenetic=0.5)

## 05. Generate initial values based on observation data
  gem_format_init(ygf=ygf, ygm=ygm)

## 06. Run model
  
  # Parameters for the model to keep track of
  params <- c("nf","nm", "be", "wm", "wf", "wnm", "wnf", "wtot", "z", "lam", 
              "p", "pgenetic", "s.surv")

  skip_to_next <- FALSE
  tryCatch(gem_run_model(yp = yp, yc = yc, yg = yg, ygf = ygf, params = params,
                         n.group = n.group, s.group = s.group,
                         n.timestep = n.timestep, n.visits = n.visits,
                         n.litter = n.litter, sr = 0.5, n.iter = 400000,
                         n.burnin = 10000, s.surv.init = 0.7),
  error = function(e) {suppressWarnings(skip_to_next <<- TRUE)})
  if(skip_to_next) { next }

  finalout <-list(out, N_long, TPM_long, yp, yc, yg, ygf, ygm)
  name <-paste("rep",q,"_", n.timestep, "t", "_", n.group,"g", "_", s.group,"st", 
               "_", n.states, "states", "_", s.surv, "surv", "_", p.litter, 
               "plitter", "_", n.litter, "perlitter","_", sr.litter, "sr","_",
               p,"p", "_", n.visits, "visits", sep="")
  save(finalout,file= paste0 ("./", name))
}

## 07. Format population prediction results
gem_format_results(n.group = n.group, s.group = s.group, n.timestep = n.timestep, 
                   n.states = n.states, s.surv = s.surv, p.litter = p.litter, 
                   n.litter = n.litter, sr.litter = sr.litter, p = p, 
                   pgenetic = pgenetic, n.visits = n.visits, n.reps = 1, 
                   s.surv.init = 0.7)

## 08. Summarize population diagnostics 
gem_comp_diagnostics(p_results = p_results, pg_results = pg_results, 
                     nf_results = nf_results, nm_results = nm_results, 
                     be_results = be_results, wf_results = wf_results, 
                     wm_results = wm_results, wnm_results = wnm_results,
                     wnf_results = wnf_results, wtot_results = wtot_results, 
                     z_results = z_results, surv_results = surv_results)

## 09. Write tables for population diagnostics
gem_tbl_ds(p_results = p_results, pg_results = pg_results, 
           nf_results = nf_results, nm_results = nm_results, 
           be_results = be_results, wf_results = wf_results, 
           wm_results = wm_results, wnm_results = wnm_results,
           wnf_results = wnf_results, wtot_results = wtot_results, 
           z_results = z_results, s_results = s_results)

## 10. Load data for transition calculations
gem_format_tx(n.group = n.group, s.group = s.group, n.timestep = n.timestep, 
              n.states = n.states, s.surv = s.surv, p.litter = p.litter, 
              n.litter = n.litter, sr.litter = sr.litter, p = p, 
              n.visits = n.visits, n.reps = 1)

## 11. Calculate predicted transitions
gem_calc_tx(prob2 = prob2, prob3 = prob3, prob4 = prob4, wtot = wtot, wnf = wnf,
            wnm = wnm, s_results = s_results, n.timestep = n.timestep, n.reps = 1)

## 12. Calculate transition diagnostics
gem_calc_txd(tprob=tprob, twtot=twtot, twnf = twnf, twnm = twnm, 
             psi21p = psi21p, psi22p = psi22p, psi31p = psi31p, 
             psi32p = psi32p, psi33p = psi33p, psi41p = psi41p, psi42p = psi42p, 
             psi43p = psi43p, psi44p = psi44p)

## 13. Write transition diagnostics to tables 
gem_tbl_txd(psi44_results = psi44_results, psi43_results = psi43_results, 
            psi42_results = psi42_results, psi41_results = psi41_results,
            psi33_results = psi33_results, psi32_results = psi32_results, 
            psi31_results = psi31_results, psi22_results = psi22_results,
            psi21_results = psi21_results)


