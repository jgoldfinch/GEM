# Goal Efficient Monitoring (GEM)
### Simulation of populations, observations, and predictions for GEM approach

Code for Golding 2022. Goal efficient monitoring: an approach to monitoring as information changes – dissertation Chapter 1 available at https://scholarworks.umt.edu/cgi/viewcontent.cgi?article=13112&context=etd.

The example provided below (in the master file) is based on a hypothetical Canada lynx (*Lynx canadensis*) population. 

## Programs needed to run:
R, JAGS

## Files needed to run:
gem_00_master.R  
gem_01_load_packages.R  
gem_02_sim_bio.R  
gem_03_sim_tx.R  
gem_04_sim_obs.R  
gem_05_format_init.R  
gem_06_run_model.R  
gem_07_format_results.R  
gem_08_comp_diagnostics.R   
gem_09_tbl_diagnostics.R  
gem_10_format_tx.R  
gem_11_calc_tx.R  
gem_12_calc_tx_diagnostics.R  
gem_13_tbl_tx_diagnostics.R  


## Complete the following steps:

### 1) Open gem_00_master.R 
All R code will be run from here and line numbers below refer to this file.  

### 2) Modify the appropriate working directory 
Edit line 11 (or set your working directory with the method you prefer) to the location that contains the files listed under “Files needed to run simulation” above.  

### 3) Modify n.reps 
Edit line 29 according to how many times you want each simulation to repeat (the paper uses 100, but we recommend testing the functionality with a low # first and then increasing it).   
  * Note that the replications of the simulation are set to loop over multiple simulations, observations, and model runs, as well as loop past model runs that result in an error. This is to ensure starting values of the simulated population and starting values in the model are compatible (i.e., avoid invalid parent node error from JAGS). To get a full 100 simulations, you may need to run 400-500 replications. If you run more replicates as suggested, be sure to adjust the beginning of the file names (rep1_…) to start at 1 and be consecutive, as that is required for the code. Addressing this issue is in process.  

### 4) Adjust values in biological process - gem_02_sim_bio.R
Simulate populations to observe. If desired, modify the following input values within the gem_sim_bio function (lines 34 – 36) – values provided in the gem_00_master.R file are for the main simulation in the publication:  
  * n.group = number of populations (groups) (values 1 and above)  
  * s.group = number of females and males within each population (group) (values 1 and above)  
  * n.timestep = number of time steps the population is simulated for (values 1 and above)  
  * n.states = number of GEM population states (values 1 and above)  
  * s.surv = individual survival probability (values between 0 and 1)  
  * p.litter = probability of having a litter (values between 0 and 1)  
  * n.litter = number of individuals in a litter (values 1 and above)  
  * sr.litter = sex ratio of litter (values between 0 and 1)    

### 5) Calculate the GEM state transition probabilities for each time step - gem_03_sim_tx.R
Probabilities are caclculated from the values for the populations simulated with gem_02_sim_bio.R.  

### 6) Generate observations of simulated populations using GEM sampling rules - gem_04_sim_obs.R
If desired, modify the following input values within the gem_sim_obs function (lines 42 – 43) – values provided in the gem_00_master.R file are for the main simulation in the publication:
  * p = detection probability (values between 0 and 1)  
  * n.visits = number of repeat visits to each population (group) to survey during a single season (values 1 and above)  
  * pgenetic = detection probability of genetic sign (values between 0 and 1)  

### 7) Format initial values for the GEM model - gem_05_format_init.R
Generate initial values for the GEM model based on observed data (generated by the gem_04_sim_obs function)  

### 8) Run the GEM model - gem_06_run_model.R
Run model Bayesian integrated population model (BIPM) for the GEM approach, with input derived from observation (generated by the gem_04_sim_obs function). If desired, modify the following input values within the gem_06_run_model function (lines 55 – 59) – values provided in the gem_00_master.R file are for the main simulation in the publication:  
  * sr = sex ratio of the population (values between 0 and 1)
  * n.iter = number of iterations to run each MCMC chain
  * n.burnin = number of burn iterations discarded in MCMC chain runs
  * s.surv.init = initial value for survival (values between 0 and 1)  

### 9) Format GEM model results - gem_07_format_results.R
Format model results. All values inputted into the function should be consistent with that was using in steps 4 – 8 of this document.  

### 10) GEM model diagnostics - gem_08_comp_diagnostics.R
Compute model diagnostics with input from the gem_07_format_results function using output from gem_02_sim_bio.R as comparison for the true values.    
 
### 11) GEM model diagnostic tables - gem_09_comp_diagnostics.R
Write diagnostics tables from gem_08_comp_diagnostics.R to files.   

### 12) Examine file output  
At this point you should examine the file output through step 11. Starting with the summary tables may be most useful (e.g., nf_results_summary.csv). Look at the convergence column in the table. If all of the replicates, including all time steps within each replicate, have an average Rhat value of less than 1.06, then that variable is considered converged and given a 1. If not, convergence is listed as 0. If you are interested in what values did not converge you can look at the detailed information in all results files (e.g., nf_results_all.csv). If most of the variables in a replicate did not converge, you should re-run replicates (repeat steps 2 – 11) and discard that non-converged replicate. Code for automating this process is in progress.    

### 13) Format GEM state transition output - gem_10_format_tx.R
Format GEM state transition outputs from the model calculations from gem_06_run_model.R. All values inputted into the function should be consistent with values used in steps 4 – 8 of this document.  

### 14) Calculate predicted GEM state transitions - gem_11_calc_tx.R
Calculate predicted transitions from formatted model output from gem_10_format_tx.R.   

### 15) Calculate GEM state transitions diagnostics - gem_12_calc_tx_diagnostics.R
Calculate diagnostics for predicted GEM transitions from formatted from gem_11_calc_tx.R using output from gem_03_sim_tx.R as comparison for the true values.    

### 16) Write transitions diagnostics tables - gem_13_tbl_tx_diagnostics.R
Write diagnostic tables from gem_12_calc_tx_diagnostics.R to files.  

