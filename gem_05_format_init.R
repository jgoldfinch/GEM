# Goal Efficient Monitoring (GEM)
# Author: Jessie Golding, Jamie Sanderlin
# Date: 07/12/2022

#################################### Intro #####################################

# Function name: gem_format_init
# Description: function to format and generate intial values for GEM model

################################# Arguments ####################################

# ygf: 
#       Number of females observed with genetic sign. Array generated with 
#       gem_sim_obs function.
# ymf: 
#       Number of males observed with genetic sign. Array generated with 
#       gem_sim_obs function.

################################## Output ######################################

# initdata: list saved to global environment. Contains model initial values.

################################# Function #####################################

gem_format_init <- function(ygf, ygm){

# Max presence or count for each time step to start initial values
nm.in <-apply(yp,c(2),max)*10
nm.init <-min(nm.in, na.rm = TRUE)

nf.in <-apply(yp,c(2),max)*10
nf.init <-min(nf.in, na.rm = TRUE)

initdata <- list("nf.init" = nf.init,"nm.init" = nm.init)
list2env(initdata,.GlobalEnv)

}