# Goal Efficient Monitoring
# Author: Jessie Golding, Jamie Sanderlin
# Date: 8/31/2021

# Function purpose: load packages (check if installed and if not install them)

################################################################################
## Load required packages
# Function from https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
## First specify the packages of interest
packages = c("jagsUI", "reshape2", "dplyr", "truncnorm")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)