##############################################################################
#
# Runs the particle Metropolis-Hastings algorithm from different number 
# of observations generated from a LGSS model.
#
# Copyright (C) 2017 Johan Dahlin < liu (at) johandahlin.com.nospam >
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
##############################################################################

# Import helper
source("../helpers/dataGeneration.R")
source("../helpers/stateEstimation.R")
source("../helpers/parameterEstimation.R")

# Should the results be loaded from file (to quickly generate plots)
loadSavedWorkspace <- FALSE


##############################################################################
# Define the model
##############################################################################

# Here, we use the following model
#
# x[tt+1] = phi   * x[tt] + sigmav * v[tt]
# y[tt]   = x[tt]         + sigmae * e[tt]
#
# where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

# Set the parameters of the model
phi <- 0.75
sigmav <- 1.00
sigmae <- 0.10

# Set the number of time steps to simulate
T <- 250

# Set the initial state
initialState <- 0


##############################################################################
# PMH settings
##############################################################################

# The inital guess of the parameter
initialPhi <- 0.50

# No. particles in the particle filter
noParticles <- 100

# The length of the burn-in and the no. iterations of PMH ( noBurnInIterations < noIterations )
noBurnInIterations <- 1000
noIterations <- 5000

# Step size in the random walk proposal
stepSize <- 0.10


##############################################################################
# Loop over different data lengths
##############################################################################

TT <- c(10, 20, 50, 100, 200, 500)
Tmean <- matrix(0, nrow = length(TT), ncol = 1)
Tvar <- matrix(0, nrow = length(TT), ncol = 1)

if (loadSavedWorkspace) {
  load("../savedWorkspaces/example2-lgss-varyingT.RData")
} else {
  for (i in 1:length(TT)) {
    
    # Set the random seed to replicate results in tutorial
    set.seed(10)
    
    # Generate data
    data <- generateData(c(phi, sigmav, sigmae), TT[i], initialState)
    
    # Run the PMH algorithm
    res <-
      particleMetropolisHastings(
        data$y,
        initialPhi,
        sigmav,
        sigmae,
        noParticles,
        initialState,
        noIterations,
        stepSize
      )
    
    Tmean[i] <- mean(res[noBurnInIterations:noIterations])
    Tvar[i]  <- var(res[noBurnInIterations:noIterations])
  }
}


##############################################################################
# Save workspace and print results
##############################################################################

# Save workspace
if (!loadSavedWorkspace) {
  save.image("../savedWorkspaces/example2-lgss-varyingT.RData")
}

# Print the results to screen (no. observations, posterior mean, posterior variance)
cbind(TT, Tmean, Tvar)

# [1,]  10 0.5955020 0.0399332238
# [2,]  20 0.7943218 0.0127682838
# [3,]  50 0.7649620 0.0089581720
# [4,] 100 0.7269762 0.0060643002
# [5,] 200 0.6960883 0.0026939445
# [6,] 500 0.7185719 0.0009992732

##############################################################################
# End of file
##############################################################################
