##############################################################################
#
# Example of particle Metropolis-Hastings
# in a stochastic volatility model
#
#
# Copyright (C) 2017 Johan Dahlin < liu (at) johandahlin.se >
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

# Import library
library("Quandl")

# Import helpers
source("helpers/stateEstimation.R")
source("helpers/parameterEstimation.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

# Should the results be loaded from file (to quickly generate plots)
loadSavedWorkspace <- FALSE

# Save plot to file
savePlotToFile <- TRUE


##############################################################################
# Define the model
##############################################################################

# Here, we use the following model
#
# x[tt+1] = phi  * x[tt] + sigma   * v[tt]
# y[tt]   = beta * exp( xt[tt]/2 ) * e[tt]
#
# where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

# Set the number of time steps to simulate
T <- 500


##############################################################################
# Load data
##############################################################################
d <-
  Quandl(
    "NASDAQOMX/OMXS30",
    start_date = "2012-01-02",
    end_date = "2014-01-02",
    type = "zoo"
  )
y <- as.numeric(100 * diff(log(d$"Index Value")))


##############################################################################
# Parameter estimation using PMH
##############################################################################

# The inital guess of the parameter
initialTheta <- c(0, 0.9, 0.2)

# No. particles in the particle filter ( choose noParticles ~ T )
noParticles <- 500

# The length of the burn-in and the no. iterations of PMH ( noBurnInIterations < noIterations )
noBurnInIterations <- 2500
noIterations <- 7500

# The standard deviation in the random walk proposal
stepSize <- matrix(
  c(
    0.137255431,-0.0016258103,
    0.0015047492,-0.0016258103,
    0.0004802053,-0.0009973058,
    0.0015047492,-0.0009973058,
    0.0031307062
  ),
  ncol = 3,
  nrow = 3
)
stepSize <- 0.8 * stepSize

# Run the PMH algorithm
if (loadSavedWorkspace) {
  load("savedWorkspaces/example4-sv.RData")
} else {
  res <-
    particleMetropolisHastingsSVmodel(y, initialTheta, noParticles, noIterations, stepSize)
}

##############################################################################
# Plot the results
##############################################################################

# Export plot to file
if (savePlotToFile) {
  cairo_pdf("figures/example4-sv.pdf",
            height = 10,
            width = 8)
}

makePlotsParticleMetropolisHastingsSVModel(y, res, noBurnInIterations, noIterations, nPlot)

# Close the plotting device
if (savePlotToFile) {
  dev.off()
}


##############################################################################
# Compute and save the results
##############################################################################

# Print the estimate of the posterior mean and standard deviation
print(thhat)
print(thhatSD)

#[1] -0.1087619  0.9668493  0.1593125
#[1] 0.24976843 0.02232583 0.05356500

# Compute an estimate of the IACT using the first 100 ACF coefficients
(iact <-
   1 + 2 * c(sum(muACF$acf), sum(phiACF$acf), sum(sigmavACF$acf)))
# [1] 13.28575 26.50253 23.31947

# Estimate the covariance of the posterior to tune the proposal
estCov <- var(resTh)

# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("savedWorkspaces/example4-sv.RData")
}


##############################################################################
# End of file
##############################################################################
