##############################################################################
#
# Example of particle Metropolis-Hastings in a stochastic volatility model
# The effect on mixing while varying N.
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

# Import libraries
library("Quandl")
library("mvtnorm")

# Import helpers
source("../helpers/stateEstimation.R")
source("../helpers/parameterEstimation.R")
source("../helpers/plotting.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

# Should the results be loaded from file (to quickly generate plots)
loadSavedWorkspace <- FALSE

# Should the proposals be tuned by a pilot run
tuneProposals <- FALSE

# Should we use the tuned proposals (requires "../savedWorkspaces/example4-sv-varyingN-proposals.RData")
useTunedProposals <- TRUE


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
# Likelihood estimation using particle filter
##############################################################################

# True parameters estimated in example5-sv.R
theta <- c(-0.12, 0.96, 0.17)

# No. particles in the particle filter to try out
noParticles <- c(50, 100, 200, 300, 400, 500)

# No. repetitions of log-likelihood estimate
noSimulations <- 1000

# Pre-allocate vectors
logLikelihoodEstimates <- matrix(0, nrow = length(noParticles), ncol = noSimulations)
logLikelihoodVariance <- rep(0, length(noParticles))
computationalTimePerSample <- rep(0, length(noParticles))

# Main loop
if (!loadSavedWorkspace) {
  for (k in 1:length(noParticles)) {
    # Save the current time
    ptm <- proc.time()
    
    for (i in 1:noSimulations) {
      # Run the particle filter
      res <- particleFilterSVmodel(y, theta, noParticles[k])
      
      # Save the log-Likelihood estimate
      logLikelihoodEstimates[k, i] <- res$logLikelihood
    }
    
    # Compute the variance of the log-likelihood and computational time per sample
    logLikelihoodVariance[k] <- var(logLikelihoodEstimates[k, ])
    computationalTimePerSample[k] <- (proc.time() - ptm)[3] / noSimulations
    
    # Print to screen
    print(paste(paste(paste(paste("Simulation: ", k, sep = ""), " of ", sep = ""), length(noParticles), sep = ""), " completed.", sep = ""))
    print(paste(paste(paste(paste("No. particles: ", noParticles[k], sep = ""), " requires ", sep = ""), computationalTimePerSample[k], sep = ""), " seconds for computing one sample.", sep = ""))
  }
}


##############################################################################
# Parameter estimation using PMH
##############################################################################

# The inital guess of the parameter (use the estimate of the posterior mean to
# accelerated the algorithm, i.e., so less PMH iterations can be used).
initialTheta <- theta

# The length of the burn-in and the no. iterations of PMH ( noBurnInIterations < noIterations )
noBurnInIterations <- 2500
noIterations <- 7500

# The standard deviation in the random walk proposal
if (useTunedProposals) {
  load(file = "../savedWorkspaces/example4-sv-varyingN-proposals.RData")
} else {
  proposals <- array(0, dim = c(length(noParticles), 3, 3))
  for (k in 1:length(noParticles)) {
    proposals[k, , ] <- diag(c(0.10, 0.01, 0.05) ^ 2)
  }
}

# Main loop
if (loadSavedWorkspace) {
  load("../savedWorkspaces/example4-sv-varyingN.RData")
} else {
  resTheta <- array(0, dim = c(length(noParticles), noIterations - noBurnInIterations + 1, 3))
  computationalTimePerIteration <- rep(0, length(noParticles))
  acceptProbability <- rep(0, length(noParticles))
  
  for (k in 1:length(noParticles)) {
    # Save the current time
    ptm <- proc.time()
    
    # Run the PMH algorithm
    res <- particleMetropolisHastingsSVmodel(y, initialTheta, noParticles[k], noIterations, stepSize = proposals[k, ,])
    
    # Save the parameter trace
    resTheta[k, ,] <- res$theta[noBurnInIterations:noIterations,]
    
    # Compute acceptance probability and computational time per sample
    computationalTimePerIteration[k] <- (proc.time() - ptm)[3] / noIterations    
    acceptProbability[k] <- mean(res$proposedThetaAccepted[noBurnInIterations:noIterations])
    
    # Print to screen
    print(paste(paste(paste(paste("Simulation: ", k, sep = ""), " of ", sep = ""), length(noParticles), sep = ""), " completed.", sep = ""))
  }
}


##############################################################################
# Post-processing (computing IACT and IACT * time)
##############################################################################

resThetaIACT <- matrix(0, nrow = length(noParticles), ncol = 3)
resThetaIACTperSecond <- matrix(0, nrow = length(noParticles), ncol = 3)

for (k in 1:length(noParticles)) {
  acf_mu <- acf(resTheta[k, , 1], plot = FALSE, lag.max = 100)
  acf_phi <- acf(resTheta[k, , 2], plot = FALSE, lag.max = 100)
  acf_sigmav <- acf(resTheta[k, , 3], plot = FALSE, lag.max = 100)
  
  resThetaIACT[k, ] <- 1 + 2 * c(sum(acf_mu$acf), sum(acf_phi$acf), sum(acf_sigmav$acf))
  resThetaIACTperSecond[k, ] <- resThetaIACT[k, ] / computationalTimePerIteration[k]
}

table <- rbind(noParticles, logLikelihoodVariance, 100 * acceptProbability, apply(resThetaIACT, 1, max), apply(resThetaIACT, 1, max) * computationalTimePerIteration, computationalTimePerIteration)
table <- round(table, 2)
print(table)


##############################################################################
# Tune the PMH proposal using a pilot run
##############################################################################

if (tuneProposals) {
  proposals <- array(0, dim = c(length(noParticles), 3, 3))
  
  for (k in 1:length(noParticles)) {
    proposals[k, , ] <- cov(resTheta[k, , ]) * 2.562^2 / 3
  }
  save(proposals, file = "../savedWorkspaces/example4-sv-varyingN-proposals.RData")
}


##############################################################################
# Compute and save the results
##############################################################################


# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("../savedWorkspaces/example4-sv-varyingN.RData")
}


##############################################################################
# End of file
##############################################################################
