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

# Save plot to file
savePlotToFile <- FALSE


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
noParticles <- c(10, 20, 50, 100, 200, 500, 1000)

# No. repetitions of log-likelihood estimate
noSimulations <- 1000

# Pre-allocate vectors
logLikelihoodEstimates <- matrix(0, nrow = length(noParticles), ncol = noSimulations)
logLikelihoodVariance <- rep(0, length(noParticles))
computationalTimePerSample <- rep(0, length(noParticles))

# Main loop
if (loadSavedWorkspace) {
  load("../savedWorkspaces/example4-sv-varyingN.RData")
  } else {
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
stepSize <- matrix(
  c(
    0.137255431,
    -0.0016258103,
    0.0015047492,
    -0.0016258103,
    0.0004802053,
    -0.0009973058,
    0.0015047492,
    -0.0009973058,
    0.0031307062
  ),
  ncol = 3,
  nrow = 3
)
stepSize <- 0.8 * stepSize

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
    res <- particleMetropolisHastingsSVmodel(y, initialTheta, noParticles[k], noIterations, stepSize)
    
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
  acf_mu <- acf(resTheta[k, , 1], plot = FALSE, lag.max = 200)
  acf_phi <- acf(resTheta[k, , 2], plot = FALSE, lag.max = 200)
  acf_sigmav <- acf(resTheta[k, , 3], plot = FALSE, lag.max = 200)
  
  resThetaIACT[k, ] <- 1 + 2 * c(sum(acf_mu$acf), sum(acf_phi$acf), sum(acf_sigmav$acf))
  resThetaIACTperSecond[k, ] <- resThetaIACT[k, ] / computationalTimePerIteration[k]
}

print(rbind(noParticles, apply(resThetaIACT, 1, max), apply(resThetaIACTperSecond, 1, max)))



# # Export plot to file
# if (savePlotToFile) {
#   cairo_pdf("figures/example4-sv-varyingN.pdf",
#             height = 10,
#             width = 8)
# }
# 
# layout(matrix(1:14, 7, 2, byrow = FALSE))
# par(mar = c(4, 5, 0, 0))
# 
# for (k in 1:length(noParticles)) {
#   xlab = ""
#   if (k == length(noParticles)) { xlab = "log-likelihood"}
#   hist(
#     logLikelihoodEstimates[k, ],
#     breaks = floor(sqrt(noSimulations)),
#     col = rgb(t(col2rgb("#1B9E77")) / 256, alpha = 0.25),
#     border = NA,
#     xlab = xlab,
#     ylab = "",
#     xlim = c(-730, -680),
#     main = "",
#     freq = FALSE
#   )
# }
# 
# # Close the plotting device
# if (savePlotToFile) {
#   dev.off()
# }


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
