##############################################################################
# Parameter estimation using particle Metropolis-Hastings in a reparameterised version of a
# stochastic volatility model with a proposal adapted from a pilot run.
# (c) Johan Dahlin 2017 under MIT license <liu@johandahlin.com.nospam>
##############################################################################

library("Quandl")
library("mvtnorm")
source("helpers/stateEstimation.R")
source("helpers/parameterEstimation.R")
source("helpers/plotting.R")

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
# PMH
##############################################################################
initialTheta <- c(0, 0.9, 0.2)
noParticles <- 500
noBurnInIterations <- 2500
noIterations <- 7500
stepSize <- matrix(
  c(
    0.041871682,-0.001200581,-0.002706803,-0.001200581,
    0.054894707,-0.056321320,-0.002706803,-0.056321320,
    0.087342276
  ),
  ncol = 3,
  nrow = 3
)
stepSize <- 2.562^2 / 3 * stepSize

if (loadSavedWorkspace) {
  load("savedWorkspaces/example5-sv.RData")
} else {
  res <- particleMetropolisHastingsSVmodelReparameterised(y, initialTheta, noParticles, noIterations, stepSize)
}

##############################################################################
# Plot the results
##############################################################################
if (savePlotToFile) {
  cairo_pdf("figures/example5-sv.pdf", height = 10, width = 8)
}

iact <- makePlotsParticleMetropolisHastingsSVModel(y, res, noBurnInIterations, noIterations, nPlot)

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

#[1] -0.1466918  0.9577250  0.1813333
#[1] 0.21236553 0.02287443 0.05967315

# Compute an estimate of the IACT using the first 100 ACF coefficients
print(iact)
# [1] 12.55590 20.63091 16.94274

# Estimate the covariance of the posterior to tune the proposal
resThTransformed <- res$thetaTransformed[noBurnInIterations:noIterations,]
estCov <- var(resThTransformed)

# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("savedWorkspaces/example5-sv.RData")
}