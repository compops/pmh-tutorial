##############################################################################
# Parameter estimation using particle Metropolis-Hastings in a SV model
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
stepSize <- diag(c(0.10, 0.01, 0.05) ^ 2)

if (loadSavedWorkspace) {
  load("savedWorkspaces/example3-sv.RData")
} else {
  res <- particleMetropolisHastingsSVmodel(y, initialTheta, noParticles, noIterations, stepSize)
}

##############################################################################
# Plot the results
##############################################################################
if (savePlotToFile) {
  cairo_pdf("figures/example3-sv.pdf",
            height = 10,
            width = 8)
}

iact <- makePlotsParticleMetropolisHastingsSVModel(y, res, noBurnInIterations, noIterations, nPlot)

# Close the plotting device
if (savePlotToFile) {
  dev.off()
}

# Print the estimate of the posterior mean and standard deviation
print(thhat)
print(thhatSD)

#[1] -0.2337134  0.9708399  0.1498914
#[1] 0.37048000 0.02191359 0.05595271

# Compute an estimate of the IACT using the first 100 ACF coefficients
print(iact)
# [1] 135.19084  85.98935  65.80120

# Estimate the covariance of the posterior to tune the proposal
estCov <- var(resTh)
#               [,1]          [,2]          [,3]
# [1,]  0.137255431 -0.0016258103  0.0015047492
# [2,] -0.001625810  0.0004802053 -0.0009973058
# [3,]  0.001504749 -0.0009973058  0.0031307062

# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("savedWorkspaces/example3-sv.RData")
}