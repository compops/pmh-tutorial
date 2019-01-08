##############################################################################
# Parameter estimation using particle Metropolis-Hastings in a reparameterised version of a
# stochastic volatility model with a proposal adapted from a pilot run.
#
# Johan Dahlin <uni (at) johandahlin.com.nospam>
# Documentation at https://github.com/compops/pmh-tutorial
# Published under GNU General Public License
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
nPlot <- 2500

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
resTh <- res$theta[noBurnInIterations:noIterations, ]
thhat   <- colMeans(resTh)
thhatSD <- apply(resTh, 2, sd)

#[1] -0.1550373  0.9601144  0.1742736
#[1] 0.23637116 0.02239614 0.05701460

# Compute an estimate of the IACT using the first 100 ACF coefficients
print(iact)
# [1] 21.93670 28.96783 16.65938

# Estimate the covariance of the posterior to tune the proposal
resThTransformed <- res$thetaTransformed[noBurnInIterations:noIterations,]
estCov <- var(resThTransformed)

# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("savedWorkspaces/example5-sv.RData")
}