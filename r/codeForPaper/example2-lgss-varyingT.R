##############################################################################
# Runs the particle Metropolis-Hastings algorithm from different number 
# of observations generated from a LGSS model.
# (c) Johan Dahlin 2017 under MIT license <liu@johandahlin.com.nospam>
##############################################################################

source("../helpers/dataGeneration.R")
source("../helpers/stateEstimation.R")
source("../helpers/parameterEstimation.R")

# Should the results be loaded from file (to quickly generate plots)
loadSavedWorkspace <- FALSE

##############################################################################
# Define the model and generate data
# x[t + 1] = phi * x[t] + sigmav * v[t],    v[t] ~ N(0, 1)
# y[t] = x[t] + sigmae * e[t],              e[t] ~ N(0, 1)
##############################################################################
phi <- 0.75
sigmav <- 1.00
sigmae <- 0.10
T <- 250
initialState <- 0

##############################################################################
# PMH
##############################################################################

initialPhi <- 0.50
noParticles <- 100
noBurnInIterations <- 1000
noIterations <- 5000
stepSize <- 0.10

# Loop over different data lengths
TT <- c(10, 20, 50, 100, 200, 500)
Tmean <- matrix(0, nrow = length(TT), ncol = 1)
Tvar <- matrix(0, nrow = length(TT), ncol = 1)

if (loadSavedWorkspace) {
  load("../savedWorkspaces/example2-lgss-varyingT.RData")
} else {
  for (i in 1:length(TT)) {
    
    set.seed(10)
    data <- generateData(c(phi, sigmav, sigmae), TT[i], initialState)
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