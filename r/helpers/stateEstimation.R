##############################################################################
# State estimation in LGSS and SV models using Kalman and particle filters.
#
# Johan Dahlin <liu (at) johandahlin.com.nospam>
# Documentation at https://github.com/compops/pmh-tutorial
# Published under GNU General Public License
##############################################################################

##############################################################################
# Fully-adapted particle filter for the linear Gaussian SSM
##############################################################################
particleFilter <- function(y, theta, noParticles, initialState) {
  
  T <- length(y) - 1
  phi <- theta[1] 
  sigmav <- theta[2]
  sigmae <- theta[3]
  
  # Initialise variables
  particles <- matrix(0, nrow = noParticles, ncol = T + 1)
  ancestorIndices <- matrix(0, nrow = noParticles, ncol = T + 1)
  weights <- matrix(1, nrow = noParticles, ncol = T + 1)
  normalisedWeights <- matrix(0, nrow = noParticles, ncol = T + 1)
  xHatFiltered <- matrix(0, nrow = T, ncol = 1)
  logLikelihood <- 0
  
  ancestorIndices[, 1] <- 1:noParticles
  particles[ ,1] <- initialState
  xHatFiltered[ ,1] <- initialState  
  normalisedWeights[, 1] = 1 / noParticles
  
  for (t in 2:T) {
    # Resample ( multinomial )
    newAncestors <- sample(noParticles, replace = TRUE, prob = normalisedWeights[, t - 1])
    ancestorIndices[, 1:(t - 1)] <- ancestorIndices[newAncestors, 1:(t - 1)]
    ancestorIndices[, t] <- newAncestors
    
    # Propagate
    part1 <- (sigmav^(-2) + sigmae^(-2))^(-1)
    part2 <- sigmae^(-2) * y[t]
    part2 <- part2 + sigmav^(-2) * phi * particles[newAncestors, t - 1]
    particles[, t] <- part1 * part2 + rnorm(noParticles, 0, sqrt(part1))
    
    # Compute weights
    yhatMean <- phi * particles[, t]
    yhatVariance <- sqrt(sigmae^2 + sigmav^2)
    weights[, t] <- dnorm(y[t + 1], yhatMean, yhatVariance, log = TRUE)
    
    maxWeight <- max(weights[, t])
    weights[, t] <- exp(weights[, t] - maxWeight)
    
    sumWeights <- sum(weights[, t])
    normalisedWeights[, t] <- weights[, t] / sumWeights
    
    # Estimate the state
    xHatFiltered[t] <- mean(particles[, t])
    
    # Estimate the log-likelihood
    predictiveLikelihood <- maxWeight + log(sumWeights) - log(noParticles)
    logLikelihood <- logLikelihood + predictiveLikelihood
    
  }

  list(xHatFiltered = xHatFiltered,
       logLikelihood = logLikelihood,
       particles = particles,
       weights = normalisedWeights)
  
}

##############################################################################
# Kalman filter for the linear Gaussian SSM
##############################################################################
kalmanFilter <- function(y, theta, initialState, initialStateCovariance) {
  
  T <- length(y)
  yHatPredicted <- matrix(initialState, nrow = T, ncol = 1)
  xHatFiltered <- matrix(initialState, nrow = T, ncol = 1)
  xHatPredicted <- matrix(initialState, nrow = T + 1, ncol = 1)
  predictedStateCovariance <- initialStateCovariance
  logLikelihood <- 0
  
  A <- theta[1] 
  C <- 1
  Q <- theta[2] ^ 2
  R <- theta[3] ^ 2
  
  for (t in 2:T) {
    # Correction step
    S <- C * predictedStateCovariance * C + R
    kalmanGain <- predictedStateCovariance * C / S
    filteredStateCovariance <- predictedStateCovariance - kalmanGain * S * kalmanGain
    
    yHatPredicted[t] <- C * xHatPredicted[t]
    xHatFiltered[t] <- xHatPredicted[t] + kalmanGain * (y[t] - yHatPredicted[t])
    
    # Prediction step
    xHatPredicted[t + 1] <- A * xHatFiltered[t]
    predictedStateCovariance <- A * filteredStateCovariance * A + Q
    
    # Estimate loglikelihood (not in the last iteration, to be able to compare with faPF)
    if (t < T) {
      logLikelihood = logLikelihood + dnorm(y[t], yHatPredicted[t], sqrt(S), log = TRUE)
    }
  }
  
  list(xHatFiltered = xHatFiltered, logLikelihood = logLikelihood)
}

##############################################################################
# Bootstrap particle filter for the stochastic volatility model
##############################################################################
particleFilterSVmodel <- function(y, theta, noParticles) {
  
  T <- length(y) - 1
  mu <- theta[1] 
  phi <- theta[2]
  sigmav <- theta[3]  
  
  particles <- matrix(0, nrow = noParticles, ncol = T + 1)
  ancestorIndices <- matrix(0, nrow = noParticles, ncol = T + 1)
  weights <- matrix(1, nrow = noParticles, ncol = T + 1)
  normalisedWeights <- matrix(0, nrow = noParticles, ncol = T + 1)
  xHatFiltered <- matrix(0, nrow = T, ncol = 1)
  logLikelihood <- 0
  
  ancestorIndices[, 1] <- 1:noParticles
  normalisedWeights[, 1] = 1 / noParticles

  # Generate initial state
  particles[, 1] <- rnorm(noParticles, mu, sigmav / sqrt(1 - phi^2))
  
  for (t in 2:(T + 1)) {
    # Resample ( multinomial )
    newAncestors <- sample(noParticles, replace = TRUE, prob = normalisedWeights[, t - 1])
    ancestorIndices[, 1:(t - 1)] <- ancestorIndices[newAncestors, 1:(t - 1)]
    ancestorIndices[, t] <- newAncestors
    
    # Propagate
    part1 <- mu + phi * (particles[newAncestors, t - 1] - mu)
    particles[, t] <- part1 + rnorm(noParticles, 0, sigmav)
    
    # Compute weights
    yhatMean <- 0
    yhatVariance <- exp(particles[, t] / 2)
    weights[, t] <- dnorm(y[t - 1], yhatMean, yhatVariance, log = TRUE)
    
    maxWeight <- max(weights[, t])
    weights[, t] <- exp(weights[, t] - maxWeight)
    
    sumWeights <- sum(weights[, t])
    normalisedWeights[, t] <- weights[, t] / sumWeights
    
    # Estimate the log-likelihood
    logLikelihood <- logLikelihood + maxWeight + log(sumWeights) - log(noParticles)
    
  }

  # Sample the state estimate using the weights at t=T
  ancestorIndex  <- sample(noParticles, 1, prob = normalisedWeights[, T])
  xHatFiltered <- particles[cbind(ancestorIndices[ancestorIndex, ], 1:(T + 1))]
  
  list(xHatFiltered = xHatFiltered, logLikelihood = logLikelihood)
}