##############################################################################
#
# State estimation in LGSS and SV models using Kalman and particle filters.
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


##############################################################################
# Fully-adapted particle filter (LGSS)
##############################################################################

particleFilter <- function(y, theta, noParticles, initialState) {
  #
  # Fully-adapted particle filter for the linear Gaussian SSM
  #
  # Inputs:
  # y:                   observations from the system for t=1,...,T.
  #
  # theta:               the persistence of the state and the
  # phi, sigmav, sigmae  standard deviations of the state innovations and
  #                      observation noise.
  #
  # noParticles:         number of particles (N)
  #
  # initialState:        initial state
  #
  # Outputs:
  # xHatFiltered:        vector with T elements
  #                      estimates of the filtered state
  #                      for each t=0,1,...,T-1.
  #
  # logLikelihood:       estimate of the log-likelihood at T-1
  #
  #
  
  T <- length(y) - 1
  phi <- theta[1] 
  sigmav <- theta[2]
  sigmae <- theta[3]
  
  #===========================================================
  # Initialise variables
  #===========================================================
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
  
  #===========================================================
  # Run main loop
  #===========================================================
  for (t in 2:T) {
    
    #=========================================================
    # Resample ( multinomial )
    #=========================================================
    newAncestors <- sample(noParticles, replace = TRUE, prob = normalisedWeights[, t - 1])
    ancestorIndices[, 1:(t - 1)] <- ancestorIndices[newAncestors, 1:(t - 1)]
    ancestorIndices[, t] <- newAncestors
    
    #=========================================================
    # Propagate
    #=========================================================
    part1 <- (sigmav^(-2) + sigmae^(-2))^(-1)
    part2 <- sigmae^(-2) * y[t]
    part2 <- part2 + sigmav^(-2) * phi * particles[newAncestors, t - 1]
    particles[, t] <- part1 * part2 + rnorm(noParticles, 0, sqrt(part1))
    
    #=========================================================
    # Compute weights
    #=========================================================
    yhatMean <- phi * particles[, t]
    yhatVariance <- sqrt(sigmae^2 + sigmav^2)
    weights[, t] <- dnorm(y[t + 1], yhatMean, yhatVariance, log = TRUE)
    
    # Rescale log-weights and recover weights
    maxWeight <- max(weights[, t])
    weights[, t] <- exp(weights[, t] - maxWeight)
    
    # Normalize the weights
    sumWeights <- sum(weights[, t])
    normalisedWeights[, t] <- weights[, t] / sumWeights
    
    # Estimate the state
    xHatFiltered[t] <- mean(particles[, t])
    
    # Estimate the log-likelihood
    predictiveLikelihood <- maxWeight + log(sumWeights) - log(noParticles)
    logLikelihood <- logLikelihood + predictiveLikelihood
    
  }
  
  #===========================================================
  # Return state estimate and log-likelihood estimate
  #===========================================================
  list(xHatFiltered = xHatFiltered,
       logLikelihood = logLikelihood,
       particles = particles,
       weights = normalisedWeights)
  
}


###################################################################################
# Kalman filter (LGSS)
###################################################################################
kalmanFilter <- function(y, theta, initialState, initialStateCovariance) {
  #
  # Kalman filter for the linear Gaussian SSM
  #
  # Inputs:
  # y:                      observations from the system for t=1,...,T.
  #
  # theta:                  the persistence of the state and the
  # phi, sigmav, sigmae     standard deviations of the state innovations and
  #                         observation noise.
  #
  # initialStat:            the initial state
  # initialStateCovariance: the covariance of the initial state
  #
  # Outputs:
  # xHatFiltered:           vector with T elements
  #                         estimates of the filtered state
  #                         for each t=0,1,...,T-1.
  #
  # logLikelihood:         estimate of the log-likelihood at T-1
  #
  #
  
  T <- length(y)
  yHatPredicted <- matrix(initialState, nrow = T, ncol = 1)
  xHatFiltered <- matrix(initialState, nrow = T, ncol = 1)
  xHatPredicted <- matrix(initialState, nrow = T + 1, ncol = 1)
  predictedStateCovariance <- initialStateCovariance
  logLikelihood <- 0
  
  # Set parameters
  A <- theta[1] 
  C <- 1
  Q <- theta[2] ^ 2
  R <- theta[3] ^ 2
  
  for (t in 2:T) {
    # Compute Kalman Gain
    S <- C * predictedStateCovariance * C + R
    kalmanGain <- predictedStateCovariance * C / S
    
    # Compute state estimate
    yHatPredicted[t] <- C * xHatPredicted[t]
    xHatFiltered[t] <- xHatPredicted[t] + kalmanGain * (y[t] - yHatPredicted[t])
    xHatPredicted[t + 1] <- A * xHatFiltered[t]
    
    # Update covariance
    filteredStateCovariance <- predictedStateCovariance - kalmanGain * S * kalmanGain
    predictedStateCovariance <- A * filteredStateCovariance * A + Q
    
    # Estimate loglikelihood (not in the last iteration, to be able to compare with faPF)
    if (t < T) {
      logLikelihood = logLikelihood + dnorm(y[t], yHatPredicted[t], sqrt(S), log = TRUE)
    }
  }
  
  list(xHatFiltered = xHatFiltered, logLikelihood = logLikelihood)
}


##############################################################################
# Bootstrap particle filter (SV model)
##############################################################################
particleFilterSVmodel <- function(y, theta, noParticles) {
  #
  # Bootstrap particle filter for the stochastic volatility model
  #
  # Inputs:
  # y:                   observations from the system for t=1,...,T.
  #
  # theta:               the mean and persistence of the state and the
  # mu, phi, sigmav      standard deviations of the state innovations.
  #
  # noParticles:         number of particles (N)
  #
  # Outputs:
  # xHatFiltered:        vector with T+1 elements
  #                      estimates of the smoothed state
  #                      for each t=0,1,...,T.
  #
  # logLikelihood:       estimate of the log-likelihood at T
  #
  #

  T <- length(y) - 1
  mu <- theta[1] 
  phi <- theta[2]
  sigmav <- theta[3]  
    
  #===========================================================
  # Initialise variables
  #===========================================================
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
  
  #===========================================================
  # Run main loop
  #===========================================================
  for (t in 2:(T + 1)) {
    #=========================================================
    # Resample ( multinomial )
    #=========================================================
    newAncestors <- sample(noParticles, replace = TRUE, prob = normalisedWeights[, t - 1])
    
    # Resample the ancestory linage
    ancestorIndices[, 1:(t - 1)] <- ancestorIndices[newAncestors, 1:(t - 1)]
    
    # Add the most recent ancestors
    ancestorIndices[, t] <- newAncestors
    
    #=========================================================
    # Propagate
    #=========================================================
    part1 <- mu + phi * (particles[newAncestors, t - 1] - mu)
    particles[, t] <- part1 + rnorm(noParticles, 0, sigmav)
    
    #=========================================================
    # Compute weights
    #=========================================================
    yhatMean <- 0
    yhatVariance <- exp(particles[, t] / 2)
    weights[, t] <- dnorm(y[t - 1], yhatMean, yhatVariance, log = TRUE)
    
    # Rescale log-weights and recover weights
    maxWeight <- max(weights[, t])
    weights[, t] <- exp(weights[, t] - maxWeight)
    
    # Normalize the weights
    sumWeights <- sum(weights[, t])
    normalisedWeights[, t] <- weights[, t] / sumWeights
    
    # Estimate the log-likelihood
    logLikelihood <- logLikelihood + maxWeight + log(sumWeights) - log(noParticles)
    
  }
  #===========================================================
  # Return state estimate and log-likelihood estimate
  #===========================================================
  
  # Sample the state estimate using the weights at t=T
  ancestorIndex  <- sample(noParticles, 1, prob = normalisedWeights[, T])
  xHatFiltered <- particles[cbind(ancestorIndices[ancestorIndex, ], 1:(T + 1))]
  
  list(xHatFiltered = xHatFiltered, logLikelihood = logLikelihood)
}


##############################################################################
# End of file
##############################################################################