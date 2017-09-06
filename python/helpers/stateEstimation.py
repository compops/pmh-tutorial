##############################################################################
# State estimation in LGSS and SV models using Kalman and particle filters
# (c) Johan Dahlin 2017 under MIT license <liu@johandahlin.com.nospam>
##############################################################################

from __future__ import print_function, division
import numpy as np
from numpy.random import randn, choice
from scipy.stats import norm

##############################################################################
# Kalman filter for the linear Gaussian SSM
##############################################################################
def kalmanFilter(observations, parameters, initialState, initialStateCov):
    
    noObservations = len(observations)
    A = parameters[0]
    C = 1
    Q = parameters[1]**2
    R = parameters[2]**2

    predictiveCov = initialStateCov
    xHatPredicted = initialState * np.ones((noObservations + 1, 1))
    xHatFiltered = initialState * np.ones((noObservations, 1))

    for t in range(0, noObservations):
        # Correction step
        S = C * predictiveCov * C + R
        kalmanGain = predictiveCov * C / S
        filteredCovariance = predictiveCov - kalmanGain * S * kalmanGain
        yHatPredicted = C * xHatPredicted[t]    
        xHatFiltered[t] = xHatPredicted[t] + kalmanGain * (observations[t - 1] - yHatPredicted)

        # Prediction step
        xHatPredicted[t + 1] = A * xHatFiltered[t]
        predictiveCov = A * filteredCovariance * A + Q

    return xHatFiltered

##############################################################################
# Fully-adapted particle filter for the linear Gaussian SSM
##############################################################################
def particleFilter(observations, parameters, noParticles, initialState):
        
    noObservations = len(observations) - 1
    phi = parameters[0]
    sigmav = parameters[1]
    sigmae = parameters[2]

    particles = np.zeros((noParticles, noObservations))
    ancestorIndices = np.zeros((noParticles, noObservations))
    weights = np.zeros((noParticles, noObservations))
    normalisedWeights = np.zeros((noParticles, noObservations))
    xHatFiltered = np.zeros((noObservations, 1))

    # Set the initial state and weights
    ancestorIndices[: , 0] = range(noParticles)
    particles[:, 0] = initialState
    xHatFiltered[0] = initialState
    normalisedWeights[:, 0] = 1.0 / noParticles
    logLikelihood = 0

    for t in range(1, noObservations):
        # Resample (multinomial)
        newAncestors = choice(noParticles, noParticles, p=normalisedWeights[:, t - 1], replace=True)
        ancestorIndices[:, 1:t - 1] = ancestorIndices[newAncestors, 1:t - 1]
        ancestorIndices[:, t] = newAncestors

        # Propagate
        part1 = (sigmav**(-2) + sigmae**(-2))**(-1)
        part2 = sigmae**(-2) * observations[t]
        part2 = part2 + sigmav**(-2) * phi * particles[newAncestors, t - 1]
        particles[:, t] = part1 * part2 + np.sqrt(part1) * randn(1, noParticles)

        # Compute weights
        yhatMean = phi * particles[:, t]
        yhatVariance = np.sqrt(sigmav**2 + sigmae**2)
        weights[:, t] = norm.logpdf(observations[t + 1], yhatMean, yhatVariance)

        maxWeight = np.max(weights[:, t])
        weights[:, t] = np.exp(weights[:, t] - maxWeight)
        sumWeights = np.sum(weights[:, t])
        normalisedWeights[:, t] = weights[:, t] / sumWeights

        # Estimate the state
        xHatFiltered[t] = np.sum(normalisedWeights[:, t] * particles[:, t])

        # Estimate log-likelihood
        predictiveLikelihood = maxWeight + np.log(sumWeights) - np.log(noParticles)
        logLikelihood += predictiveLikelihood

    return xHatFiltered, logLikelihood

##############################################################################
# Bootstrap particle filter for the stochastic volatility model
##############################################################################
def particleFilterSVmodel(observations, parameters, noParticles):

    noObservations = len(observations)
    mu = parameters[0]
    phi = parameters[1]
    sigmav = parameters[2]

    particles = np.zeros((noParticles, noObservations))
    ancestorIndices = np.zeros((noParticles, noObservations))
    weights = np.zeros((noParticles, noObservations))
    normalisedWeights = np.zeros((noParticles, noObservations))
    xHatFiltered = np.zeros((noObservations, 1))

    # Set the initial state and weights
    particles[:, 0] = mu + sigmav / np.sqrt(1.0 - phi**2) * randn(1, noParticles)
    normalisedWeights[:, 0] = 1.0 / noParticles
    weights[:, 0] = 1.0
    logLikelihood = 0
    
    for t in range(1, noObservations):
        # Resample particles
        newAncestors = choice(noParticles, noParticles, p=normalisedWeights[:, t - 1], replace=True)
        ancestorIndices[:, 1:t - 1] = ancestorIndices[newAncestors, 1:t - 1]
        ancestorIndices[:, t] = newAncestors

        # Propagate particles
        particles[:, t] = mu + phi * (particles[newAncestors, t - 1] - mu) + sigmav * randn(1, noParticles)

        # Weight particles
        weights[:, t] = norm.logpdf(observations[t - 1], 0, np.exp(particles[:, t] / 2))

        maxWeight = np.max(weights[:, t])
        weights[:, t] = np.exp(weights[:, t] - maxWeight)
        sumWeights = np.sum(weights[:, t])
        normalisedWeights[:, t] = weights[:, t] / sumWeights

        # Estimate the filtered state
        xHatFiltered[t] = np.sum(normalisedWeights[:, t] * particles[:, t])

        # Estimate log-likelihood
        predictiveLikelihood = maxWeight + np.log(sumWeights) - np.log(noParticles)
        logLikelihood += predictiveLikelihood

    
    # Sample the state estimate using the weights at t=T
    ancestorIndex = choice(noParticles, 1, p=normalisedWeights[:, noObservations - 1])
    stateTrajectory = particles[ancestorIndices[ancestorIndex, noObservations - 1].astype(int), :]

    return stateTrajectory, logLikelihood