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

from __future__ import print_function, division

import numpy as np
from numpy.random import randn, choice
from scipy.stats import norm


##############################################################################
##########################################################################
#
# Kalman filter for the linear Gaussian SSM
#
# Inputs:
# y:                     observations from the system for t=1,...,T.
#
# theta:                 the persistence of the state theta[0] and the
#                        standard deviations of the state noise theta[1] and
#                        observation noise theta[2].
#
# initialState           the initial state 
#
# initialStateCovariance the initial state  covariance
#
# Outputs:
# xHatFiltered:          vector with T elements
#                        estimates of the filtered state
#                        for each t=0,1,...,T-1.
#
# loglikelihood:         estimate of the log-likelihood at T-1
#
##########################################################################
##############################################################################

def kalmanFilter(y, theta, initialState, initialStateCovariance):

    T = len(y)

    # Initalise variables
    predictiveCovariance = initialStateCovariance
    xHatPredicted = initialState * np.ones((T + 1, 1))
    xHatFiltered = initialState * np.ones((T, 1))

    # Set parameters
    A = theta[0]
    C = 1
    Q = theta[1]**2
    R = theta[2]**2

    # Run main loop
    for t in range(0, T):

        # Calculate the Kalman Gain
        S = C * predictiveCovariance * C + R
        kalmanGain = predictiveCovariance * C / S

        # Compute the state estimate
        yHatPredicted = C * xHatPredicted[t]
        xHatFiltered[t] = xHatPredicted[t] + kalmanGain * (y[t - 1] - yHatPredicted)
        xHatPredicted[t + 1] = A * xHatFiltered[t]

        # Update covariance
        filteredCovariance = predictiveCovariance - kalmanGain * S * kalmanGain
        predictiveCovariance = A * filteredCovariance * A + Q

    return xHatFiltered


##############################################################################
##############################################################################
#
# Fully-adapted particle filter for the linear Gaussian SSM
#
# Inputs:
# y:                   observations from the system for t=1,...,T.
#
# theta:               the persistence of the state theta[0] and the
#                      standard deviations of the state noise theta[1] and
#                      observation noise theta[2].
#
# noParticles:         number of particles (N)
#
# initialState:        the initial state
#
# Outputs:
# xHatFiltered:        vector with T elements
#                      estimates of the filtered state
#                      for each t=0,1,...,T-1.
#
# logLikelihood:       estimate of the log-likelihood at T-1
#
##############################################################################
##############################################################################

def particleFilter(y, theta, noParticles, initialState):

    T = len(y)

    # Initalise variables
    particles = np.zeros((noParticles, T))
    ancestorIndices = np.zeros((noParticles, T))
    weights = np.zeros((noParticles, T))
    normalisedWeights = np.zeros((noParticles, T))
    xHatFiltered = np.zeros((T, 1))
    logLikelihood = 0

    # Set the initial state and weight
    particles[:, 0] = initialState
    normalisedWeights[:, 0] = 1.0 / noParticles
    weights[:, 0] = 1.0

    #=====================================================================
    # Run main loop
    #=====================================================================
    for t in range(1, T):

        #=============================================================
        # Resample particles
        #=============================================================
        newAncestors = choice(noParticles, noParticles, p=normalisedWeights[:, t - 1], replace=True)

        # Resample the ancestory linage
        ancestorIndices[:, 1:t - 1] = ancestorIndices[newAncestors, 1:t - 1]

        # Add the most recent ancestors
        ancestorIndices[:, t] = newAncestors

        #=============================================================
        # Propagate particles
        #=============================================================
        particles[:, t] = theta[0] * particles[newAncestors, t - 1] + theta[1] * randn(1, noParticles)

        #=================================================================
        # Weight particles
        #=================================================================

        # Compute log-weights
        weights[:, t] = norm.logpdf(y[t - 1], particles[:, t], theta[2])

        # Rescale log-weights and recover weights
        maxWeight = np.max(weights[:, t])
        weights[:, t] = np.exp(weights[:, t] - maxWeight)

        # Normalise the weights
        sumWeights = np.sum(weights[:, t])
        normalisedWeights[:, t] = weights[:, t] / sumWeights

        # Estimate the filtered state
        xHatFiltered[t] = np.sum(normalisedWeights[:, t] * particles[:, t])

        # Estimate log-likelihood
        predictiveLikelihood = maxWeight + np.log(sumWeights) - np.log(noParticles)
        logLikelihood += predictiveLikelihood

    #=====================================================================
    # Return state estimate and log-likelihood estimate
    #=====================================================================
    return xHatFiltered, logLikelihood


##############################################################################
##############################################################################
#
# Bootstrap particle filter for the stochastic volatility model
#
# Inputs:
# y:                   observations from the system for t=1,...,T.
#
# theta:               the mean of the state theta[0]
#                      the persistence of the state theta[1] and the
#                      standard deviations of the state noise theta[2].
#
# noParticles:         number of particles (N)
#
# Outputs:
# stateTrajectory:     sample from the state filtering equation.
#
# logLikelihood:       estimate of the log-likelihood at T-1
#
##############################################################################
##############################################################################

def particleFilterSVmodel(y, theta, noParticles):

    T = len(y)

    # Initalise variables
    particles = np.zeros((noParticles, T))
    ancestorIndices = np.zeros((noParticles, T))
    weights = np.zeros((noParticles, T))
    normalisedWeights = np.zeros((noParticles, T))
    xHatFiltered = np.zeros((T, 1))
    logLikelihood = 0

    # Set the initial state and weight
    particles[:, 0] = theta[0] + theta[2] / np.sqrt(1.0 - theta[1]**2) * randn(1, noParticles)
    normalisedWeights[:, 0] = 1.0 / noParticles
    weights[:, 0] = 1.0

    #=====================================================================
    # Run main loop
    #=====================================================================
    for t in range(1, T):

        #=============================================================
        # Resample particles
        #=============================================================
        newAncestors = choice(noParticles, noParticles, p=normalisedWeights[:, t - 1], replace=True)

        # Resample the ancestory linage
        ancestorIndices[:, 1:t - 1] = ancestorIndices[newAncestors, 1:t - 1]

        # Add the most recent ancestors
        ancestorIndices[:, t] = newAncestors

        #=============================================================
        # Propagate particles
        #=============================================================
        particles[:, t] = theta[0] + theta[1] * (particles[newAncestors, t - 1] - theta[0]) + theta[2] * randn(1, noParticles)

        #=================================================================
        # Weight particles
        #=================================================================

        # Compute log-weights
        weights[:, t] = norm.logpdf(y[t - 1], 0, np.exp(particles[:, t] / 2))

        # Rescale log-weights and recover weights
        maxWeight = np.max(weights[:, t])
        weights[:, t] = np.exp(weights[:, t] - maxWeight)

        # Normalise the weights
        sumWeights = np.sum(weights[:, t])
        normalisedWeights[:, t] = weights[:, t] / sumWeights

        # Estimate the filtered state
        xHatFiltered[t] = np.sum(normalisedWeights[:, t] * particles[:, t])

        # Estimate log-likelihood
        predictiveLikelihood = maxWeight + np.log(sumWeights) - np.log(noParticles)
        logLikelihood += predictiveLikelihood

    #=====================================================================
    # Return state estimate and log-likelihood estimate
    #=====================================================================

    # Sample the state estimate using the weights at t=T
    ancestorIndex = choice(noParticles, 1, p=normalisedWeights[:, T - 1])
    stateTrajectory = particles[ancestorIndices[ancestorIndex, T - 1].astype(int), :]

    return stateTrajectory, logLikelihood