##############################################################################
#
# Example of particle Metropolis-Hastings
#
# Copyright (C) 2017 Johan Dahlin < liu (at) johandahlin.com >
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

import numpy as np
from numpy.random import randn, uniform, multivariate_normal
from scipy.stats import gamma, norm

#############################################################################
# Particle Metropolis-Hastings (PMH) for the LGSS model
#
# Inputs:
# y:                   observations from the system for t=1,...,T.
#
# initialPhi:          initial value for phi (persistence of the state)
#
# theta:               the standard deviations of the state noise theta[1]
#                      and observation noise theta[2].
#
# noParticles:         number of particles (N)
#
# initialState:        the initial state.
#
# particleFilter:      function for estimating the likelihood
#
# noIterations:        the number of iterations in PMH 
# 
# stepSize:            the standard deviation of the RW proposal.
#
# Output:
# phi:                 K samples from the parameter posterior for phi.
#
#
#############################################################################


def particleMetropolisHastings(y, initialPhi, theta, noParticles, initialState, particleFilter, noIterations, stepSize):

    # Initalise variables
    phi = np.zeros(noIterations)
    phiProposed = np.zeros(noIterations)
    logLikelihood = np.zeros(noIterations)
    logLikelihoodProposed = np.zeros(noIterations)
    proposedPhiAccepted = np.zeros(noIterations)

    # Set the initial parameter and estimate the initial log-likelihood
    phi[0] = initialPhi
    _, logLikelihood[0] = particleFilter(y, (phi[0], theta[1], theta[2]), noParticles, initialState)

    #=====================================================================
    # Run main loop
    #=====================================================================
    for k in range(1, noIterations):

        # Propose a new parameter
        phiProposed[k] = phi[k - 1] + stepSize * randn()

        # Estimate the log-likelihood if the proposed phi results in a stable model
        if (np.abs(phiProposed[k]) < 1.0):
            _, logLikelihoodProposed[k] = particleFilter(y, (phiProposed[k], theta[1], theta[2]), noParticles, initialState)

        # Compute the acceptance probability
        acceptProbability = np.min((1.0, np.exp(logLikelihoodProposed[k] - logLikelihood[k - 1])))
        
        # Set the acceptance probability to zero if the proposed phi results in an unstable model
        acceptProbability *= np.abs(phiProposed[k]) < 1.0
        
        # Generate uniform random variable in U[0,1]
        uniformRandomVariable = uniform()

        # Accept / reject step
        if uniformRandomVariable < acceptProbability:
            # Accept the parameter
            phi[k] = phiProposed[k]
            logLikelihood[k] = logLikelihoodProposed[k]
            proposedPhiAccepted[k] = 1.0
        else:
            # Reject the parameter
            phi[k] = phi[k - 1]
            logLikelihood[k] = logLikelihood[k - 1]
            proposedPhiAccepted[k] = 0.0

        # Write out progress
        if np.remainder(k, 100) == 0:
            print(
                "##################################################################### ")
            print(" Iteration: " + str(k) +
                  " of : " + str(noIterations) + " completed.")
            print("")
            print(" Current state of the Markov chain:       " + "%.4f" %
                  phi[k] + ".")
            print(" Proposed next state of the Markov chain: " + "%.4f" %
                  phiProposed[k] + ".")
            print(" Current posterior mean:                  " + "%.4f" %
                  np.mean(phi[0:k]) + ".")
            print(" Current acceptance rate:                 " + "%.4f" %
                  np.mean(proposedPhiAccepted[0:k]) + ".")
            print(
                "##################################################################### ")
    
    #=====================================================================
    # Return traces of the parameters
    #=====================================================================
    return phi

#############################################################################
# Particle Metropolis-Hastings (PMH) for the SV model
#
# Inputs:
# y:                   observations from the system for t=1,...,T.
#
# initialTheta:        initial values for the parameters (mu, phi, sigmav)
#
#
# noParticles:         number of particles (N)
#
# initialState:        the initial state.
#
# particleFilter:      function for estimating the likelihood
#
# noIterations:        the number of iterations in PMH
#
# stepSize:            the standard deviation of the RW proposal.
#
# Outputs:
# xHatFiltered:        Estimate of the log-volatility
# theta:               K samples from the parameter posterior.
#
#
#############################################################################


def particleMetropolisHastingsSVModel(y, initialTheta, noParticles, particleFilter, noIterations, stepSize):

    T = len(y)

    # Initalise variables
    theta = np.zeros((noIterations, 3))
    thetaProposed = np.zeros((noIterations, 3))
    logLikelihood = np.zeros(noIterations)
    logLikelihoodProposed = np.zeros(noIterations)
    xHatFiltered = np.zeros((noIterations, T))
    xHatFilteredProposed = np.zeros((noIterations, T))
    proposedThetaAccepted = np.zeros(noIterations)

    # Set the initial parameter and estimate the initial log-likelihood
    theta[0, :] = initialTheta;
    (xHatFiltered[0, :], logLikelihood[0]) = particleFilter(y, theta[0, :], noParticles, T);

    #=====================================================================
    # Run main loop
    #=====================================================================
    for k in range(1, noIterations):

        # Propose a new parameter
        thetaProposed[k, :] = theta[k - 1, :] + \
            multivariate_normal(mean=np.zeros(3), cov=stepSize);

        # Estimate the log-likelihood if the proposed theta results in a stable model
        if ((np.abs(thetaProposed[k, 1]) < 1.0) & (thetaProposed[k, 2] > 0.0)):
            (xHatFilteredProposed[k, :], logLikelihoodProposed[k]) = particleFilter(y, thetaProposed[k, :], noParticles)

        # Compute the ratio between the prior distributions (in log-form)
        prior = norm.logpdf(thetaProposed[k, 0], 0, 1) 
        prior -= norm.logpdf(theta[k - 1, 0], 0, 1)
        
        prior += norm.logpdf(thetaProposed[k, 1], 0.95, 0.05) 
        prior -= norm.logpdf(theta[k - 1, 1], 0.95, 0.05)
        
        prior += gamma.logpdf(thetaProposed[k, 2], 2, 1.0 / 10.0) 
        prior -= gamma.logpdf(theta[k - 1, 2], 2, 1.0 / 10.0)

        # Compute the acceptance probability
        acceptProbability = np.min((1.0, np.exp(prior + logLikelihoodProposed[k] - logLikelihood[k - 1])))

        # Set the acceptance probability to zero if the proposed theta results in an unstable model
        acceptProbability *= np.abs(thetaProposed[k, 1]) < 1.0
        acceptProbability *= thetaProposed[k, 2] > 0.0

        # Generate uniform random variable in U[0,1]
        uniformRandomVariable = uniform()

        # Accept / reject step
        if (uniformRandomVariable < acceptProbability):
            # Accept the parameter
            theta[k, :] = thetaProposed[k, :]
            xHatFiltered[k, :] = xHatFilteredProposed[k, :]
            logLikelihood[k] = logLikelihoodProposed[k]
            proposedThetaAccepted[k] = 1.0
        else:
            # Reject the parameter
            theta[k, :] = theta[k - 1, :]
            xHatFiltered[k, :] = xHatFiltered[k - 1, :]
            logLikelihood[k] = logLikelihood[k - 1]
            proposedThetaAccepted[k] = 0.0

        # Write out progress
        if np.remainder(k, 100) == 0:
            print(
                "##################################################################### ")
            print(" Iteration: " + str(k) +
                  " of : " + str(noIterations) + " completed.")
            print("")
            print(" Current state of the Markov chain:       " + "%.4f" %
                  theta[k, 0] + " " + "%.4f" % theta[k, 1] + " " + "%.4f" % theta[k, 2] + ".")
            print(" Proposed next state of the Markov chain: " + "%.4f" %
                  thetaProposed[k, 0] + " " + "%.4f" % thetaProposed[k, 1] + " " + "%.4f" % thetaProposed[k, 2] + ".")
            print(" Current posterior mean:                  " + "%.4f" % np.mean(
                theta[0:k, 0]) + " " + "%.4f" % np.mean(theta[0:k, 1]) + " " + "%.4f" % np.mean(theta[0:k, 2]) + ".")
            print(" Current acceptance rate:                 " + "%.4f" %
                  np.mean(proposedThetaAccepted[0:k]) + ".")
            print(
                "##################################################################### ")
    
    #=====================================================================
    # Return traces of the parameters
    #=====================================================================
    return (xHatFiltered, theta)
