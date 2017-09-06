##############################################################################
# Particle Metropolis-Hastings for LGSS and SV models
# (c) Johan Dahlin 2017 under MIT license <liu@johandahlin.com.nospam>
##############################################################################

from __future__ import print_function, division
import numpy as np
from numpy.random import randn, uniform, multivariate_normal
from scipy.stats import gamma, norm

##############################################################################
# Particle Metropolis-Hastings (PMH) for the LGSS model
##############################################################################
def particleMetropolisHastings(observations, initialPhi, parameters, noParticles, 
        initialState, particleFilter, noIterations, stepSize):

    phi = np.zeros(noIterations)
    phiProposed = np.zeros(noIterations)
    logLikelihood = np.zeros(noIterations)
    logLikelihoodProposed = np.zeros(noIterations)
    proposedPhiAccepted = np.zeros(noIterations)

    # Set the initial parameter and estimate the initial log-likelihood
    phi[0] = initialPhi
    _, logLikelihood[0] = particleFilter(observations, (phi[0], parameters[1], parameters[2]), noParticles, initialState)

    for k in range(1, noIterations):
        # Propose a new parameter
        phiProposed[k] = phi[k - 1] + stepSize * randn()

        # Estimate the log-likelihood if the proposed phi results in a stable model
        if (np.abs(phiProposed[k]) < 1.0):
            _, logLikelihoodProposed[k] = particleFilter(observations, (phiProposed[k], parameters[1], parameters[2]), noParticles, initialState)

        # Compute the acceptance probability
        acceptProbability = np.min((1.0, np.exp(logLikelihoodProposed[k] - logLikelihood[k - 1])))
        acceptProbability *= np.abs(phiProposed[k]) < 1.0
        
        # Accept / reject step
        uniformRandomVariable = uniform()
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
            print("#####################################################################")
            print(" Iteration: " + str(k) + " of : " + str(noIterations) + " completed.")
            print("")
            print(" Current state of the Markov chain:       " + "%.4f" % phi[k] + ".")
            print(" Proposed next state of the Markov chain: " + "%.4f" % phiProposed[k] + ".")
            print(" Current posterior mean:                  " + "%.4f" % np.mean(phi[0:k]) + ".")
            print(" Current acceptance rate:                 " + "%.4f" % np.mean(proposedPhiAccepted[0:k]) + ".")
            print("#####################################################################")
    
    return phi

##############################################################################
# Particle Metropolis-Hastings (PMH) for the SV model
##############################################################################
def particleMetropolisHastingsSVModel(observations, initialTheta, 
        noParticles, particleFilter, noIterations, stepSize):

    noObservations = len(observations)

    theta = np.zeros((noIterations, 3))
    thetaProposed = np.zeros((noIterations, 3))
    logLikelihood = np.zeros(noIterations)
    logLikelihoodProposed = np.zeros(noIterations)
    xHatFiltered = np.zeros((noIterations, noObservations))
    xHatFilteredProposed = np.zeros((noIterations, noObservations))
    proposedThetaAccepted = np.zeros(noIterations)

    # Set the initial parameter and estimate the initial log-likelihood
    theta[0, :] = initialTheta
    (xHatFiltered[0, :], logLikelihood[0]) = particleFilter(observations, theta[0, :], noParticles)

    for k in range(1, noIterations):

        # Propose a new parameter
        thetaProposed[k, :] = theta[k - 1, :] + multivariate_normal(mean = np.zeros(3), cov = stepSize)

        # Estimate the log-likelihood if the proposed theta results in a stable model
        if ((np.abs(thetaProposed[k, 1]) < 1.0) & (thetaProposed[k, 2] > 0.0)):
            (xHatFilteredProposed[k, :], logLikelihoodProposed[k]) = particleFilter(observations, thetaProposed[k, :], noParticles)

        # Compute the ratio between the prior distributions (in log-form)
        prior = norm.logpdf(thetaProposed[k, 0], 0, 1) 
        prior -= norm.logpdf(theta[k - 1, 0], 0, 1)

        prior += norm.logpdf(thetaProposed[k, 1], 0.95, 0.05) 
        prior -= norm.logpdf(theta[k - 1, 1], 0.95, 0.05)

        prior += gamma.logpdf(thetaProposed[k, 2], 2, 1.0 / 10.0) 
        prior -= gamma.logpdf(theta[k - 1, 2], 2, 1.0 / 10.0)

        # Compute the acceptance probability
        acceptProbability = np.min((1.0, np.exp(prior + logLikelihoodProposed[k] - logLikelihood[k - 1])))
        acceptProbability *= np.abs(thetaProposed[k, 1]) < 1.0
        acceptProbability *= thetaProposed[k, 2] > 0.0

        # Accept / reject step
        uniformRandomVariable = uniform()
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
            print("#####################################################################")
            print(" Iteration: " + str(k) + " of : " + str(noIterations) + " completed.")
            print("")
            print(" Current state of the Markov chain:       " + "%.4f" % theta[k, 0] + " " + "%.4f" % theta[k, 1] + " " + "%.4f" % theta[k, 2] + ".")
            print(" Proposed next state of the Markov chain: " + "%.4f" % thetaProposed[k, 0] + " " + "%.4f" % thetaProposed[k, 1] + " " + "%.4f" % thetaProposed[k, 2] + ".")
            print(" Current posterior mean:                  " + "%.4f" % np.mean( theta[0:k, 0]) + " " + "%.4f" % np.mean(theta[0:k, 1]) + " " + "%.4f" % np.mean(theta[0:k, 2]) + ".")
            print(" Current acceptance rate:                 " + "%.4f" % np.mean(proposedThetaAccepted[0:k]) + ".")
            print("#####################################################################")
    
    return (xHatFiltered, theta)
