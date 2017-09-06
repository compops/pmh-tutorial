##############################################################################
# Parameter estimation using particle Metropolis-Hastings in a LGSS model.
#
# Johan Dahlin <liu (at) johandahlin.com.nospam>
# Documentation at https://github.com/compops/pmh-tutorial
# Published under GNU General Public License
##############################################################################

from __future__ import print_function, division
import matplotlib.pylab as plt
import numpy as np

from helpers.dataGeneration import generateData
from helpers.stateEstimation import particleFilter, kalmanFilter
from helpers.parameterEstimation import particleMetropolisHastings

# Set the random seed to replicate results in tutorial
np.random.seed(10)

##############################################################################
# Define the model and generate data
# x[t + 1] = phi * x[t] + sigmav * v[t],    v[t] ~ N(0, 1)
# y[t] = x[t] + sigmae * e[t],              e[t] ~ N(0, 1)
##############################################################################
parameters = np.zeros(3)    # theta = (phi, sigmav, sigmae)
parameters[0] = 0.75
parameters[1] = 1.00
parameters[2] = 0.10
noObservations = 250
initialState = 0

state, observations = generateData(parameters, noObservations, initialState)

##############################################################################
# PMH
##############################################################################
initialPhi = 0.50
noParticles = 500           # Use noParticles ~ noObservations
noBurnInIterations = 1000
noIterations = 5000
stepSize = 0.10

phiTrace = particleMetropolisHastings(
    observations, initialPhi, parameters, noParticles, 
    initialState, particleFilter, noIterations, stepSize)

##############################################################################
# Plot the results
##############################################################################
noBins = int(np.floor(np.sqrt(noIterations - noBurnInIterations)))
grid = np.arange(noBurnInIterations, noIterations, 1)
phiTrace = phiTrace[noBurnInIterations:noIterations]

# Plot the parameter posterior estimate (solid black line = posterior mean)
plt.subplot(3, 1, 1)
plt.hist(phiTrace, noBins, normed=1, facecolor='#7570B3')
plt.xlabel("phi")
plt.ylabel("posterior density estimate")
plt.axvline(np.mean(phiTrace), color='k')

# Plot the trace of the Markov chain after burn-in (solid black line = posterior mean)
plt.subplot(3, 1, 2)
plt.plot(grid, phiTrace, color='#7570B3')
plt.xlabel("iteration")
plt.ylabel("phi")
plt.axhline(np.mean(phiTrace), color='k')

# Plot the autocorrelation function
plt.subplot(3, 1, 3)
macf = np.correlate(phiTrace - np.mean(phiTrace), phiTrace - np.mean(phiTrace), mode='full')
idx = int(macf.size/2)
macf = macf[idx:]
macf = macf[0:100]
macf /= macf[0]
grid = range(len(macf))
plt.plot(grid, macf, color='#7570B3')
plt.xlabel("lag")
plt.ylabel("ACF of phi")

plt.show()