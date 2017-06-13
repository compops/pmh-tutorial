##############################################################################
#
# Example of particle Metropolis-Hastings in a LGSS model.
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

# Import libraries
import matplotlib.pylab as plt
import numpy as np

# Import helpers
from helpers.dataGeneration import generateData
from helpers.stateEstimation import particleFilter, kalmanFilter
from helpers.parameterEstimation import particleMetropolisHastings

# Set the random seed
np.random.seed(10)


#=============================================================================
# Define the model
#=============================================================================

# Here, we use the following model
#
# x[t + 1] = phi * x[t] + sigmav * v[t]
# y[t] = x[t] + sigmae * e[t]
#
# where v[t] ~ N(0, 1) and e[t] ~ N(0, 1)

# Set the parameters of the model (phi, sigmav, sigmae)
theta = np.zeros(3)
theta[0] = 0.75
theta[1] = 1.00
theta[2] = 0.10

# Set the number of time steps to simulate
T = 250

# Set the initial state
initialState = 0


#=============================================================================
# Generate data
#=============================================================================

x, y = generateData(theta, T, initialState)


#=============================================================================
# Parameter estimation using PMH
#=============================================================================

# The inital guess of the parameter
initialPhi = 0.50

# No. particles in the particle filter ( choose noParticles ~ T )
noParticles = 500

# The length of the burn-in and the no. iterations of PMH 
# ( noBurnInIterations < noIterations )
noBurnInIterations = 1000
noIterations = 5000

# The standard deviation in the random walk proposal
stepSize = 0.10

# Run the PMH algorithm
res = particleMetropolisHastings(y, initialPhi, theta, noParticles, initialState,
                                 particleFilter, noIterations, stepSize)


#=============================================================================
# Plot the results
#=============================================================================

noBins = int(np.floor(np.sqrt(noIterations - noBurnInIterations)))
grid = np.arange(noBurnInIterations, noIterations, 1)
resPhi = res[noBurnInIterations:noIterations]

# Plot the parameter posterior estimate
# Solid black line indicate posterior mean
plt.subplot(3, 1, 1)
plt.hist(resPhi, noBins, normed=1, facecolor='#7570B3')
plt.xlabel("phi")
plt.ylabel("posterior density estimate")
plt.axvline(np.mean(resPhi), color='k')

# Plot the trace of the Markov chain after burn-in
# Solid black line indicate posterior mean
plt.subplot(3, 1, 2)
plt.plot(grid, resPhi, color='#7570B3')
plt.xlabel("iteration")
plt.ylabel("phi")
plt.axhline(np.mean(resPhi), color='k')

# Plot the autocorrelation function
plt.subplot(3, 1, 3)
macf = np.correlate(resPhi - np.mean(resPhi), resPhi - np.mean(resPhi), mode='full')
macf = macf[macf.size/2:]
macf = macf[0:100]
macf /= macf[0]
grid = range(len(macf))
plt.plot(grid, macf, color='#7570B3')
plt.xlabel("lag")
plt.ylabel("ACF of phi")


##############################################################################
# End of file
##############################################################################
