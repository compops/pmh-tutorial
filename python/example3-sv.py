##############################################################################
#
# Example of particle Metropolis-Hastings in a stochastic volatility model
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

# Import libraries
import matplotlib.pylab as plt
import quandl
import numpy as np

# Import helpers
from helpers.stateEstimation import particleFilterSVmodel
from helpers.parameterEstimation import particleMetropolisHastingsSVModel

# Set the random seed to replicate results in tutorial
np.random.seed(10)


#=============================================================================
# Load data
#=============================================================================
d = quandl.get("NASDAQOMX/OMXS30", trim_start="2012-01-02", trim_end="2014-01-02")
y = 100 * np.diff(np.log(d['Index Value']))
T = len(y)


#=============================================================================
# Parameter estimation using PMH
#=============================================================================

# The inital guess of the parameter (mu, phi, sigmav)
initialTheta = np.array((0.0, 0.9, 0.2))

# No. particles in the particle filter ( choose noParticles ~ T )
noParticles = 500

# The length of the burn-in and the no. iterations of PMH ( nBurnIn < nRuns )
noBurnInIterations = 2500
noIterations = 7500

# The standard deviation in the random walk proposal
stepSize = np.diag((0.01**2, 0.05**2, 0.05**2))

# Run the PMH algorithm
xHatFiltered, theta = particleMetropolisHastingsSVModel(y, initialTheta, noParticles, particleFilterSVmodel, noIterations, stepSize)


#=============================================================================
# Plot the results
#=============================================================================
noBins = int(np.floor(np.sqrt(noIterations - noBurnInIterations)))
grid = np.arange(noBurnInIterations, noIterations, 1)
resXhat = xHatFiltered[noBurnInIterations:noIterations, :]
resTheta = theta[noBurnInIterations:noIterations, :]

plt.figure(1)

## Plot the log-returns with volatility estimate
plt.subplot(4, 2, (1, 2))
plt.plot(y, color='#1B9E77', linewidth=1.5)
plt.plot(np.mean(resXhat, axis=0), color='#D95F02', linewidth=1.5)
plt.xlabel("time")
plt.ylabel("log-return")

## Plot the parameter posterior estimate (trace after burn-in and means)

# Mu
plt.subplot(4, 2, 3)
plt.hist(resTheta[:, 0], noBins, normed=1, facecolor='#7570B3')
plt.xlabel("mu")
plt.ylabel("posterior density estimate")
plt.axvline(np.mean(resTheta[:, 0]), linewidth=1.5, color='k')

plt.subplot(4, 2, 4)
plt.plot(grid, resTheta[:, 0], color='#7570B3')
plt.xlabel("iteration")
plt.ylabel("trace of mu")
plt.axhline(np.mean(resTheta[:, 0]), linewidth=1.5, color='k')

# Phi
plt.subplot(4, 2, 5)
plt.hist(resTheta[:, 1], noBins, normed=1, facecolor='#E7298A')
plt.xlabel("phi")
plt.ylabel("posterior density estimate")
plt.axvline(np.mean(resTheta[:, 1]), linewidth=1.5, color='k')

plt.subplot(4, 2, 6)
plt.plot(grid, resTheta[:, 1], color='#E7298A')
plt.xlabel("iteration")
plt.ylabel("trace of phi")
plt.axhline(np.mean(resTheta[:, 1]), linewidth=1.5, color='k')

# Sigma
plt.subplot(4, 2, 7)
plt.hist(resTheta[:, 2], noBins, normed=1, facecolor='#66A61E')
plt.xlabel("sigmav")
plt.ylabel("posterior density estimate")
plt.axvline(np.mean(resTheta[:, 2]), linewidth=1.5, color='k')

plt.subplot(4, 2, 8)
plt.plot(grid, resTheta[:, 2], color='#66A61E')
plt.xlabel("iteration")
plt.ylabel("trace of sigmav")
plt.axhline(np.mean(resTheta[:, 2]), linewidth=1.5, color='k')


##############################################################################
# End of file
##############################################################################
