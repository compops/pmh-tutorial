# Parameter estimation using particle Metropolis-Hastings 
# in a stochastic volatility model

from __future__ import print_function, division
import matplotlib.pylab as plt
import quandl
import numpy as np

from helpers.stateEstimation import particleFilterSVmodel
from helpers.parameterEstimation import particleMetropolisHastingsSVModel


# Set the random seed to replicate results in tutorial
np.random.seed(10)

# Load data
data = quandl.get("NASDAQOMX/OMXS30", trim_start="2012-01-02", trim_end="2014-01-02")
logReturns = 100 * np.diff(np.log(data['Index Value']))
noLogReturns = len(logReturns)

# Settings for PMH
initialTheta = np.array((0.0, 0.9, 0.2))    # Inital guess of theta = (mu, phi, sigmav)
noParticles = 500                           # Choose noParticles ~ noLogReturns
noBurnInIterations = 2500
noIterations = 7500
stepSize = np.diag((0.10**2, 0.01**2, 0.05**2))

# Run the PMH algorithm
logVolatilityEst, parameterTrace = particleMetropolisHastingsSVModel(
    logReturns, initialTheta, noParticles, 
    particleFilterSVmodel, noIterations, stepSize)


# Plot the results
noBins = int(np.floor(np.sqrt(noIterations - noBurnInIterations)))
grid = np.arange(noBurnInIterations, noIterations, 1)
logVolatilityEst = logVolatilityEst[noBurnInIterations:noIterations, :]
parameterEst = parameterTrace[noBurnInIterations:noIterations, :]

plt.figure(1)

plt.subplot(5, 3, (1, 3))
plt.plot(logReturns, color='#1B9E77', linewidth=1.5)
plt.xlabel("time")
plt.ylabel("log-return")

plt.subplot(5, 3, (4, 6))
plt.plot(np.mean(logVolatilityEst, axis=0), color='#D95F02', linewidth=1.5)
plt.xlabel("time")
plt.ylabel("log-volatility estimate")

# Histogram of marginal parameter posterior of mu
plt.subplot(5, 3, 7)
plt.hist(parameterEst[:, 0], noBins, normed=1, facecolor='#7570B3')
plt.xlabel("mu")
plt.ylabel("posterior density estimate")
plt.axvline(np.mean(parameterEst[:, 0]), linewidth=1.5, color='k')

# Trace plot of mu
plt.subplot(5, 3, 8)
plt.plot(grid, parameterEst[:, 0], color='#7570B3')
plt.xlabel("iteration")
plt.ylabel("trace of mu")
plt.axhline(np.mean(parameterEst[:, 0]), linewidth=1.5, color='k')

# Autocorrelation function for mu
plt.subplot(5, 3, 9)
detrended_trace = parameterEst[:, 0] - np.mean(parameterEst[:, 0])
macf = np.correlate(detrended_trace, detrended_trace, mode='full')
idx = int(macf.size/2)
macf = macf[idx:]
macf = macf[0:100]
macf /= macf[0]
grid_acf = range(len(macf))
plt.plot(grid_acf, macf, color='#7570B3')
plt.xlabel("lag")
plt.ylabel("ACF of mu")

# Histogram of marginal parameter posterior of phi
plt.subplot(5, 3, 10)
plt.hist(parameterEst[:, 1], noBins, normed=1, facecolor='#E7298A')
plt.xlabel("phi")
plt.ylabel("posterior density estimate")
plt.axvline(np.mean(parameterEst[:, 1]), linewidth=1.5, color='k')

# Trace plot of phi
plt.subplot(5, 3, 11)
plt.plot(grid, parameterEst[:, 1], color='#E7298A')
plt.xlabel("iteration")
plt.ylabel("trace of phi")
plt.axhline(np.mean(parameterEst[:, 1]), linewidth=1.5, color='k')

# Autocorrelation function for phi
plt.subplot(5, 3, 12)
detrended_trace = parameterEst[:, 1] - np.mean(parameterEst[:, 1])
macf = np.correlate(detrended_trace, detrended_trace, mode='full')
idx = int(macf.size/2)
macf = macf[idx:]
macf = macf[0:100]
macf /= macf[0]
grid_acf = range(len(macf))
plt.plot(grid_acf, macf, color='#E7298A')
plt.xlabel("lag")
plt.ylabel("ACF of phi")

# Histogram of marginal parameter posterior of sigma
plt.subplot(5, 3, 13)
plt.hist(parameterEst[:, 2], noBins, normed=1, facecolor='#66A61E')
plt.xlabel("sigmav")
plt.ylabel("posterior density estimate")
plt.axvline(np.mean(parameterEst[:, 2]), linewidth=1.5, color='k')

# Trace plot of sigma
plt.subplot(5, 3, 14)
plt.plot(grid, parameterEst[:, 2], color='#66A61E')
plt.xlabel("iteration")
plt.ylabel("trace of sigmav")
plt.axhline(np.mean(parameterEst[:, 2]), linewidth=1.5, color='k')

# Autocorrelation function for sigma
plt.subplot(5, 3, 15)
detrended_trace = parameterEst[:, 2] - np.mean(parameterEst[:, 2])
macf = np.correlate(detrended_trace, detrended_trace, mode='full')
idx = int(macf.size/2)
macf = macf[idx:]
macf = macf[0:100]
macf /= macf[0]
grid_acf = range(len(macf))
plt.plot(grid_acf, macf, color='#66A61E')
plt.xlabel("lag")
plt.ylabel("ACF of sigmav")

plt.show()