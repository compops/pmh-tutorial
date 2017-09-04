##############################################################################
#
# Example of state estimation in a LGSS model 
# using particle filters and Kalman filters
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

# Set the random seed to replicate results in tutorial
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

# Plot the measurement
plt.subplot(3, 1, 1)
plt.plot(y, color='#1B9E77', linewidth=1.5)
plt.xlabel("time")
plt.ylabel("measurement")

# Plot the states
plt.subplot(3, 1, 2)
plt.plot(x, color='#D95F02', linewidth=1.5)
plt.xlabel("time")
plt.ylabel("latent state")


#=============================================================================
# State estimation
#=============================================================================

# Using N = 100 particles and plot the estimate of the latent state
xHatFilteredParticleFilter, _ = particleFilter(y, theta, 100, initialState)

# Using the Kalman filter
xHatFilteredKalmanFilter = kalmanFilter(y, theta, initialState, 0.01)

plt.subplot(3, 1, 3)
plt.plot(xHatFilteredKalmanFilter[1:T] - xHatFilteredParticleFilter[0:T-1], color='#7570B3', linewidth=1.5)
plt.xlabel("time")
plt.ylabel("difference in estimate")
plt.show()

##############################################################################
# End of file
##############################################################################
