##############################################################################
#
# Example of particle filtering
# in a linear Gaussian state space model
#
# Copyright (C) 2015 Johan Dahlin < johan.dahlin (at) liu.se >
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

# Import some libraries
import matplotlib.pylab as plt
import numpy as np;
import stateEstimationHelper as helpState

# Set the random seed to replicate results in tutorial
np.random.seed( 10 );

##############################################################################
# Define the model
##############################################################################

# Here, we use the following model
#
# x[tt+1] = phi   * x[tt] + sigmav * v[tt]
# y[tt]   = x[tt]         + sigmae * e[tt]
#
# where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

# Set the parameters of the model (par[0],par[1],par[2])
par = np.zeros(3)
par[0] = 0.75
par[1] = 1.00
par[2] = 1.00

# Set the number of time steps to simulate
T      = 250;

# Set the initial state
x0     = 0;

##############################################################################
# Generate data
##############################################################################

(x, y) = helpState.generateData(par, T, x0)

# Plot the measurement
plt.subplot(3,1,1);
plt.plot(y,color = '#1B9E77', linewidth=1.5);
plt.xlabel("time"); plt.ylabel("measurement")

# Plot the states
plt.subplot(3,1,2);
plt.plot(x,color = '#D95F02', linewidth=1.5);
plt.xlabel("time"); plt.ylabel("latent state")

##############################################################################
# State estimation using the particle filter
##############################################################################

plt.subplot(3,1,3);

# Using N = 100 particles and plot the estimate of the latent state
( xhat, ll ) = helpState.pf(y,par,100,T,x0);
plt.plot(xhat,color = '#7570B3', linewidth=1.5)
plt.xlabel("time"); plt.ylabel("state estimate")

###################################################################################
# State estimation using the Kalman filter
###################################################################################

xhat = helpState.kf(y,par,T,x0,0.01);
plt.plot(xhat,'k', linewidth=1.5)

##############################################################################
# End of file
##############################################################################
