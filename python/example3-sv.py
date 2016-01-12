##############################################################################
#
# Example of particle Metropolis-Hastings in a stochastic volatility model
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
import matplotlib.pylab          as plt
import Quandl
import numpy                     as np
import stateEstimationHelper     as helpState
import parameterEstimationHelper as helpParam

# Set the random seed to replicate results in tutorial
np.random.seed( 10 );

##############################################################################
# Define the model
##############################################################################

# Here, we use the following model
#
# x[tt+1] = par[0] * x[tt] + par[1]  * v[tt]
# y[tt]   = par[2] * exp( xt[tt]/2 ) * e[tt]
#
# where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

# Set the number of time steps
T      = 500;

##############################################################################
# Load data
##############################################################################
d = Quandl.get("NASDAQOMX/OMXS30", trim_start="2012-01-02", trim_end="2014-01-02")
y = 100 * np.diff(np.log(d['Index Value']))

##############################################################################
# Parameter estimation using PMH
##############################################################################

# The inital guess of the parameter (mu, phi, sigmav)
initPar  = np.array(( 0.0, 0.9, 0.2));

# No. particles in the particle filter ( choose nPart ~ T )
nPart    = 500;

# The length of the burn-in and the no. iterations of PMH ( nBurnIn < nRuns )
nBurnIn  = 2500;
nRuns    = 7500;

# The standard deviation in the random walk proposal
stepSize = np.diag( ( 0.01**2, 0.05**2, 0.05**2) );

# Run the PMH algorithm
( xhat, thhat) = helpParam.pmh_sv(y,initPar,nPart,T,helpState.pf_sv,nRuns,stepSize)

#=============================================================================
# Plot the results
#=============================================================================

plt.figure(1);

# Plot the log-returns with volatility estimate
plt.subplot(4,2,(1,2));
plt.plot(y,color = '#1B9E77', linewidth=1.5);
plt.plot(np.mean(xhat[nBurnIn:nRuns,:],axis=0),color = '#D95F02', linewidth=1.5);
plt.xlabel("time"); plt.ylabel("log-return")

# Plot the parameter posterior estimate
# Plot the trace of the Markov chain after burn-in
# Solid black line indicate posterior mean
plt.subplot(4,2,3);
plt.hist(thhat[nBurnIn:nRuns,0], np.floor(np.sqrt(nRuns-nBurnIn)), normed=1, facecolor='#7570B3')
plt.xlabel("mu"); plt.ylabel("posterior density estimate")
plt.axvline( np.mean(thhat[nBurnIn:nRuns,0]), linewidth=1.5, color = 'k' )

plt.subplot(4,2,4);
plt.plot(np.arange(nBurnIn,nRuns,1),thhat[nBurnIn:nRuns,0],color='#7570B3')
plt.xlabel("iteration"); plt.ylabel("trace of mu")
plt.axhline( np.mean(thhat[nBurnIn:nRuns,0]), linewidth=1.5, color = 'k' )

plt.subplot(4,2,5);
plt.hist(thhat[nBurnIn:nRuns,1], np.floor(np.sqrt(nRuns-nBurnIn)), normed=1, facecolor='#E7298A')
plt.xlabel("phi"); plt.ylabel("posterior density estimate")
plt.axvline( np.mean(thhat[nBurnIn:nRuns,1]), linewidth=1.5, color = 'k' )

plt.subplot(4,2,6);
plt.plot(np.arange(nBurnIn,nRuns,1),thhat[nBurnIn:nRuns,1],color='#E7298A')
plt.xlabel("iteration"); plt.ylabel("trace of phi")
plt.axhline( np.mean(thhat[nBurnIn:nRuns,1]), linewidth=1.5, color = 'k' )

plt.subplot(4,2,7);
plt.hist(thhat[nBurnIn:nRuns,2], np.floor(np.sqrt(nRuns-nBurnIn)), normed=1, facecolor='#66A61E')
plt.xlabel("sigmav"); plt.ylabel("posterior density estimate")
plt.axvline( np.mean(thhat[nBurnIn:nRuns,2]), linewidth=1.5, color = 'k' )

plt.subplot(4,2,8);
plt.plot(np.arange(nBurnIn,nRuns,1),thhat[nBurnIn:nRuns,2],color='#66A61E')
plt.xlabel("iteration"); plt.ylabel("trace of sigmav")
plt.axhline( np.mean(thhat[nBurnIn:nRuns,2]), linewidth=1.5, color = 'k' )

##############################################################################
# End of file
##############################################################################
