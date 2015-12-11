##############################################################################
#
# Example of particle Metropolis-Hastings
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

import numpy as np
from scipy.stats import gamma
from scipy.stats import norm

#############################################################################
# Particle Metropolis-Hastings (PMH) for the LGSS model
#
# Inputs:
# y:                   observations from the system for t=1,...,T.
#
# initPar:             initial value for phi (persistence of the state)
#
# par:                 the standard deviations of the state noise par[1]
#                      and observation noise par[2].
#
# nPart:               number of particles (N)
#
# T and xo:            the no. observations and initial state.
#
# sm:                  function for estimating the likelihood
#
# nIter and stepSize:  the number of iterations in PMH and the
#                      standard deviation of the RW proposal.
#
# Outputs:
# th:                  K samples from the parameter posterior.
#
#
#############################################################################

def pmh(y,initPar,par,nPart,T,xo,sm,nIter,stepSize):

    # Initalise variables
    th     = np.zeros(nIter);
    thp    = np.zeros(nIter);
    ll     = np.zeros(nIter);
    llp    = np.zeros(nIter);
    accept = np.zeros(nIter);

    # Set the initial parameter and estimate the initial log-likelihood
    th[0]  = initPar;
    ( xhat, ll[0] ) = sm(y,(th[0],par[1],par[2]),nPart,T,xo);

    #=====================================================================
    # Run main loop
    #=====================================================================
    for kk in range(1, nIter):

        # Propose a new parameter
        thp[kk] = th[kk-1] + stepSize * np.random.randn();

        # Check if | par[0] | > 1.0, in that case reject the parameter.
        # This is due to that the model
        # is only stable when | par[0] | is smaller than 1
        if ( np.abs( thp[kk] ) < 1.0 ):

            # Estimate the log-likelihood
            ( xhat, llp[kk] ) = sm(y,(thp[kk],par[1],par[2]),nPart,T,xo);

            # Compute the acceptance probability
            aprob = np.min( (1.0, np.exp( llp[kk] - ll[kk-1] ) ) );

            # Generate uniform random variable in U[0,1]
            u = np.random.uniform()

        # Accept / reject step
        if ( u < aprob ):
            # Accept the parameter
            th[kk]     = thp[kk]
            ll[kk]     = llp[kk]
            accept[kk] = 1.0;
        else:
            # Reject the parameter
            th[kk]     = th[kk-1]
            ll[kk]     = ll[kk-1]
            accept[kk] = 0.0;

        # Write out progress
        if np.remainder(kk,100) == 0:
            print("##################################################################### ");
            print(" Iteration: " + str(kk) + " of : " + str(nIter) + " completed.")
            print("");
            print(" Current state of the Markov chain:       " + "%.4f" % th[kk] + "." )
            print(" Proposed next state of the Markov chain: " + "%.4f" % thp[kk] + "." )
            print(" Current posterior mean:                  " + "%.4f" % np.mean(th[0:kk]) + "." )
            print(" Current acceptance rate:                 " + "%.4f" % np.mean(accept[0:kk]) + "." )
            print("##################################################################### ");
    #=====================================================================
    # Return traces of the parameters
    #=====================================================================
    return th;

#############################################################################
# Particle Metropolis-Hastings (PMH) for the SV model
#
# Inputs:
# y:                   observations from the system for t=1,...,T.
#
# initPar:             initial values for the parameters (mu,phi,sigmav)
#
#
# nPart:               number of particles (N)
#
# T and xo:            the no. observations and initial state.
#
# sm:                  function for estimating the likelihood
#
# nIter and stepSize:  the number of iterations in PMH and the
#                      standard deviation of the RW proposal.
#
# Outputs:
# x:                   Estimate of the log-volatility
# th:                  K samples from the parameter posterior.
#
#
#############################################################################

def pmh_sv(y,initPar,nPart,T,sm,nIter,stepSize):

    # Initalise variables
    th     = np.zeros((nIter,3));
    thp    = np.zeros((nIter,3));
    ll     = np.zeros(nIter);
    llp    = np.zeros(nIter);
    x      = np.zeros((nIter,T));
    xp     = np.zeros((nIter,T));
    accept = np.zeros(nIter);

    # Set the initial parameter and estimate the initial log-likelihood
    th[0,:]         = initPar;
    ( x[0,:], ll[0] ) = sm(y,th[0,:],nPart,T);

    #=====================================================================
    # Run main loop
    #=====================================================================
    for kk in range(1, nIter):

        # Propose a new parameter
        thp[kk,:] = th[kk-1,:] + np.random.multivariate_normal( mean=np.zeros(3), cov=stepSize );

        # Estimate the log-likelihood
        # Dont run if system is unstable
        if ( ( np.abs( thp[kk,1] ) < 1.0 ) & ( thp[kk,2] > 0.0 ) ):
            ( xp[kk,:], llp[kk] ) = sm(y,thp[kk,:],nPart,T);

        # Compute the ratio between the prior distributions (in log-form)
        prior =  norm.logpdf(  thp[kk,0], 0, 1)       - norm.logpdf(  th[kk-1,0], 0, 1);
        prior += norm.logpdf(  thp[kk,1], 0.95, 0.05) - norm.logpdf(  th[kk-1,1], 0.95, 0.05)
        prior += gamma.logpdf( thp[kk,2], 2, 1.0/10.0)- gamma.logpdf( th[kk-1,2], 2, 1.0/10.0);

        # Compute the acceptance probability
        aprob = np.min( (1.0, np.exp( prior + llp[kk] - ll[kk-1] ) ) );

        # Generate uniform random variable in U[0,1]
        u = np.random.uniform()

        # Check if | par[0] | > 1.0, in that case set u = 1.0;
        # Check if par[1] < 0.0, in that case set u = 1.0;
        # Check if par[2] < 0.0, in that case set u = 1.0;
        # So that the parameter is rejected. This is due to that the model
        # is only stable when | par[0] | is smaller than 1
        if ( ( np.abs( thp[kk,1] ) > 1.0 ) | ( thp[kk,2] < 0.0 ) ):
            u = 1.0;

        # Accept / reject step
        if ( u < aprob ):
            # Accept the parameter
            th[kk,:]   = thp[kk,:]
            x[kk,:]    = xp[kk,:]
            ll[kk]     = llp[kk]
            accept[kk] = 1.0;
        else:
            # Reject the parameter
            th[kk,:]   = th[kk-1,:]
            x[kk,:]    = x[kk-1,:]
            ll[kk]     = ll[kk-1]
            accept[kk] = 0.0;

        # Write out progress
        if np.remainder(kk,100) == 0:
            print("##################################################################### ");
            print(" Iteration: " + str(kk) + " of : " + str(nIter) + " completed.")
            print("");
            print(" Current state of the Markov chain:       " + "%.4f" % th[kk,0] + " " + "%.4f" % th[kk,1] + " " + "%.4f" % th[kk,2] + "." )
            print(" Proposed next state of the Markov chain: " + "%.4f" % thp[kk,0] + " " + "%.4f" % thp[kk,1] + " " + "%.4f" % thp[kk,2] + "." )
            print(" Current posterior mean:                  " + "%.4f" % np.mean(th[0:kk,0]) + " " + "%.4f" % np.mean(th[0:kk,1]) + " " + "%.4f" % np.mean(th[0:kk,2]) + "." )
            print(" Current acceptance rate:                 " + "%.4f" % np.mean(accept[0:kk]) + "." )
            print("##################################################################### ");
    #=====================================================================
    # Return traces of the parameters
    #=====================================================================
    return ( x, th );
