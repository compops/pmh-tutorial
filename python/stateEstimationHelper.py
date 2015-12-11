##############################################################################
#
# Example of particle filtering
#
# Subroutine for data generation and particle filtering
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
from   scipy.stats import norm


##############################################################################
#
# Generates data from the LGSS model with parameters (phi,sigmav,sigmae)
#
# Inputs:
# par:                 the persistence of the state par[0] and the
#                      standard deviations of the state noise par[1] and
#                      observation noise par[2].
#
# T and x0:            the no. observations and initial state.
#
# Outputs:
# x,y:                 the latent state and observations
#
#
##############################################################################

def generateData( par, T, x0 ):
    # Pre-allocate vectors for state (x)
    # and measurements (y)
    x = np.zeros( T+1 )
    y = np.zeros( T   )

    # Set the initial state
    x[0] = x0;

    # Simulate for each time step
    for tt in range(1,T):
        x[tt] = par[0] * x[tt-1] + par[1] * np.random.randn();
        y[tt] =          x[tt]   + par[2] * np.random.randn();

    return( x, y)


##############################################################################
#
# Fully-adapted particle filter for the linear Gaussian SSM
#
# Inputs:
# y:                   observations from the system for t=1,...,T.
#
# par:                 the persistence of the state par[0] and the
#                      standard deviations of the state noise par[1] and
#                      observation noise par[2].
#
# nPart:               number of particles (N)
#
# T and xo:            the no. observations and initial state.
#
# Outputs:
# xh:                  vector with T elements
#                      estimates of the filtered state
#                      for each t=0,1,...,T-1.
#
# ll:                  estimate of the log-likelihood at T-1
#
##############################################################################

def pf(y,par,nPart,T,x0):

    # Initalise variables
    p   = np.zeros((nPart,T));
    v   = np.zeros((nPart,T));
    w   = np.zeros((nPart,T));
    xh  = np.zeros((T,1));
    ll  = np.zeros(T);

    # Set the initial state and weight
    p[:,0] = x0;
    w[:,0] = 1.0 / nPart;
    v[:,0] = 1.0;

    #=====================================================================
    # Run main loop
    #=====================================================================
    for tt in range(1, T):

        #=============================================================
        # Resample particles
        #=============================================================
        nIdx     = np.random.choice(nPart, nPart, p=w[:,tt-1], replace=True );

        #=============================================================
        # Propagate particles
        #=============================================================
        p[:,tt] = par[0] * p[nIdx,tt-1] + par[1] * np.random.randn(1,nPart)

        #=================================================================
        # Weight particles
        #=================================================================

        # Compute log-weights
        v[:,tt] = norm.logpdf( y[tt-1], p[:,tt], par[2] );

        # Rescale log-weights and recover weights
        vmax    = np.max( v[:,tt] );
        v[:,tt] = np.exp( v[:,tt] - vmax );

        # Normalise the weights
        w[:,tt]  = v[:,tt] / np.sum( v[:,tt] );

        # Estimate the filtered state
        xh[tt]  = np.sum( w[:,tt] * p[:,tt] );

         # Estimate log-likelihood
        ll[tt]   = vmax + np.log( np.sum( v[:,tt] ) ) - np.log(nPart);

    #=====================================================================
    # Return state estimate and log-likelihood estimate
    #=====================================================================
    return ( xh, np.sum(ll) )


###################################################################################
#
# Kalman filter for the linear Gaussian SSM
#
# Inputs:
# y:                   observations from the system for t=1,...,T.
#
# par:                 the persistence of the state par[0] and the
#                      standard deviations of the state noise par[1] and
#                      observation noise par[2].
#
# x0 and P0:           the initial state and the corresponding covariance
#
# Outputs:
# xh:                  vector with T elements
#                      estimates of the filtered state
#                      for each t=0,1,...,T-1.
#
# ll:                  estimate of the log-likelihood at T-1
#
###################################################################################

def kf(y,par,T,x0,P0):

    # Initalise variables
    Pp  = P0;
    xp  = x0 * np.ones((T+1,1));
    xf  = x0 * np.ones((T,1));

    # Set parameters
    A = par[0];
    C = 1;
    Q = par[1]**2;
    R = par[2]**2;

    # Run main loop
    for tt in range(0, T):

        # Calculate the Kalman Gain
        S = C  * Pp * C + R;
        K = Pp * C  / S;

        # Compute the state estimate
        yp       = C * xp[tt];
        xf[tt]   = xp[tt] + K * ( y[tt-1] - yp );
        xp[tt+1] = A * xf[tt];

        # Update covariance
        Pf    = Pp - K * S * K;
        Pp    = A * Pf * A + Q;

    return xf;

##############################################################################
#
# Fully-adapted particle filter for the stochastic volatility model
#
# Inputs:
# y:                   observations from the system for t=1,...,T.
#
# par:                 the mean of the state par[0]
#                      the persistence of the state par[1] and the
#                      standard deviations of the state noise par[2].
#
# nPart:               number of particles (N)
#
# T:                   the no. observations.
#
# Outputs:
# xh:                  vector with T elements
#                      estimates of the filtered state
#                      for each t=0,1,...,T-1.
#
# ll:                  estimate of the log-likelihood at T-1
#
##############################################################################

def pf_sv(y,par,nPart,T):

    # Initalise variables
    a   = np.zeros((nPart,T));
    p   = np.zeros((nPart,T));
    v   = np.zeros((nPart,T));
    w   = np.zeros((nPart,T));
    ll  = np.zeros(T);

    # Set the initial state and weight
    p[:,0] = par[0] + par[2] / np.sqrt( 1.0 - par[1]**2 ) * np.random.normal( size=nPart );
    w[:,0] = 1.0 / nPart;
    v[:,0] = 1.0;

    #=====================================================================
    # Run main loop
    #=====================================================================
    for tt in range(1, T):

        #=============================================================
        # Resample particles
        #=============================================================
        idx         = np.random.choice(nPart, nPart, p=w[:,tt-1], replace=True );

        # Resample the ancestory linage
        a[:,1:tt-1]  = a[idx,1:tt-1];

        # Add the most recent ancestors
        a[:,tt]      = idx;

        #=============================================================
        # Propagate particles
        #=============================================================
        p[:,tt] = par[0] + par[1] * ( p[idx,tt-1] - par[0] ) + par[2] * np.random.randn(1,nPart)

        #=================================================================
        # Weight particles
        #=================================================================

        # Compute log-weights
        v[:,tt] = norm.logpdf( y[tt-1], 0, np.exp(p[:,tt] / 2) );

        # Rescale log-weights and recover weights
        vmax    = np.max( v[:,tt] );
        v[:,tt] = np.exp( v[:,tt] - vmax );

        # Normalise the weights
        w[:,tt]  = v[:,tt] / np.sum( v[:,tt] );


        # Estimate log-likelihood
        ll[tt]   = vmax + np.log( np.sum( v[:,tt] ) ) - np.log(nPart);

    #=====================================================================
    # Return state estimate and log-likelihood estimate
    #=====================================================================

    # Sample the state estimate using the weights at tt=T
    nIdx  = np.random.choice( nPart, 1, p=w[:,T-1] );
    xhatf = p[ a[nIdx,T-1].astype(int), : ]

    return ( xhatf, np.sum(ll) )
