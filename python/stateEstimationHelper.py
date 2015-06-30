##############################################################################
#
# Example of particle filtering
#
# Subroutine for data generation and particle filtering
#
# (c) 2015 Johan Dahlin
# johan.dahlin (at) liu.se
#
##############################################################################

import numpy as np
from   scipy.stats import norm


##############################################################################
# Generate data for LGSS model
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
# Bootstrap particle filter (LGSS)
##############################################################################

def pf(y,par,nPart,T,xo):

    # Initalise variables
    p   = np.zeros((nPart,T));
    w   = np.zeros((nPart,T));
    xh  = np.zeros((T,1));
    ll  = np.zeros(T);

    # Set the initial state and weight
    p[:,0] = xo;
    w[:,0] = 1.0 / nPart;

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
        w[:,tt] = norm.logpdf( y[tt-1], p[:,tt], par[2] );

        # Rescale log-weights and recover weights
        wmax    = np.max( w[:,tt] );
        w[:,tt] = np.exp( w[:,tt] - wmax );

        # Estimate log-likelihood
        ll[tt]   = wmax + np.log( np.sum( w[:,tt] ) ) - np.log(nPart);

        # Normalise the weights
        w[:,tt] /= np.sum( w[:,tt] );

        # Estimate the filtered state
        xh[tt]  = np.sum( w[:,tt] * p[:,tt] );

    #=====================================================================
    # Return state estimate and log-likelihood estimate
    #=====================================================================
    return ( xh, np.sum(ll) )


###################################################################################
# Kalman filter (LGSS)
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

#############################################################################
# Particle filter routine (SV model)
#############################################################################

def pf_sv(y,par,nPart,T):

    # Initalise variables
    p   = np.zeros((nPart,T));
    w   = np.zeros((nPart,T));
    xh  = np.zeros((T,1));
    ll  = np.zeros(T);

    # Set the initial state and weight
    p[:,0] = par[1] / np.sqrt( 1.0 - par[0]**2 ) * np.random.normal( size=nPart );
    xh[0]  = np.mean( p[:,0] )
    w[:,0] = 1.0 / nPart;

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
        w[:,tt] = norm.logpdf( y[tt-1], 0, par[2] * np.exp(p[:,tt]/2) );

        # Rescale log-weights and recover weights
        wmax    = np.max( w[:,tt] );
        w[:,tt] = np.exp( w[:,tt] - wmax );

        # Estimate log-likelihood
        ll[tt]   = wmax + np.log( np.sum( w[:,tt] ) ) - np.log(nPart);

        # Normalise the weights
        w[:,tt] /= np.sum( w[:,tt] );

        # Estimate the filtered state
        xh[tt]  = np.sum( w[:,tt] * p[:,tt] );

    #=====================================================================
    # Return state estimate and log-likelihood estimate
    #=====================================================================
    return ( xh, np.sum(ll) )
