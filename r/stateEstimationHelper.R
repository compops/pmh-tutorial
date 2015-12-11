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

##############################################################################
# Generate data for LGSS model
##############################################################################

generateData <- function(phi,sigmav,sigmae,T,xo)
{
  #
  # Generates data from the LGSS model with parameters (phi,sigmav,sigmae)
  #
  # Inputs:
  # phi, sigmav, sigmae: the persistence of the state and the 
  #                      standard deviations of the state innovations and 
  #                      observation noise.
  #
  # T and xo:            the no. observations and initial state.
  #
  # Outputs:
  # x,y:                 the latent state and observations
  #
  #
  
  # Pre-allocate vectors for log-volatility/state (x) 
  # and log-returns/observations (y)
  x    = matrix( 0, nrow=T+1, ncol=1 );
  y    = matrix( 0, nrow=T+1, ncol=1 );
  
  # Set the initial state
  x[1] = xo;
  y[1] = NA;
  
  # Simulate the system for each time step
  for ( tt in 2:(T+1) ) {
    x[tt] = phi  * x[tt-1] + sigmav * rnorm(1);
    y[tt] =        x[tt]   + sigmae * rnorm(1);
  }
  
  list(x=x, y=y);
}


##############################################################################
# Fully-adapted particle filter (LGSS)
##############################################################################

sm <- function(y,phi,sigmav,sigmae,nPart,T,x0)
{
  #
  # Fully-adapted particle filter for the linear Gaussian SSM
  #
  # Inputs:
  # y:                   observations from the system for t=1,...,T.
  #  
  # phi, sigmav, sigmae: the persistence of the state and the 
  #                      standard deviations of the state innovations and 
  #                      observation noise.
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
  #
  
  #===========================================================
  # Initialise variables
  #===========================================================
  xhatf = matrix( x0,      nrow=T,     ncol=1   );
  p     = matrix( x0,      nrow=nPart, ncol=T+1 );
  w     = matrix( 1/nPart, nrow=nPart, ncol=T+1 );
  v     = matrix( 1      , nrow=nPart, ncol=T+1 );
  ll    = 0;

  #===========================================================
  # Run main loop
  #===========================================================
  for ( tt in 2:T )
  {
    #=========================================================
    # Resample ( multinomial )
    #=========================================================
    nIdx   = sample( 1:nPart, nPart, replace=T, prob = w[,tt-1] );
    
    #=========================================================
    # Propagate
    #=========================================================
    Delta  = ( sigmav^(-2) + sigmae^(-2) )^(-1);
    mup    = sigmae^(-2) * y[tt] + sigmav^(-2) * phi * p[nIdx,tt-1];
    p[,tt] = Delta * mup + rnorm( nPart, 0, sqrt( Delta ) );
    
    #=========================================================
    # Compute weights
    #=========================================================
    v[,tt] = dnorm( y[tt+1], phi * p[,tt], sqrt( sigmae^2 + sigmav^2), log=T );
    
    # Rescale log-weights and recover weight
    vmax   = max( v[,tt] );
    v[,tt] = exp( v[,tt] - vmax );
    
    # Normalize the weights
    w[,tt] = v[,tt] / sum( v[,tt] );
        
    # Estimate the state
    xhatf[tt] = mean( p[,tt] );
    
    # Estimate the log-likelihood
    ll     = ll + vmax + log( sum(v[,tt]) ) - log(nPart);    
    
  }
  
  #===========================================================
  # Return state estimate and log-likelihood estimate
  #===========================================================
  list( xh = xhatf, ll=ll, p=p, w=w);
  
}

###################################################################################
# Kalman filter (LGSS)
###################################################################################
kf <- function(y,phi,sigmav,sigmae,x0,P0)
{
  #
  # Kalman filter for the linear Gaussian SSM
  #
  # Inputs:
  # y:                   observations from the system for t=1,...,T.
  #  
  # phi, sigmav, sigmae: the persistence of the state and the 
  #                      standard deviations of the state innovations and 
  #                      observation noise.
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
  #
    
  T     = length(y);
  yhatp = matrix( x0, nrow=T,   ncol=1 );
  xhatf = matrix( x0, nrow=T,   ncol=1 );
  xhatp = matrix( x0, nrow=T+1, ncol=1 );
  Pp       = P0;
  ll       = 0;
  
  # Set parameters 
  A = phi;
  C = 1;
  Q = sigmav^2;
  R = sigmae^2;
  
  for ( tt in 2:T)
  {
    # Compute Kalman Gain
    S = C * Pp * C + R;
    K = Pp * C / S;
    
    # Compute state estimate
    yhatp[tt]   = C * xhatp[tt];
    xhatf[tt]   = xhatp[tt] + K * ( y[tt] - yhatp[tt] );
    xhatp[tt+1] = A * xhatf[tt];
    
    # Update covariance
    Pf = Pp - K * S * K;
    Pp = A * Pf * A + Q;
    
    # Estimate loglikelihood (not in the last iteration, to be able to compare with faPF)
    if ( tt < T ) { ll = ll + dnorm( y[tt], yhatp[tt], sqrt(S), log=T) };
  }
  
  list( xh = xhatf, ll=ll )
}

##############################################################################
# Bootstrap particle filter (SV model)
##############################################################################
sm_sv <- function(y,mu,phi,sigmav,N,T)
{
  
  #
  # Bootstrap particle filter for the stochastic volatility model
  #
  # Inputs:
  # y:                   observations from the system for t=1,...,T.
  #
  # mu, phi, sigmav:     the mean and persistence of the state and the 
  #                      standard deviations of the state innovations.
  #
  # nPart:               number of particles (N)
  #
  # T and xo:            the no. observations and initial state.
  #
  # Outputs:
  # xh:                  vector with T+1 elements
  #                      estimates of the smoothed state
  #                      for each t=0,1,...,T.
  #
  # ll:                  estimate of the log-likelihood at T
  #
  #

  #===========================================================
  # Initialise variables
  #===========================================================
  a     = matrix( 0,       nrow=nPart, ncol=T+1 );
  p     = matrix( 0,       nrow=nPart, ncol=T+1);
  w     = matrix( 1/nPart, nrow=nPart, ncol=T+1);
  v     = matrix( 1      , nrow=nPart, ncol=T+1 );
  ll    = 0;
  
  # Generate initial state
  p[,1]     = rnorm(nPart, mu, sigmav / sqrt( 1 - phi * phi ) );
  a[,1]     = 1:N;
  
  #===========================================================
  # Run main loop
  #===========================================================
  for ( tt in 2:(T+1) )
  {
    #=========================================================
    # Resample ( multinomial )
    #=========================================================
    idx = sample(1:nPart, nPart, replace=T, prob = w[,tt-1] );
    
    # Resample the ancestory linage
    a[,1:tt-1]  = a[idx,1:tt-1];
    
    # Add the most recent ancestors
    a[,tt]      = idx;
    
    #=========================================================
    # Propagate
    #=========================================================
    p[,tt] = mu + phi * ( p[idx,tt-1] - mu ) + rnorm(nPart, 0, sigmav ) ;
    
    #=========================================================
    # Compute weights
    #=========================================================
    v[,tt] = dnorm( y[tt-1], 0, exp( p[,tt] / 2 ), log=T);
    
    # Rescale log-weights and recover weight
    vmax   = max(v[,tt]);
    v[,tt] = exp(v[,tt] - vmax);
    
    # Normalize the weights
    w[,tt] = v[,tt] / sum( v[,tt] );
    
    # Estimate the log-likelihood
    ll     = ll + vmax + log( sum(v[,tt]) ) - log(nPart);    
    
  }
  #===========================================================
  # Return state estimate and log-likelihood estimate
  #===========================================================

  # Sample the state estimate using the weights at tt=T
  nIdx  = sample( 1:nPart, 1, prob=w[,T] )
  xhatf = p[ cbind( a[ nIdx, ] , 1:(T+1) ) ]  
  
  list( xh = xhatf, ll=ll)
}
