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

##############################################################################
# Generate data for LGSS model
##############################################################################

generateData <- function(phi,sigmaV,sigmaE,T,xo)
{
  # Pre-allocate vectors for log-volatility/state (x) 
  # and log-returns/observations (y)
  x    = matrix(0, nrow=T, ncol=1)
  y    = matrix(0, nrow=T, ncol=1)
  
  # Set the initial state
  x[1] = xo;
  y[1] = xo + sigmaE * rnorm(1)
  
  # Simulate the system for each time step
  for ( tt in 2:T ) {
    x[tt] = phi  * x[tt-1] + sigmaV * rnorm(1);
    y[tt] =        x[tt]   + sigmaE * rnorm(1);
  }
  
  data.frame(x=x, y=y)
}


##############################################################################
# Bootstrap particle filter (LGSS)
##############################################################################

sm <- function(y,phi,sigmav,sigmae,nPart,T,x0)
{
  # Initalise variables
  xhatf = matrix( x0,      nrow=T+1, ncol=1)
  p     = matrix( x0,      nrow=nPart, ncol=T+1)
  w     = matrix( 1/nPart, nrow=nPart, ncol=T+1)
  ll    = 0;

  #===========================================================
  # Run main loop
  #===========================================================
  for ( tt in 2:(T+1) )
  {
    #=========================================================
    # Resample ( multinomial )
    #=========================================================
    nIdx   = sample(1:nPart, nPart, replace=T, prob = w[,tt-1] )
    
    #=========================================================
    # Propagate
    #=========================================================
    p[,tt] = phi * p[nIdx,tt-1] + rnorm(nPart, 0, sigmav ) 
    
    #=========================================================
    # Compute weights
    #=========================================================
    w[,tt] = dnorm( y[tt-1], p[,tt], sigmae, log=T)
    
    # Rescale log-weights and recover weight
    wmax   = max(w[,tt]);
    w[,tt] = exp(w[,tt] - wmax)
    
    # Estimate the log-likelihood
    ll     = ll + wmax + log(sum(w[,tt])) - log(nPart);
    
    # Normalize the weights
    w[,tt] = w[,tt] / sum(w[,tt])
    
    # Estimate the state
    xhatf[tt] = sum( w[,tt] * p[,tt] )
    
  }
  #===========================================================
  # Return state estimate and log-likelihood estimate
  #===========================================================
  output = list( xh = xhatf, ll=ll)
}

###################################################################################
# Kalman filter (LGSS)
###################################################################################
kf <- function(y,phi,sigmaV,sigmaE,x0,P0)
{
  xhatf = matrix( x0, nrow=T, ncol=1)
  xhatp = matrix( x0, nrow=T, ncol=1)
  Pp       = P0;
  T        = length(y);
  
  # Set parameters 
  A = phi;
  C = 1;
  Q = sigmaV^2;
  R = sigmaE^2;
  
  for ( tt in 1:T )
  {
    # Compute Kalman Gain
    S = C * Pp * C + R;
    K = Pp * C / S;
    
    # Compute state estimate
    xhatf[tt]   = xhatp[tt] + K * ( y[tt] - C * xhatp[tt] );
    xhatp[tt+1] = A * xhatf[tt]; 
    
    # Update covariance
    Pf = Pp - K * S * K;
    Pp = A * Pf * A + Q;
  }
  output = list( xh = xhatf )
}

##############################################################################
# Bootstrap particle filter (SV model)
##############################################################################

sm_sv <- function(y,phi,sigmav,beta,N,T)
{
  # Initalise variables
  xhatf = matrix( 0,       nrow=T+1,   ncol=1)
  p     = matrix( 0,       nrow=nPart, ncol=T+1)
  w     = matrix( 1/nPart, nrow=nPart, ncol=T+1)
  ll    = 0;
  
  # Generate initial state
  p[,1]     = rnorm(nPart, 0, sigmav / sqrt( 1 - phi*phi ) );
  xhatf[,1] = mean(p[,1]);
  
  #===========================================================
  # Run main loop
  #===========================================================
  for ( tt in 2:(T+1) )
  {
    #=========================================================
    # Resample ( multinomial )
    #=========================================================
    idx = sample(1:nPart, nPart, replace=T, prob = w[,tt-1] )
    
    #=========================================================
    # Propagate
    #=========================================================
    p[,tt] = phi * p[idx,tt-1] + rnorm(nPart, 0, sigmav ) 
    
    #=========================================================
    # Compute weights
    #=========================================================
    w[,tt] = dnorm( y[tt-1], 0, beta*exp(p[,tt]/2), log=T)
    
    # Rescale log-weights and recover weight
    wmax   = max(w[,tt]);
    w[,tt] = exp(w[,tt] - wmax)
    
    # Estimate the log-likelihood
    ll     = ll + wmax + log(sum(w[,tt])) - log(nPart);
    
    # Normalize the weights
    w[,tt] = w[,tt] / sum(w[,tt])
    
    # Estimate the state
    xhatf[tt] = sum( w[,tt] * p[,tt] )
    
  }
  #===========================================================
  # Return state estimate and log-likelihood estimate
  #===========================================================
  output = list( xh = xhatf, ll=ll)
}
