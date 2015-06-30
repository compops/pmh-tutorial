##############################################################################
#
# Example of particle Metropolis-Hastings 
#
# Subroutine for particle Metropolis-Hastings
#
# (c) 2015 Johan Dahlin
# johan.dahlin (at) liu.se
#
##############################################################################

##############################################################################
# Particle Metropolis-Hastings (LGSS model)
##############################################################################

pmh <- function(y,initPar,sigmav,sigmae,nPart,T,xo,nIter,stepSize) {
  
  # Initalise variables
  th     = matrix(0,nrow=nIter,ncol=1);
  thp    = matrix(0,nrow=nIter,ncol=1);
  ll     = matrix(0,nrow=nIter,ncol=1);
  llp    = matrix(0,nrow=nIter,ncol=1);
  accept = matrix(0,nrow=nIter,ncol=1);
  
  # Set the initial parameter and estimate the initial log-likelihood
  th[1]  = initPar;  
  ll[1]  = sm(y,th[1],sigmav,sigmae,nPart,T,xo)$ll
  
  #=====================================================================
  # Run main loop
  #=====================================================================
  for ( kk in 2:nIter ) {
    
    # Propose a new parameter
    thp[kk] = th[kk-1] + stepSize * rnorm(1);
  
    # Estimate the log-likelihood (don't run if unstable system)
    if ( abs( thp[kk] ) < 1.0 ) {
      llp[kk] <- sm(y,thp[kk],sigmav,sigmae,nPart,T,xo)$ll
    }
    
    # Compute the acceptance probability
    aprob = exp( llp[kk] - ll[kk-1] );
    
    # Generate uniform random variable in U[0,1]
    u = runif(1);
    
    # Accept / reject step
    # Check if | phi | > 1.0, in that case always reject.
    if ( (u < aprob) && ( abs(thp[kk]) < 1.0 ) ) {
      
      # Accept the parameter
      th[kk]     = thp[kk]
      ll[kk]     = llp[kk]
      accept[kk] = 1.0;
      
    } else {
      # Reject the parameter
      th[kk]     = th[kk-1]
      ll[kk]     = ll[kk-1]        
      accept[kk] = 0.0;
    }
    
    # Write out progress
    if ( kk%%100 == 0 ) {
      cat(sprintf("#####################################################################\n"))
      cat(sprintf(" Iteration: %d of : %d completed.\n \n", kk, nIter))
      cat(sprintf(" Current state of the Markov chain:       %.4f \n", th[kk] ))
      cat(sprintf(" Proposed next state of the Markov chain: %.4f \n", thp[kk] ))
      cat(sprintf(" Current posterior mean:                  %.4f \n", mean(th[0:kk]) ))
      cat(sprintf(" Current acceptance rate:                 %.4f \n", mean(accept[0:kk]) ))
      cat(sprintf("#####################################################################\n"))
    }
  }
  
  #=====================================================================
  # Return traces of the parameters
  #=====================================================================
  th;
}


##############################################################################
# Particle Metropolis-Hastings (SV model)
##############################################################################

pmh_sv <- function(y,initPar,nPart,T,nIter,stepSize) {
  
  # Initalise variables
  th     = matrix(0,nrow=nIter,ncol=3);
  thp    = matrix(0,nrow=nIter,ncol=3);
  ll     = matrix(0,nrow=nIter,ncol=1);
  llp    = matrix(0,nrow=nIter,ncol=1);
  accept = matrix(0,nrow=nIter,ncol=1);
  
  # Set the initial parameter and estimate the initial log-likelihood
  th[1,]  = initPar;  
  ll[1]   = sm_sv(y,th[1,1],th[1,2],th[1,3],nPart,T)$ll
  
  #=====================================================================
  # Run main loop
  #=====================================================================
  for ( kk in 2:nIter ) {
    
    # Propose a new parameter
    thp[kk,] = th[kk-1,] + stepSize * rnorm(3);
    
    # Estimate the log-likelihood (don't run if unstable system)
    if ( ( abs(thp[kk,1]) < 1.0 ) && ( thp[kk,2] > 0.0 ) && ( thp[kk,3] > 0.0 ) ) {
      llp[kk] <- sm_sv(y,thp[kk,1],thp[kk,2],thp[kk,3],nPart,T)$ll
    }
    
    # Compute the acceptance probability
    aprob = exp( llp[kk] - ll[kk-1] );
    
    # Generate uniform random variable in U[0,1]
    u = runif(1);
    
    # Accept / reject step
    # Check if | phi | > 1.0, in that case always reject.
    # Check if sigmav < 0.0, in that case always reject.
    # Check if sigmae < 0.0, in that case always reject.
    if ( (u < aprob) && ( abs(thp[kk,1]) < 1.0 ) && ( thp[kk,2] > 0.0 ) && ( thp[kk,3] > 0.0 ) ) {
      
      # Accept the parameter
      th[kk,]     = thp[kk,]
      ll[kk]     = llp[kk]
      accept[kk] = 1.0;
      
    } else {
      # Reject the parameter
      th[kk,]     = th[kk-1,]
      ll[kk]     = ll[kk-1]        
      accept[kk] = 0.0;
    }
    
    # Write out progress
    if ( kk%%100 == 0 ) {
      cat(sprintf("#####################################################################\n"))
      cat(sprintf(" Iteration: %d of : %d completed.\n \n", kk, nIter))
      cat(sprintf(" Current state of the Markov chain:       %.4f %.4f %.4f \n", th[kk,1], th[kk,2], th[kk,3] ))
      cat(sprintf(" Proposed next state of the Markov chain: %.4f %.4f %.4f \n", thp[kk,1], thp[kk,2], thp[kk,3] ))
      cat(sprintf(" Current posterior mean:                  %.4f %.4f %.4f \n", mean(th[0:kk,1]), mean(th[0:kk,2]), mean(th[0:kk,3]) ))
      cat(sprintf(" Current acceptance rate:                 %.4f \n", mean(accept[0:kk]) ))
      cat(sprintf("#####################################################################\n"))
    }
  }
  
  #=====================================================================
  # Return traces of the parameters
  #=====================================================================
  th;
}