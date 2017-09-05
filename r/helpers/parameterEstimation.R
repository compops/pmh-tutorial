# Particle Metropolis-Hastings (LGSS model)
particleMetropolisHastings <-
  function(y,
           initialPhi,
           sigmav,
           sigmae,
           noParticles,
           initialState,
           noIterations,
           stepSize) {
    
    phi <- matrix(0, nrow = noIterations, ncol = 1)
    phiProposed <- matrix(0, nrow = noIterations, ncol = 1)
    logLikelihood <- matrix(0, nrow = noIterations, ncol = 1)
    logLikelihoodProposed <- matrix(0, nrow = noIterations, ncol = 1)
    proposedPhiAccepted <- matrix(0, nrow = noIterations, ncol = 1)
    
    # Set the initial parameter and estimate the initial log-likelihood
    phi[1] <- initialPhi
    theta <- c(phi[1], sigmav, sigmae)
    outputPF <- particleFilter(y, theta, noParticles, initialState)
    logLikelihood[1]<- outputPF$logLikelihood
    
    # Main loop
    for (k in 2:noIterations) {
      # Propose a new parameter
      phiProposed[k] <- phi[k - 1] + stepSize * rnorm(1)
      
      # Estimate the log-likelihood (don't run if unstable system)
      if (abs(phiProposed[k]) < 1.0) {
        theta <- c(phiProposed[k], sigmav, sigmae)
        outputPF <- particleFilter(y, theta, noParticles, initialState)
        logLikelihoodProposed[k] <- outputPF$logLikelihood
      }
      
      # Compute the acceptance probability
      priorPart <- dnorm(phiProposed[k], log = TRUE)
      priorPart <- priorPart - dnorm(phi[k - 1], log = TRUE)
      likelihoodDifference <- logLikelihoodProposed[k] - logLikelihood[k - 1]
      acceptProbability <- exp(priorPart + likelihoodDifference)
    
      acceptProbability <- acceptProbability * (abs(phiProposed[k]) < 1.0)

      # Accept / reject step
      uniformRandomVariable <- runif(1)

      if (uniformRandomVariable < acceptProbability) {
        # Accept the parameter
        phi[k] <- phiProposed[k]
        logLikelihood[k] <- logLikelihoodProposed[k]
        proposedPhiAccepted[k] <- 1
      } else {
        # Reject the parameter
        phi[k] <- phi[k - 1]
        logLikelihood[k] <- logLikelihood[k - 1]
        proposedPhiAccepted[k] <- 0
      }
      
      # Write out progress
      if (k %% 100 == 0) {
        cat(
          sprintf(
            "#####################################################################\n"
          )
        )
        cat(sprintf(" Iteration: %d of : %d completed.\n \n", k, noIterations))
        cat(sprintf(" Current state of the Markov chain:       %.4f \n", phi[k]))
        cat(sprintf(" Proposed next state of the Markov chain: %.4f \n", phiProposed[k]))
        cat(sprintf(" Current posterior mean:                  %.4f \n", mean(phi[0:k])))
        cat(sprintf(" Current acceptance rate:                 %.4f \n", mean(proposedPhiAccepted[0:k])))
        cat(
          sprintf(
            "#####################################################################\n"
          )
        )
      }
    }
    
    phi
  }


# Particle Metropolis-Hastings (SV model)
particleMetropolisHastingsSVmodel <- function(y, initialTheta, noParticles, noIterations, stepSize) {
  
  T <- length(y) - 1

  xHatFiltered <- matrix(0, nrow = noIterations, ncol = T + 1)
  xHatFilteredProposed <- matrix(0, nrow = noIterations, ncol = T + 1)
  theta <- matrix(0, nrow = noIterations, ncol = 3)
  thetaProposed <- matrix(0, nrow = noIterations, ncol = 3)
  logLikelihood <- matrix(0, nrow = noIterations, ncol = 1)
  logLikelihoodProposed <- matrix(0, nrow = noIterations, ncol = 1)
  proposedThetaAccepted <- matrix(0, nrow = noIterations, ncol = 1)
  
  # Set the initial parameter and estimate the initial log-likelihood
  theta[1, ] <- initialTheta
  res <- particleFilterSVmodel(y, theta[1, ], noParticles)
  logLikelihood[1] <- res$logLikelihood
  xHatFiltered[1, ] <- res$xHatFiltered
  
  # Main loop
  for (k in 2:noIterations) {
    # Propose a new parameter
    thetaProposed[k, ] <- rmvnorm(1, mean = theta[k - 1, ], sigma = stepSize)
    
    # Estimate the log-likelihood (don't run if unstable system)
    if ((abs(thetaProposed[k, 2]) < 1.0) && (thetaProposed[k, 3] > 0.0)) {
      res <- particleFilterSVmodel(y, thetaProposed[k, ], noParticles)
      logLikelihoodProposed[k]  <- res$logLikelihood
      xHatFilteredProposed[k, ] <- res$xHatFiltered
    }
    
    # Compute difference in the log-priors
    priorMu <- dnorm(thetaProposed[k, 1], 0, 1, log = TRUE) 
    priorMu <- priorMu - dnorm(theta[k - 1, 1], 0, 1, log = TRUE)
    priorPhi <- dnorm(thetaProposed[k, 2], 0.95, 0.05, log = TRUE) 
    priorPhi <- priorPhi - dnorm(theta[k - 1, 2], 0.95, 0.05, log = TRUE)
    priorSigmaV <- dgamma(thetaProposed[k, 3], 2, 10, log = TRUE)
    priorSigmaV <- priorSigmaV - dgamma(theta[k - 1, 3], 2, 10, log = TRUE)
    prior <- priorMu + priorPhi + priorSigmaV
    
    # Compute the acceptance probability
    likelihoodDifference <- logLikelihoodProposed[k] - logLikelihood[k - 1]
    acceptProbability <- exp(prior + likelihoodDifference)
    
    acceptProbability <- acceptProbability * (abs(thetaProposed[k, 2]) < 1.0)
    acceptProbability <- acceptProbability * (thetaProposed[k, 3] > 0.0)
    
    # Accept / reject step
    uniformRandomVariable <- runif(1)

    if (uniformRandomVariable < acceptProbability) {
      # Accept the parameter
      theta[k, ] <- thetaProposed[k, ]
      logLikelihood[k] <- logLikelihoodProposed[k]
      xHatFiltered[k, ] <- xHatFilteredProposed[k, ]
      proposedThetaAccepted[k] <- 1
    } else {
      # Reject the parameter
      theta[k, ]  <- theta[k - 1, ]
      logLikelihood[k] <- logLikelihood[k - 1]
      xHatFiltered[k, ] <- xHatFiltered[k - 1, ]
      proposedThetaAccepted[k] <- 0
    }
    
    # Write out progress
    if (k %% 100 == 0) {
      cat(
        sprintf(
          "#####################################################################\n"
        )
      )
      cat(sprintf(" Iteration: %d of : %d completed.\n \n", k, noIterations))
      
      cat(sprintf(
        " Current state of the Markov chain:       %.4f %.4f %.4f \n",
        theta[k, 1],
        theta[k, 2],
        theta[k, 3]
      ))
      cat(
        sprintf(
          " Proposed next state of the Markov chain: %.4f %.4f %.4f \n",
          thetaProposed[k, 1],
          thetaProposed[k, 2],
          thetaProposed[k, 3]
        )
      )
      cat(sprintf(
        " Current posterior mean:                  %.4f %.4f %.4f \n",
        mean(thetaProposed[0:k, 1]),
        mean(thetaProposed[0:k, 2]),
        mean(thetaProposed[0:k, 3])
      ))
      cat(sprintf(" Current acceptance rate:                 %.4f \n", mean(proposedThetaAccepted[0:k])))
      cat(
        sprintf(
          "#####################################################################\n"
        )
      )
    }
  }
  
  list(theta = theta, xHatFiltered = xHatFiltered, proposedThetaAccepted = proposedThetaAccepted)
}


# Particle Metropolis-Hastings (reparameterised SV model)
particleMetropolisHastingsSVmodelReparameterised <-
    function(y, initialTheta, noParticles, noIterations, stepSize) {
          
    T <- length(y) - 1
    
    xHatFiltered <- matrix(0, nrow = noIterations, ncol = T + 1)
    xHatFilteredProposed <- matrix(0, nrow = noIterations, ncol = T + 1)
    theta <- matrix(0, nrow = noIterations, ncol = 3)
    thetaProposed <- matrix(0, nrow = noIterations, ncol = 3)
    thetaTransformed <- matrix(0, nrow = noIterations, ncol = 3)
    thetaTransformedProposed <- matrix(0, nrow = noIterations, ncol = 3)
    logLikelihood <- matrix(0, nrow = noIterations, ncol = 1)
    logLikelihoodProposed <- matrix(0, nrow = noIterations, ncol = 1)
    proposedThetaAccepted <- matrix(0, nrow = noIterations, ncol = 1)      
    
    # Set the initial parameter and estimate the initial log-likelihood
    theta[1, ] <- initialTheta
    res <- particleFilterSVmodel(y, theta[1, ], noParticles)
    thetaTransformed[1, ] <- c(theta[1, 1], atanh(theta[1, 2]), log(theta[1, 3]))
    logLikelihood[1] <- res$logLikelihood
    xHatFiltered[1, ] <- res$xHatFiltered
    
    # Main loop
    for (k in 2:noIterations) {
      # Propose a new parameter
      thetaTransformedProposed[k, ] <- rmvnorm(1, mean = thetaTransformed[k - 1, ], sigma = stepSize)
      
      # Run the particle filter
      thetaProposed[k, ] <- c(thetaTransformedProposed[k, 1], tanh(thetaTransformedProposed[k, 2]), exp(thetaTransformedProposed[k, 3]))
      res <- particleFilterSVmodel(y, thetaProposed[k, ], noParticles)
      xHatFilteredProposed[k, ] <- res$xHatFiltered
      logLikelihoodProposed[k] <- res$logLikelihood
      
      # Compute the acceptance probability
      logPrior1 <- dnorm(thetaProposed[k, 1], log = TRUE) - dnorm(theta[k - 1, 1], log = TRUE)
      logPrior2 <-dnorm(thetaProposed[k, 2], 0.95, 0.05, log = TRUE) - dnorm(theta[k - 1, 2], 0.95, 0.05, log = TRUE)
      logPrior3 <- dgamma(thetaProposed[k, 3], 3, 10, log = TRUE) - dgamma(theta[k - 1, 3], 3, 10, log = TRUE)
      logPrior <- logPrior1 + logPrior2 + logPrior3
      
      logJacob1 <- log(abs(1 - thetaProposed[k, 2]^2)) -log(abs(1 - theta[k - 1, 2]^2))
      logJacob2 <- log(abs(thetaProposed[k, 3])) - log(abs(theta[k - 1, 3]))
      logJacob <- logJacob1 + logJacob2
      
      acceptProbability <- exp(logPrior + logLikelihoodProposed[k] - logLikelihood[k - 1] + logJacob)

      # Accept / reject step
      uniformRandomVariable <- runif(1)

      if (uniformRandomVariable < acceptProbability) {
        # Accept the parameter
        theta[k, ] <- thetaProposed[k, ]
        thetaTransformed[k, ] <- thetaTransformedProposed[k, ]
        logLikelihood[k] <- logLikelihoodProposed[k]
        xHatFiltered[k, ] <- xHatFilteredProposed[k, ]
        proposedThetaAccepted[k] <- 1
      } else {
        # Reject the parameter
        theta[k, ] <- theta[k - 1, ]
        thetaTransformed[k, ] <- thetaTransformed[k - 1, ]
        logLikelihood[k] <- logLikelihood[k - 1]
        xHatFiltered[k, ] <- xHatFiltered[k - 1, ]
        proposedThetaAccepted[k]  <- 0
      }
      
      # Write out progress
      if (k %% 100 == 0) {
        cat(
          sprintf(
            "#####################################################################\n"
          )
        )
        cat(sprintf(" Iteration: %d of : %d completed.\n \n", k, noIterations))
        cat(sprintf(
          " Current state of the Markov chain:       %.4f %.4f %.4f \n",
          thetaTransformed[k, 1],
          thetaTransformed[k, 2],
          thetaTransformed[k, 3]
        ))
        cat(
          sprintf(
            " Proposed next state of the Markov chain: %.4f %.4f %.4f \n",
            thetaTransformedProposed[k, 1],
            thetaTransformedProposed[k, 2],
            thetaTransformedProposed[k, 3]
          )
        )
        cat(sprintf(
          " Current posterior mean:                  %.4f %.4f %.4f \n",
          mean(theta[0:k, 1]),
          mean(theta[0:k, 2]),
          mean(theta[0:k, 3])
        ))
        cat(sprintf(" Current acceptance rate:                 %.4f \n", mean(proposedThetaAccepted[0:k])))
        cat(
          sprintf(
            "#####################################################################\n"
          )
        )
        
      }
    }
    
    list(theta = theta,
         xHatFiltered = xHatFiltered,
         thetaTransformed = thetaTransformed)
    }