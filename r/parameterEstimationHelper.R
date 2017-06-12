##############################################################################
#
# Example of particle Metropolis-Hastings
#
# Subroutine for particle Metropolis-Hastings
#
# Copyright (C) 2017 Johan Dahlin < liu (at) johandahlin.se >
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
# Particle Metropolis-Hastings (LGSS model)
##############################################################################

particleMetropolisHastings <-
  function(y,
           initialPhi,
           sigmav,
           sigmae,
           noParticles,
           initialState,
           noIterations,
           stepSize) {
    #
    # Particle Metropolis-Hastings (PMH) for the LGSS model
    #
    # Inputs:
    # y:                   observations from the system for t=1,...,T.
    #
    # initPar:             initial value for phi (persistence of the state)
    #
    # sigmav, sigmae:      the standard deviations of the state innovations
    #                      and observation noise.
    #
    # noParticles:         number of particles (N)
    #
    # initialState:        initial state
    #
    # noIterations:        the number of iterations in PMH.
    #
    # stepSize:            the standard deviation of the RW proposal.
    #
    # Outputs:
    # phi:                 K samples from the parameter posterior
    #
    #
    
    #===========================================================
    # Initialise variables
    #===========================================================
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
    
    #=====================================================================
    # Run main loop
    #=====================================================================
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
      
      # Always reject if parameter results in an unstable system
      acceptProbability <- acceptProbability * (abs(phiProposed[k]) < 1.0)
       
      # Generate uniform random variable in U[0,1]
      uniformRandomVariable <- runif(1)
      
      # Accept / reject step
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
    
    #=====================================================================
    # Return traces of the parameter phi
    #=====================================================================
    phi
  }


##############################################################################
# Particle Metropolis-Hastings (SV model)
##############################################################################

particleMetropolisHastingsSVmodel <- function(y, initialTheta, noParticles, noIterations, stepSize) {
  #
  # Particle Metropolis-Hastings (PMH) for the SV model
  #
  # Inputs:
  # y:                   observations from the system for t=1,...,T.
  #
  # initialTheta:        initial values of the parameters
  #                      ( mu, phi, sigmav )
  #
  # noParticles:         number of particles (N)
  #
  # noIterations:        the number of iterations in PMH
  #
  # stepSize:            the standard deviation of the RW proposal.
  #
  # Outputs:
  # theta:               K samples from the parameter posterior.
  #
  #
  
  T <- length(y) - 1
  #===========================================================
  # Initialise variables
  #===========================================================
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
  
  #=====================================================================
  # Run main loop
  #=====================================================================
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
    priorMu <- dnorm(thetaProposed[k, 1], 0, 1, log = TRUE) - dnorm(theta[k - 1, 1], 0, 1, log = TRUE)
    priorPhi <- dnorm(thetaProposed[k, 2], 0.95, 0.05, log = TRUE) - dnorm(theta[k - 1, 2], 0.95, 0.05, log = TRUE)
    priorSigmaV <- dgamma(thetaProposed[k, 3], 2, 10, log = TRUE) - dgamma(theta[k - 1, 3], 2, 10, log = TRUE)
    
    # Compute the acceptance probability
    likelihoodDifference <- logLikelihoodProposed[k] - logLikelihood[k - 1]
    acceptProbability <- exp(priorMu + priorPhi + priorSigmaV + likelihoodDifference)
    
    # Always reject if parameter results in an unstable system
    acceptProbability <- acceptProbability * (abs(thetaProposed[k, 2]) < 1.0)
    
    # Generate uniform random variable in U[0,1]
    uniformRandomVariable <- runif(1)
    
    # Accept / reject step
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
  
  #=====================================================================
  # Return traces of the parameters
  #=====================================================================
  list(theta = theta, xHatFiltered = xHatFiltered)
}


##############################################################################
# Particle Metropolis-Hastings (SV model)
##############################################################################

particleMetropolisHastingsSVmodelReparameterised <-
    function(y, initialTheta, noParticles, noIterations, stepSize) {
    #
    # Particle Metropolis-Hastings (PMH) for the SV model
    #
    # Inputs:
    # y:                   observations from the system for t=1,...,T.
    #
    # initialTheta:        initial values of the parameters
    #                      ( mu, phi, sigmav )
    #
    # noParticles:         number of particles (N)
    #
    # noIterations:        the number of iterations in PMH
    #
    # stepSize:            the standard deviation of the RW proposal.
    #
    # Outputs:
    # theta:               K samples from the parameter posterior.
    #
    #
      
    T <- length(y) - 1
    
    #===========================================================
    # Initialise variables
    #===========================================================
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
    
    #=====================================================================
    # Run main loop
    #=====================================================================
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
      
      # Generate uniform random variable in U[0,1]
      uniformRandomVariable <- runif(1)
      
      # Accept / reject step
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
    
    #=====================================================================
    # Return traces of the parameters
    #=====================================================================
    list(theta = theta,
         xHatFiltered = xHatFiltered,
         thetaTransformed = thetaTransformed)
    }


##############################################################################
# Particle Metropolis-Hastings (SV model)
##############################################################################

makePlotsParticleMetropolisHastingsSVModel <- function(y, res, noBurnInIterations, noIterations, nPlot) {
  
  # Extract the states after burn-in
  resTh <- res$theta[noBurnInIterations:noIterations, ]
  resXh <- res$xHatFiltered[noBurnInIterations:noIterations, ]
  
  # Estimate the posterior mean and the corresponding standard deviation
  thhat   <- colMeans(resTh)
  thhatSD <- apply(resTh, 2, sd)
  
  # Estimate the log-volatility and the corresponding standad deviation
  xhat    <- colMeans(resXh)
  xhatSD  <- apply(resXh, 2, sd)
      
  # Plot the parameter posterior estimate, solid black line indicate posterior mean
  # Plot the trace of the Markov chain after burn-in, solid black line indicate posterior mean
  layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), 5, 3, byrow = TRUE))
  par(mar = c(4, 5, 0, 0))
  
  # Grid for plotting the data and log-volatility
  gridy <- seq(1, length(y))
  gridx <- seq(1, length(y) - 1)
  
  #---------------------------------------------------------------------------
  # Observations
  #---------------------------------------------------------------------------
  plot(
    y,
    col = "#1B9E77",
    lwd = 1,
    type = "l",
    xlab = "time",
    ylab = "log-returns",
    ylim = c(-5, 5),
    bty = "n"
  )
  polygon(
    c(gridy, rev(gridy)),
    c(y, rep(-5, length(gridy))),
    border = NA,
    col = rgb(t(col2rgb("#1B9E77")) / 256, alpha = 0.25)
  )

  #---------------------------------------------------------------------------
  # Log-volatility
  #---------------------------------------------------------------------------  
  plot(
    xhat[-1],
    col = "#D95F02",
    lwd = 1.5,
    type = "l",
    xlab = "time",
    ylab = "log-volatility estimate",
    ylim = c(-2, 2),
    bty = "n"
  )
  xhat_upper <- xhat[-1] + 1.96 * xhatSD[-1]
  xhat_lower <- xhat[-1] - 1.96 * xhatSD[-1]
    
  polygon(
    c(gridx, rev(gridx)),
    c(xhat_upper, rev(xhat_lower)),
    border = NA,
    col = rgb(t(col2rgb("#D95F02")) / 256, alpha = 0.25)
  )
  
  #---------------------------------------------------------------------------
  # Parameter posteriors
  #---------------------------------------------------------------------------
  
  grid  <- seq(noBurnInIterations, noBurnInIterations + nPlot - 1, 1)
  parameterNames <- c(expression(mu), expression(phi), expression(sigma[v]))
  parameterACFnames <- c(expression("ACF of " * mu), expression("ACF of " * phi), expression("ACF of " * sigma[v]))
  parameterScales <- c(-1, 1, 0.88, 1.0, 0, 0.4)
  parameterScales <- matrix(parameterScales, nrow = 3, ncol = 2, byrow = TRUE)
  parameterColors <- c("#7570B3", "#E7298A", "#66A61E")
  
  for (k in 1:3) {
    
    # Histogram of the posterior
    hist(
      resTh[, k],
      breaks = floor(sqrt(noIterations - noBurnInIterations)),
      col = rgb(t(col2rgb(parameterColors[k])) / 256, alpha = 0.25),
      border = NA,
      xlab = parameterNames[k],
      ylab = "posterior estimate",
      main = "",
      xlim = parameterScales[k,],
      freq = FALSE
    )
    
    # Add lines for the kernel density estimate of the posterior
    kde <- density(resTh[, k], kernel = "e", from = parameterScales[k, 1], to = parameterScales[k, 2])
    lines(kde, lwd = 2, col = parameterColors[k])
    
    # Plot the estimate of the posterior mean
    abline(v = thhat[k], lwd = 1, lty = "dotted")
    
    # Plot trace of the Markov chain
    plot(
      grid,
      resTh[1:nPlot, k],
      col = parameterColors[k],
      type = "l",
      xlab = "iteration",
      ylab = parameterNames[k],
      ylim = parameterScales[k,],
      xlim = c(2500, 4000),
      bty = "n"
    )
    polygon(
      c(grid, rev(grid)),
      c(resTh[1:nPlot, k], rep(-1, length(grid))),
      border = NA,
      col = rgb(t(col2rgb(parameterColors[k])) / 256, alpha = 0.25)
    )
    abline(h = thhat[k], lwd = 1, lty = "dotted")
    
    # Plot the autocorrelation function
    acf_res <- acf(resTh[, k], plot = FALSE, lag.max = 100)
    plot(
      acf_res$lag,
      acf_res$acf,
      col = parameterColors[k],
      type = "l",
      xlab = "iteration",
      ylab = parameterACFnames[k],
      lwd = 2,
      ylim = c(-0.2, 1),
      bty = "n"
    )
    polygon(
      c(acf_res$lag, rev(acf_res$lag)),
      c(acf_res$acf, rep(0, length(acf_res$lag))),
      border = NA,
      col = rgb(t(col2rgb(parameterColors[k])) / 256, alpha = 0.25)
    )
    abline(h = 1.96 / sqrt(noIterations - noBurnInIterations), lty = "dotted")
    abline(h = -1.96 / sqrt(noIterations - noBurnInIterations), lty = "dotted")
  }
}


##############################################################################
# End of file
##############################################################################