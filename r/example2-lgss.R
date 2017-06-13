##############################################################################
#
# Example of particle Metropolis-Hastings
# in a linear Gaussian state space model
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

# Import helpers
source("helpers/stateEstimation.R")
source("helpers/parameterEstimation.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

# Should the results be loaded from file (to quickly generate plots)
loadSavedWorkspace <- FALSE

# Save plot to file
savePlotToFile <- TRUE


##############################################################################
# Define the model
##############################################################################

# Here, we use the following model
#
# x[tt+1] = phi   * x[tt] + sigmav * v[tt]
# y[tt]   = x[tt]         + sigmae * e[tt]
#
# where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

# Set the parameters of the model
phi <- 0.75
sigmav <- 1.00
sigmae <- 0.10

# Set the number of time steps to simulate
T <- 250

# Set the initial state
initialState <- 0


##############################################################################
# Generate data
##############################################################################

data <- generateData(c(phi, sigmav, sigmae), T, initialState)


##############################################################################
# Parameter estimation using PMH
##############################################################################

# The inital guess of the parameter
initialPhi <- 0.50

# No. particles in the particle filter
noParticles <- 100

# The length of the burn-in and the no. iterations of PMH ( noBurnInIterations < noIterations )
noBurnInIterations <- 1000
noIterations <- 5000

# Run the PMH algorithm
if (loadSavedWorkspace) {
  load("savedWorkspaces/example2-lgss.RData")
} else {
  res1 <-
    particleMetropolisHastings(
      data$y,
      initialPhi,
      sigmav,
      sigmae,
      noParticles,
      initialState,
      noIterations,
      stepSize = 0.01
    )
  res2 <-
    particleMetropolisHastings(
      data$y,
      initialPhi,
      sigmav,
      sigmae,
      noParticles,
      initialState,
      noIterations,
      stepSize = 0.10
    )
  res3 <-
    particleMetropolisHastings(
      data$y,
      initialPhi,
      sigmav,
      sigmae,
      noParticles,
      initialState,
      noIterations,
      stepSize = 0.50
    )
}

##############################################################################
# Plot the results
##############################################################################

# Extract the states after burn-in
resTh1 <- res1[noBurnInIterations:noIterations,]
resTh2 <- res2[noBurnInIterations:noIterations,]
resTh3 <- res3[noBurnInIterations:noIterations,]

# Estimate the KDE of the marginal posteriors
kde1  <- density(resTh1,
                 kernel = "e",
                 from = 0.5,
                 to = 0.8)
kde2  <- density(resTh2,
                 kernel = "e",
                 from = 0.5,
                 to = 0.8)
kde3  <- density(resTh3,
                 kernel = "e",
                 from = 0.5,
                 to = 0.8)

# Export plot to file
if (savePlotToFile) {
  cairo_pdf("figures/example2-lgss.pdf",
            height = 10,
            width = 8)
}

layout(matrix(1:9, 3, 3, byrow = TRUE))
par   (mar = c(4, 5, 0, 0))

# Plot the parameter posterior estimate
hist(
  resTh1,
  breaks = floor(sqrt(noIterations - noBurnInIterations)),
  col = rgb(t(col2rgb("#7570B3")) / 256, alpha = 0.25),
  border = NA,
  xlab = expression(phi),
  ylab = "posterior estimate",
  main = "",
  xlim = c(0.5, 0.8),
  ylim = c(0, 12),
  freq = FALSE
)
lines(kde1, lwd = 2, col = "#7570B3")
abline(v = mean(resTh1),
       lwd = 1,
       lty = "dotted")

hist(
  resTh2,
  breaks = floor(sqrt(noIterations - noBurnInIterations)),
  col = rgb(t(col2rgb("#E7298A")) / 256, alpha = 0.25),
  border = NA,
  xlab = expression(phi),
  ylab = "posterior estimate",
  main = "",
  xlim = c(0.5, 0.8),
  ylim = c(0, 12),
  freq = FALSE
)
lines(kde2, lwd = 2, col = "#E7298A")
abline(v = mean(resTh2),
       lwd = 1,
       lty = "dotted")

hist(
  resTh3,
  breaks = floor(sqrt(noIterations - noBurnInIterations)),
  col = rgb(t(col2rgb("#66A61E")) / 256, alpha = 0.25),
  border = NA,
  xlab = expression(phi),
  ylab = "posterior estimate",
  main = "",
  xlim = c(0.5, 0.8),
  ylim = c(0, 12),
  freq = FALSE
)
lines(kde3, lwd = 2, col = "#66A61E")
abline(v = mean(resTh3),
       lwd = 1,
       lty = "dotted")

# Plot the trace of the Markov chain during 1000 iterations after the burn-in
grid <- seq(noBurnInIterations, noBurnInIterations + 1000 - 1, 1)

plot(
  grid,
  resTh1[1:1000],
  col = '#7570B3',
  type = "l",
  xlab = "iteration",
  ylab = expression(phi),
  ylim = c(0.4, 0.8),
  bty = "n"
)
abline(h = mean(resTh1),
       lwd = 1,
       lty = "dotted")
polygon(
  c(grid, rev(grid)),
  c(resTh1[1:1000], rep(0.4, 1000)),
  border = NA,
  col = rgb(t(col2rgb("#7570B3")) / 256, alpha = 0.25)
)

plot(
  grid,
  resTh2[1:1000],
  col = '#E7298A',
  type = "l",
  xlab = "iteration",
  ylab = expression(phi),
  ylim = c(0.4, 0.8),
  bty = "n"
)
abline(h = mean(resTh2),
       lwd = 1,
       lty = "dotted")
polygon(
  c(grid, rev(grid)),
  c(resTh2[1:1000], rep(0.4, 1000)),
  border = NA,
  col = rgb(t(col2rgb("#E7298A")) / 256, alpha = 0.25)
)

plot(
  grid,
  resTh3[1:1000],
  col = '#66A61E',
  type = "l",
  xlab = "iteration",
  ylab = expression(phi),
  ylim = c(0.4, 0.8),
  bty = "n"
)
abline(h = mean(resTh3),
       lwd = 1,
       lty = "dotted")
polygon(
  c(grid, rev(grid)),
  c(resTh3[1:1000], rep(0.4, 1000)),
  border = NA,
  col = rgb(t(col2rgb("#66A61E")) / 256, alpha = 0.25)
)

# Plot the ACF of the Markov chain

res1ACF <- acf(resTh1, plot = FALSE, lag.max = 60)
plot(
  res1ACF$lag,
  res1ACF$acf,
  col = '#7570B3',
  type = "l",
  xlab = "iteration",
  ylab = "ACF",
  ylim = c(-0.2, 1),
  bty = "n"
)
polygon(
  c(res1ACF$lag, rev(res1ACF$lag)),
  c(res1ACF$acf, rep(0, length(res1ACF$lag))),
  border = NA,
  col = rgb(t(col2rgb("#7570B3")) / 256, alpha = 0.25)
)
abline(h = 1.96 / sqrt(length(grid)), lty = "dotted")
abline(h = -1.96 / sqrt(length(grid)), lty = "dotted")

res2ACF <- acf(resTh2, plot = FALSE, lag.max = 60)
plot(
  res2ACF$lag,
  res2ACF$acf,
  col = '#E7298A',
  type = "l",
  xlab = "iteration",
  ylab = "ACF",
  ylim = c(-0.2, 1),
  bty = "n"
)
polygon(
  c(res2ACF$lag, rev(res2ACF$lag)),
  c(res2ACF$acf, rep(0, length(res2ACF$lag))),
  border = NA,
  col = rgb(t(col2rgb("#E7298A")) / 256, alpha = 0.25)
)
abline(h = 1.96 / sqrt(length(grid)), lty = "dotted")
abline(h = -1.96 / sqrt(length(grid)), lty = "dotted")

res3ACF <- acf(resTh3, plot = FALSE, lag.max = 60)
plot(
  res3ACF$lag,
  res3ACF$acf,
  col = '#66A61E',
  type = "l",
  xlab = "iteration",
  ylab = "ACF",
  ylim = c(-0.2, 1),
  bty = "n"
)
polygon(
  c(res3ACF$lag, rev(res3ACF$lag)),
  c(res3ACF$acf, rep(0, length(res3ACF$lag))),
  border = NA,
  col = rgb(t(col2rgb("#66A61E")) / 256, alpha = 0.25)
)
abline(h = 1.96 / sqrt(length(grid)), lty = "dotted")
abline(h = -1.96 / sqrt(length(grid)), lty = "dotted")

# Close the plotting device
if (savePlotToFile) {
  dev.off()
}


##############################################################################
# Compute and save the results
##############################################################################

# Estimate the parameter posterior mean
mean(res1[grid])
mean(res2[grid])
mean(res3[grid])

# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("savedWorkspaces/example2-lgss.RData")
}


##############################################################################
# End of file
##############################################################################
