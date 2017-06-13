##############################################################################
#
# Example of fully-adapted particle filtering
# in a linear Gaussian state space model
#
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

# Import helper
source("stateEstimationHelper.R")

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
x <- data$x
y <- data$y

# Export plot to file
if (savePlotToFile) {
  cairo_pdf("figures/example1-lgss.pdf",
            height = 10,
            width = 8)
}

# Plot the latent state and observations
layout(matrix(c(1, 1, 2, 2, 3, 4), 3, 2, byrow = TRUE))
par   (mar = c(4, 5, 0, 0))

grid <- seq(0, T)

plot(
  grid,
  y,
  col = "#1B9E77",
  type = "l",
  xlab = "time",
  ylab = "observation",
  ylim = c(-6, 6),
  bty = "n"
)
polygon(c(grid, rev(grid)),
        c(y, rep(-6, T + 1)),
        border = NA,
        col = rgb(t(col2rgb("#1B9E77")) / 256, alpha = 0.25))


###################################################################################
# State estimation using the particle filter and Kalman filter
###################################################################################

if (loadSavedWorkspace) {
  load("savedWorkspaces/example1-lgss.RData")
} else {
  # Using noParticles = 20 particles and plot the estimate of the latent state
  noParticles <- 20
  outputPF <-
    particleFilter(y, c(phi, sigmav, sigmae), noParticles, initialState)
  outputKF <-
    kalmanFilter(y, c(phi, sigmav, sigmae), initialState, 0.01)
  difference <-
    outputPF$xHatFiltered - outputKF$xHatFiltered[-(T + 1)]
}

grid <- seq(0, T - 1)
plot(
  grid,
  difference,
  col = "#7570B3",
  type = "l",
  xlab = "time",
  ylab = "difference in state estimate",
  ylim = c(-0.1, 0.1),
  bty = "n"
)
polygon(
  c(grid, rev(grid)),
  c(difference, rep(-0.1, T)),
  border = NA,
  col = rgb(t(col2rgb("#7570B3")) / 256, alpha = 0.25)
)

# Compute bias and MSE
logBiasMSE <- matrix(0, nrow = 7, ncol = 2)
gridN <- c(10, 20, 50, 100, 200, 500, 1000)

for (ii in 1:length(gridN)) {
  pfEstimate <-
    particleFilter(y, c(phi, sigmav, sigmae), gridN[ii], initialState)
  pfEstimate <- pfEstimate$xHatFiltered
  kfEstimate <- outputKF$xHatFiltered[-(T + 1)]
  
  logBiasMSE[ii, 1] <- log(mean(abs(pfEstimate - kfEstimate)))
  logBiasMSE[ii, 2] <- log(mean((pfEstimate - kfEstimate) ^ 2))
}

# Plot the bias and MSE for comparison
plot(
  gridN,
  logBiasMSE[, 1],
  col = "#E7298A",
  type = "l",
  xlab = "no. particles (N)",
  ylab = "log-bias",
  ylim = c(-7,-3),
  bty = "n"
)
polygon(
  c(gridN, rev(gridN)),
  c(logBiasMSE[, 1], rep(-7, length(gridN))),
  border = NA,
  col = rgb(t(col2rgb("#E7298A")) / 256, alpha = 0.25)
)
points(gridN, logBiasMSE[, 1], col = "#E7298A", pch = 19)

plot(
  gridN,
  logBiasMSE[, 2],
  col = "#66A61E",
  lwd = 1.5,
  type = "l",
  xlab = "no. particles (N)",
  ylab = "log-MSE",
  ylim = c(-12,-6),
  bty = "n"
)
polygon(
  c(gridN, rev(gridN)),
  c(logBiasMSE[, 2], rep(-12, length(gridN))),
  border = NA,
  col = rgb(t(col2rgb("#66A61E")) / 256, alpha = 0.25)
)
points(gridN, logBiasMSE[, 2], col = "#66A61E", pch = 19)

# Close the plotting device
if (savePlotToFile) {
  dev.off()
}


##############################################################################
# Compute and save the results
##############################################################################

# Print a table (no. particles, log-bias, log-mse)
print(t(rbind(gridN, t(logBiasMSE))))

# gridN                     
# [1,]    10 -3.696997  -6.938594
# [2,]    20 -3.964671  -7.493297
# [3,]    50 -4.567552  -8.718346
# [4,]   100 -4.850363  -9.294468
# [5,]   200 -5.192173  -9.905719
# [6,]   500 -5.668407 -10.866745
# [7,]  1000 -6.077648 -11.671646

# Save the workspace to file
if (!loadSavedWorkspace) {
  save.image("savedWorkspaces/example1-lgss.RData")
}


###################################################################################
# End of file
###################################################################################
