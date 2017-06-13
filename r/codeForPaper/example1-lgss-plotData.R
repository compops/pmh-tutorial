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
source("../helpers/stateEstimation.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

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

##############################################################################
# Plotting
##############################################################################

# Export plot to file
if (savePlotToFile) {
  cairo_pdf("../figures/lgss-data.pdf",
            height = 3,
            width = 8)
}

grid = seq(0, T)

# Plot the latent state and observations
layout(matrix(1:3, 1, 3, byrow = TRUE))
par(mar = c(4, 5, 0, 0))

plot(
  grid,
  x,
  col = "#D95F02",
  lwd = 1,
  type = "l",
  xlab = "time",
  ylab = expression("latent state " * x[t]),
  bty = "n",
  ylim = c(-4, 6)
)
polygon(c(grid, rev(grid)),
        c(x, rep(-4, T + 1)),
        border = NA,
        col = rgb(t(col2rgb("#D95F02")) / 256, alpha = 0.25))

plot(
  grid[-1],
  y[-1],
  col = "#1B9E77",
  lwd = 1,
  type = "l",
  xlab = "time",
  ylab = expression("observation " * y[t]),
  bty = "n",
  ylim = c(-4, 6)
)
polygon(c(grid[-1], rev(grid[-1])),
        c(y[-1], rep(-4, T)),
        border = NA,
        col = rgb(t(col2rgb("#1B9E77")) / 256, alpha = 0.25))

foo = acf(y[-1], plot = F, lag.max = 25)

plot(
  foo$lag,
  foo$acf,
  col = "#66A61E",
  lwd = 1.5,
  type = "l",
  xlab = "time",
  ylab = expression("ACF of " * y[t]),
  bty = "n",
  ylim = c(-0.2, 1),
  xlim = c(0, 25)
)
polygon(
  c(foo$lag, rev(foo$lag)),
  c(foo$acf, rep(0.0, length(foo$lag))),
  border = NA,
  col = rgb(t(col2rgb("#66A61E")) / 256, alpha = 0.25)
)
abline(h = -1.96 / sqrt(T), lty = "dotted")
abline(h = 1.96 / sqrt(T), lty = "dotted")


if (savePlotToFile) {
  dev.off()
}

###################################################################################
# End of file
###################################################################################
