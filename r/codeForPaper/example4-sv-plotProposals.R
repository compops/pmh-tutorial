##############################################################################
# Ugly code to plot the estimate of the posterior distribution and the 
# proposal distribution adapted from a pilot run of particle 
# Metropolis-Hastings.
# (c) Johan Dahlin 2017 under MIT license <liu@johandahlin.com.nospam>
##############################################################################

# Import helpers
library("MASS")
library("mvtnorm")

# Save plot to file
savePlotToFile <- FALSE

# Load the run 
load("../savedWorkspaces/example3-sv.RData")

##############################################################################
# Parameter proposals
##############################################################################
# The unadapted proposal
stepSize1 <- diag(c(0.10, 0.01, 0.05) ^ 2)

# The adapted proposal
stepSize2 <- matrix(
  c(
    0.137255431,-0.0016258103,
    0.0015047492,-0.0016258103,
    0.0004802053,-0.0009973058,
    0.0015047492,-0.0009973058,
    0.0031307062
  ),
  ncol = 3,
  nrow = 3
)
stepSize2 <- 0.8 * stepSize2

##############################################################################
# Create grids
##############################################################################
# Estimate the posterior mean and covariance
resTh <- res$theta[noBurnInIterations:noIterations, ]
estThe <- colMeans(resTh)
estCov <- var(resTh)

# Create a grid for each parameter
gridth1 <- seq(-1, 1, 0.01)
gridth2 <- seq(0.90, 1.05, 0.01)
gridth3 <- seq(0.01, 0.35, 0.01)

#------------------------------------------------------------------------------
# Make a grid of all pairs of parameters
#------------------------------------------------------------------------------

grid1 <- matrix(0, length(gridth1) * length(gridth2), 2)
grid2 <- matrix(0, length(gridth1) * length(gridth3), 2)
grid3 <- matrix(0, length(gridth2) * length(gridth3), 2)

kk = 1
for (ii in 1:length(gridth1)) {
  for (jj in 1:length(gridth2)) {
    grid1[kk, ] <- c(gridth1[ii], gridth2[jj])
    kk <- kk + 1
  }
}

kk = 1
for (ii in 1:length(gridth1)) {
  for (jj in 1:length(gridth3)) {
    grid2[kk, ] <- c(gridth1[ii], gridth3[jj])
    kk <- kk + 1
  }
}

kk = 1
for (ii in 1:length(gridth2)) {
  for (jj in 1:length(gridth3)) {
    grid3[kk, ] <- c(gridth2[ii], gridth3[jj])
    kk <- kk + 1
  }
}


##############################################################################
# Evaluate the proposal distribution over the grid centered at the 
# posterior mean
##############################################################################

dgrid1 <- matrix(
  dmvnorm(grid1, mean = estThe[-3], sigma = stepSize1[-3, -3]),
  length(gridth1),
  length(gridth2),
  byrow = TRUE
)

dgrid2 <- matrix(
  dmvnorm(grid2, mean = estThe[-2], sigma = stepSize1[-2, -2]),
  length(gridth1),
  length(gridth3),
  byrow = TRUE
)

dgrid3 <- matrix(
  dmvnorm(grid3, mean = estThe[-1], sigma = stepSize1[-1, -1]),
  length(gridth2),
  length(gridth3),
  byrow = TRUE
)

dgrid4 <- matrix(
  dmvnorm(grid1, mean = estThe[-3], sigma = stepSize2[-3, -3]),
  length(gridth1),
  length(gridth2),
  byrow = TRUE
)

dgrid5 <- matrix(
  dmvnorm(grid2, mean = estThe[-2], sigma = stepSize2[-2, -2]),
  length(gridth1),
  length(gridth3),
  byrow = TRUE
)

dgrid6 <- matrix(
  dmvnorm(grid3, mean = estThe[-1], sigma = stepSize2[-1, -1]),
  length(gridth2),
  length(gridth3),
  byrow = TRUE
)


##############################################################################
# Compute the 2-dimensional kernel density estimate of the posterior
##############################################################################

foo1 <- kde2d(resTh[, 1], resTh[, 2], n = 50)
foo2 <- kde2d(resTh[, 1], resTh[, 3], n = 50)
foo3 <- kde2d(resTh[, 2], resTh[, 3], n = 50)


##############################################################################
# Greate the plot
##############################################################################

if (savePlotToFile) {
cairo_pdf("../figures/example4-sv-plotProposals.pdf",
          height = 6,
          width = 8)
}

layout(matrix(1:6, 2, 3, byrow = TRUE))
par(mar = c(4, 5, 0, 0))

#------------------------------------------------------------------------------
# Mu versus phi (old proposal)
#------------------------------------------------------------------------------

contour(
  foo1,
  xlim = c(-1, 1),
  ylim = c(0.88, 1.05),
  labels = NULL,
  bty = "n",
  col = "#7570B3",
  lwd = 1.5,
  labcex = 0.001,
  xlab = expression(mu),
  ylab = expression(phi)
)

contour(
  gridth1,
  gridth2,
  dgrid1,
  labels = NULL,
  nlevels = 5,
  add = T,
  col = "grey20",
  labcex = 0.001,
  lwd = 2
)

#------------------------------------------------------------------------------
# Mu versus sigma_v (old proposal)
#------------------------------------------------------------------------------

contour(
  foo2,
  xlim = c(-1, 1),
  ylim = c(0.00, 0.35),
  labels = NULL,
  bty = "n",
  col = "#E7298A",
  lwd = 1.5,
  labcex = 0.001,
  xlab = expression(mu),
  ylab = expression(sigma[v])
)

contour(
  gridth1,
  gridth3,
  dgrid2,
  labels = NULL,
  nlevels = 5,
  add = T,
  col = "grey20",
  labcex = 0.001,
  lwd = 2
)

#------------------------------------------------------------------------------
# Phi versus sigma_v (old proposal)
#------------------------------------------------------------------------------

contour(
  foo3,
  xlim = c(0.88, 1.05),
  ylim = c(0.00, 0.35),
  labels = NULL,
  bty = "n",
  col = "#66A61E",
  lwd = 1.5,
  labcex = 0.001,
  xlab = expression(phi),
  ylab = expression(sigma[v])
)

contour(
  gridth2,
  gridth3,
  dgrid3,
  labels = NULL,
  nlevels = 5,
  add = T,
  col = "grey20",
  labcex = 0.001,
  lwd = 2
)

#------------------------------------------------------------------------------
# Mu versus phi (new proposal)
#------------------------------------------------------------------------------
contour(
  foo1,
  xlim = c(-1, 1),
  ylim = c(0.88, 1.05),
  labels = NULL,
  bty = "n",
  col = "#7570B3",
  lwd = 1.5,
  labcex = 0.001,
  xlab = expression(mu),
  ylab = expression(phi)
)

contour(
  gridth1,
  gridth2,
  dgrid4,
  labels = NULL,
  nlevels = 5,
  add = T,
  col = "grey20",
  labcex = 0.001,
  lwd = 2
)

#------------------------------------------------------------------------------
# Mu versus sigma_v (new proposal)
#------------------------------------------------------------------------------

contour(
  foo2,
  xlim = c(-1, 1),
  ylim = c(0.00, 0.35),
  labels = NULL,
  bty = "n",
  col = "#E7298A",
  lwd = 1.5,
  labcex = 0.001,
  xlab = expression(mu),
  ylab = expression(sigma[v])
)

contour(
  gridth1,
  gridth3,
  dgrid5,
  labels = NULL,
  nlevels = 5,
  add = T,
  col = "grey20",
  labcex = 0.001,
  lwd = 2
)

#------------------------------------------------------------------------------
# Phi versus sigma_v (new proposal)
#------------------------------------------------------------------------------

contour(
  foo3,
  xlim = c(0.88, 1.05),
  ylim = c(0.00, 0.35),
  labels = NULL,
  bty = "n",
  col = "#66A61E",
  lwd = 1.5,
  labcex = 0.001,
  xlab = expression(phi),
  ylab = expression(sigma[v])
)

contour(
  gridth2,
  gridth3,
  dgrid6,
  labels = NULL,
  nlevels = 5,
  add = T,
  col = "grey20",
  labcex = 0.001,
  lwd = 2
)

if (savePlotToFile) {
  dev.off()
}