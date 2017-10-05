##############################################################################
# Generates and plots data from a LGSS model.
#
# Johan Dahlin <liu (at) johandahlin.com.nospam>
# Documentation at https://github.com/compops/pmh-tutorial
# Published under GNU General Public License
##############################################################################

source("../helpers/dataGeneration.R")
source("../helpers/stateEstimation.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

# Save plot to file
savePlotToFile <- TRUE

##############################################################################
# Define the model and generate data
# x[t + 1] = phi * x[t] + sigmav * v[t],    v[t] ~ N(0, 1)
# y[t] = x[t] + sigmae * e[t],              e[t] ~ N(0, 1)
##############################################################################
phi <- 0.75
sigmav <- 1.00
sigmae <- 0.10
T <- 250
initialState <- 0

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