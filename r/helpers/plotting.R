##############################################################################
# Make plots for tutorial
#
# Johan Dahlin <uni (at) johandahlin.com.nospam>
# Documentation at https://github.com/compops/pmh-tutorial
# Published under GNU General Public License
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
  iact <- c()

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

    # Add lines for prior
    prior_grid <- seq(parameterScales[k, 1], parameterScales[k, 2], 0.01)
    if (k==1) {prior_values = dnorm(prior_grid, 0, 1)}
    if (k==2) {prior_values = dnorm(prior_grid, 0.95, 0.05)}
    if (k==3) {prior_values = dgamma(prior_grid, 2, 10)}
    lines(prior_grid, prior_values, col = "darkgrey")

    # Plot trace of the Markov chain
    plot(
      grid,
      resTh[1:nPlot, k],
      col = parameterColors[k],
      type = "l",
      xlab = "iteration",
      ylab = parameterNames[k],
      ylim = parameterScales[k,],
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

    iact <- c(iact, 1 + 2 * sum(acf_res$acf))
  }

  iact
}