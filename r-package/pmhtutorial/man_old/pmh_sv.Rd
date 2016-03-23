\name{pmh_sv}
\alias{pmh_sv}
\title{
Particle Metropolis-Hastings algorithm for a linear Gaussian state space model
}
\description{
Estimates the parameter posterior for \eqn{\theta=\{\mu,\phi,\sigma_v\}} in a 
stochastic volatility model of the form \eqn{x_t = \mu + \phi ( x_{t-1} - \mu ) 
+ \sigma_v v_t} and \eqn{y_t = \exp(x_t/2) e_t}, where \eqn{v_t} and \eqn{e_t} denote 
independent standard Gaussian random variables, i.e. \eqn{N(0,1)}.
}
\usage{
pmh_sv(y, initPar, nPart, T, nIter, stepSize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
observations from the model for \eqn{t=1,...,T}.
}
  \item{initPar}{
An inital value for the parameters \eqn{\theta=\{\mu,\phi,\sigma_v\}}.
}
  \item{nPart}{
The number of particles in the particle filter.
}
\item{T}{
The number of observations.
}
  \item{nIter}{
The number of iterations in the PMH algorithm.
}
  \item{stepSize}{
The covariance matrix for the Gaussian random walk proposal for \eqn{\theta}.
}
}

\value{
The trace of the Markov chain exploring the posterior of \eqn{\theta}.
}
\references{
Dahlin, J. & Schoen, T. B. "Getting started with particle Metropolis-Hastings for 
inference in nonlinear dynamical models." pre-print, arXiv:1511.01707, 2015.
}
\author{
Johan Dahlin <johan.dahlin@liu.se>
}
\note{
See Section 5 in the reference above for more details.
}

\examples{
# Use a short run when testing the code on CRAN
# Set nIter manually if running on local computer, 
# nIter=1000 is usually enough
if(nchar(Sys.getenv("LONG_TEST")) != 0) {
  nIter=7500
} else {
  nIter=5
}

## Get the data
library("Quandl")
d <- Quandl("NASDAQOMX/OMXS30", start_date="2012-01-02", 
            end_date="2014-01-02", type="zoo")
y <- as.numeric(100 * diff(log(d$"Index Value")))

## Estimate the filtered state
res <- pmh_sv(y, c(0, 0.9, 0.2), nPart=500, T=500, nIter, 
             stepSize=diag(c(0.05, 0.0002, 0.002)))

## Plot the estimate
par(mfrow=c(3,1))
hist(res$thhat[,1], main="", xlab=expression(mu), 
     ylab="marginal posterior", freq=FALSE)
hist(res$thhat[,2], main="", xlab=expression(phi), 
     ylab="marginal posterior", freq=FALSE)
hist(res$thhat[,3], main="", xlab=expression(sigma[v]), 
     ylab="marginal posterior", freq=FALSE)
}
\keyword{ts}
