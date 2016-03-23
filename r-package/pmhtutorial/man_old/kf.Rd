\name{kf}
\alias{kf}
\title{
Kalman filter for state estimate in a linear Gaussian state space model
}
\description{
Estimates the filtered state and the log-likelihood for a linear Gaussian state 
space model of the form \eqn{ x_{t} = \phi x_{t-1} + \sigma_v v_t } and \eqn{ y_t = x_t + \sigma_e e_t }, where 
\eqn{v_t} and \eqn{e_t} denote independent standard Gaussian random variables, i.e.
\eqn{N(0,1)}.
}
\usage{
kf(y, phi, sigmav, sigmae, x0, P0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
observations from the model for \eqn{t=1,...,T}.
}
  \item{phi}{
The persistence of the state process \eqn{\phi}.
}
  \item{sigmav}{
The standard deviation of the state process \eqn{\sigma_v}.
}
  \item{sigmae}{
The standard deviation of the observation process \eqn{\sigma_e}.
}
  \item{x0}{
The initial state.
}
  \item{P0}{
The initial covariance of the state.
}
}
\value{
\item{xh}{The estimate of the filtered state at time \eqn{t=1,...,T}.}
\item{ll}{The estimate of the log-likelihood.}
}
\references{
Dahlin, J. & Schoen, T. B. "Getting started with particle Metropolis-Hastings for 
inference in nonlinear dynamical models." pre-print, arXiv:1511.01707, 2015.
}
\author{
Johan Dahlin <johan.dahlin@liu.se>
}
\note{
See Section 4 in the reference above for more details.
}

\examples{
## Generates T=500 observations from a linear state space model with
## ( phi, sigma_e, sigma_v ) = ( 0.5, 1.0, 0.1 ) and zero initial state.

d     <- generateData(phi=0.5, sigmav=1.0, sigmae=0.1, T=500, x0=0.0) 

## Estimate the filtered state
xhatf <- kf(d$y, phi=0.5, sigmav=1.0, sigmae=0.1, x0=0.0, P0=0.01)$xh

## Plot the estimate and the true state
par(mfrow=c(1,1))
plot(d$x, type="l", xlab="time", ylab="state")
lines(xhatf, lty="dashed",col="red")
}
\keyword{ts}
