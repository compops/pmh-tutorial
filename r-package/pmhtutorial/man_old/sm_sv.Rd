\name{sm_sv}
\alias{sm_sv}
\title{
Bootstrap particle filter for state estimate in a simple stochastic volatility 
model
}
\description{
Estimates the filtered state and the log-likelihood for a stochastic volatility 
model of the form \eqn{x_t = \mu + \phi ( x_{t-1} - \mu ) 
+ \sigma_v v_t} and \eqn{y_t = \exp(x_t/2) e_t}, where \eqn{v_t} and \eqn{e_t} denote 
independent standard Gaussian random variables, i.e. \eqn{N(0,1)}.
}
\usage{
sm_sv(y, mu, phi, sigmav, nPart, T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
observations from the model for \eqn{t=1,...,T}.
}
  \item{mu}{
The mean of the log-volatility process \eqn{\mu}.
}
  \item{phi}{
The persistence of the log-volatility process \eqn{\phi}.
}
  \item{sigmav}{
The standard deviation of the log-volatility process \eqn{\sigma_v}.
}
  \item{nPart}{
The number of particles to use in the filter.
}
  \item{T}{
The number of observations.
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
See Section 5 in the reference above for more details.
}

\examples{
## Get the data
library("Quandl")
d <- Quandl("NASDAQOMX/OMXS30", start_date="2012-01-02", 
            end_date="2014-01-02", type="zoo")
y <- as.numeric(100 * diff(log(d$"Index Value")))

## Estimate the filtered state
xhatf <- sm_sv(y, mu=-0.10, phi=0.97, sigmav=0.15, nPart=500, 
               T=500 )$xh

## Plot the data and the log-volatility estimate
par(mfrow=c(2,1))
plot(y, type="l", xlab="time", ylab="log-returns")
plot(xhatf, type="l", xlab="time", ylab="estimate of log-volatility")
}
\keyword{ts}
