\name{pmhtutorial-package}
\alias{pmhtutorial-package}
\alias{pmhtutorial}
\docType{package}
\title{
Minimal working examples for particle Metropolis-Hastings
}
\description{
The package includes routines for state estimate in a linear 
Gaussian state space model and a simple stochastic volatility model using
particle filtering. Parameter inference is also carried out in these models
using the particle Metropolis-Hastings algorithm that includes the particle
filter to provied an unbiased estimator of the likelihood. The code is developed
in the paper and all the details are explained their. This package is meant as
a minimal working example of these algorithm and is only meant for educational
use and as a start for learning to implement these algorithms on your own.
}
\details{
See the paper in the references below for background information about the 
package. The main examples in the paper are found in the five example commands
and these replicates the results in the paper and are meant for educational 
purposes to teach the reader about how to implement particle filtering and 
particle Metropolis-Hastings.
}
\author{
Johan Dahlin <johan.dahlin@liu.se>

Maintainer: Johan Dahlin <johan.dahlin@liu.se>
}
\references{
Dahlin, J. & Schoen, T. B. "Getting started with particle Metropolis-Hastings for 
inference in nonlinear dynamical models." pre-print, arXiv:1511.01707, 2015.
}
\keyword{ package }
\seealso{\code{\link{example1_lgss}}, \code{\link{example2_lgss}}, 
\code{\link{example3_sv}},\code{\link{example4_sv}},
\code{\link{example5_sv}}}