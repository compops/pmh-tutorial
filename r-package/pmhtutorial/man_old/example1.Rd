\name{example1_lgss}
\alias{example1_lgss}
\title{State estimation in a linear Gaussian state space model}
\description{
Minimal working example of state estimation in a linear Gaussian state space 
model using Kalman filtering and a fully-adapted particle filter. The code 
estimates the bias and mean squared error (compared with the Kalman estimate)
while varying the number of particles in the particle filter.
}
\usage{
example1_lgss()
}
\details{
The Kalman filter is a standard implementation without an input. The particle 
filter is fully adapted (i.e. takes the current observation into account when
proposing new particles and computing the weights).
}
\value{
Returns a plot with the generated observations y and the difference in the
state estimates obtained by the Kalman filter (the optimal solution) and the
particle filter (with 20 particles). Furthermore, the function returns plots
of the estimated bias and mean squared error of the state estimate obtained 
using the particle filter (while varying the number of particles) and the Kalman
estimates.
\item{y}{The observations generated from the model.}
\item{x}{The states generated from the model.}
\item{xhatKF}{The estimate of the state from the Kalman filter.}
\item{xhatPF}{The estimate of the state from the particle filter with 20 particles.}
}
\references{
Dahlin, J. & Schoen, T. B. "Getting started with particle Metropolis-Hastings for 
inference in nonlinear dynamical models." pre-print, arXiv:1511.01707, 2015.
}
\author{
Johan Dahlin <johan.dahlin@liu.se>
}
\note{
See Section 4.2 in the reference above for more details.
}
\keyword{misc}
