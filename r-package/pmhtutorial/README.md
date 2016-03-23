## Purpose
This package provides some minimal working examples for implementing 
the particle Metropolis-Hastings (PMH) algorithm for parameter inference 
in nonlinear state space models. The package accompanies a tutorial:

Dahlin, J. & Schön, T. B. "Getting started with particle Metropolis-Hastings for 
inference in nonlinear dynamical models." pre-print, arXiv:1511.01707, 2015.

Currently available at:

http://arxiv.org/pdf/1511.01707

## Usage
The main functions of the package are the five examples connected to the paper:

* ** example1_lgss() ** Demostrates the particle filter for estimating the
state in a linear Gaussian state space model.

* ** example2_lgss() ** Demostrates PMH for estimating the parameters in a 
linear Gaussian state space model.

* ** example3_sv() ** Demostrates PMH for estimating the parameters in a 
stochastic volatility model.

* ** example4_sv() ** Demostrates PMH for estimating the parameters in a 
stochastic volatility model using a tailored proposal distribution.

* ** example5_sv() ** Demostrates PMH for estimating the parameters in a 
stochastic volatility model using a reparameterised model.

## Simple example
The examples can be executed by e.g.

example2_lgss()

which will recreate on of the plots in the aforementioned tutorial. See
the reference for more background and details.

## How do I get it?
The package is available through CRAN and can be installed by

install.packages("pmhtutorial")
