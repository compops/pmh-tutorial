# Python code for PMH tutorial

This Python code implements the Kalman filter (KF), particle filter (PF) and particle Metropolis-Hastings (PMH) algorithm for two different dynamical models: a linear Gaussian state-space (LGSS) model and a stochastic volatility (SV) model. Note that the Kalman filter can only be employed for the first of these two models. The details of the code is described in [the tutorial paper](https://doi.org/10.18637/jss.v088.c02).

Note that the Python code in this folder covers the basic implementations in the paper. The notation of the variables has been changed slightly compared with the tutorial paper to improve readability of the code. However, it should be easy to translate between the two. See the R code in r/ for all the implementations and to recreate the results in the tutorial.

## Requirements
The code is written and tested for `Python 2.7.6/3.6` together with `NumPy 1.9.2/1.11.3`, `SciPy 0.15.1/0.18.1`, `Matplotlib 1.4.3/2.0.0` and `Quandl 2.8.9/3.1.0`. These packages are easily available via [Anaconda](https://docs.continuum.io/anaconda/install) by installing the package for your preference of Python version and then executing
``` bash
conda install numpy scipy matplotlib quandl
```
For more information about the Quandl library, see [the documentation](https://www.quandl.com/tools/python).

## Main script files
These are the main script files that implement the various algorithms discussed in the tutorial.

* **example1-lgss.py** State estimation in a LGSS model using the KM and a fully-adapted PF (faPF). The code is discussed in Section 3.1 and the results are presented in Section 3.2 as Figure 4 and Table 1.

* **example2-lgss.py** Parameter estimation of one parameter in the LGSS model using PMH with the faPF as the likelihood estimator. The code is discussed in Section 4.1 and the results are presented in Section 4.2 as Figure 5.

* **example3-sv.py** Parameter estimation of three parameters in the SV model using PMH with the bootstrap PF as the likelihood estimator. The code is discussed in Section 5.1 and the results are presented in Section 5.2 as Figure 6. The code takes about an hour to run.

## Supporting files (helpers/)
* **dataGeneration.py** Generates data from a LGSS model.

* **parameterEstimation.py** Implements the PMH algorithm for the LGSS model (particleMetropolisHastings) and the SV model (particleMetropolisHastingsSVModel).

* **stateEstimation.py** Implements the faPF for the LGSS model (particleFilter), the Kalman filter for the LGSS model (kalmanFilter) and the bPF for the SV model (paticleFilterSVmodel).


## Adapting the code for another model
See the discussion in *README.MD* in the directory *r/*.

## Generalisations
Some generalisations and improvements of this code is discussed in the tutorial, see the last paragraph in Section 7. Python code for PMH1 and PMH2 is available in the repo [pmh-stco2015](https://github.com/compops/pmh-stco2015), Python code for qPMH2 is availabe in the repo [https://github.com/compops/qpmh2-sysid2015](qpmh2-sysid2015) and Python code for correlated pseudo-marginal Metropolis-Hastings is available in the repo [https://github.com/compops/pmmh-correlated2015](pmmh-correlated2015). These are excellent resources for getting up to speed with the current frontier in research connected to PMH.

## Copyright information
``` R
##############################################################################
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
##############################################################################
```
