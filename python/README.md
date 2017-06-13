# Python code for PMH tutorial

This R code implements the Kalman filter (KF), particle filter (PF) and particle Metropolis-Hastings (PMH) algorithm for two different dynamical models: a linear Gaussian state-space (LGSS) model and a stochastic volatilty (SV) model. Note that the Kalman filter can only be employed for the first of these two models. The details of the code is described in the tutorial paper available at: < http://arxiv.org/pdf/1511.01707 >.

Note that the Python code in this folder covers the basic implementations in the paper. See the R code in r/ for all the implementations and to recreate the results in the tutorial.

Requirements
--------------
The code is written and tested for Python 2.7.6/3.6.0 togehter with NumPy 1.9.2/1.11.3, SciPy 0.15.1/0.18.1, Matplotlib 1.4.3/2.0.0 and Quandl 2.8.9/3.1.0. These packages are easily avaiable via Anaconda (https://docs.continuum.io/anaconda/install) by installing the package for your preference of Python version and then executing
``` bash
conda install numpy scipy matplotlib quandl
```
For more information about the Quandl library, see < https://www.quandl.com/tools/python >.

Main script files
--------------
These are the main script files that implement the various algorithms discussed in the tutorial.

**example1-lgss.py** State estimation in a LGSS model using the KM and a fully-adapted PF (faPF). The code is discussed in Section 3.1 and the results are presented in Section 3.2 as Figure 4 and Table 1.

**example2-lgss.py** Parameter estimation of one parameter in the LGSS model using PMH with the faPF as the likelihood estimator. The code is discussed in Section 4.1 and the results are presented in Section 4.2 as Figure 5.

**example3-sv.py** Parameter estimation of three parameters in the SV model using PMH with the bootstrap PF as the likelihood estimator. The code is discussed in Section 5.1 and the results are presented in Section 5.2 as Figure 6. The code takes about an hour to run.

Supporting files (helpers/)
--------------
**stateEstimation.py**
Implements data generation for the LGSS model (generateData), the faPF for the LGSS model (particleFilter), the Kalman filter for the LGSS model (kalmanFilter) and the bPF for the SV model (paticleFilterSVmodel).

**parameterEstimation.py**
Implements the PMH algorithm for the LGSS model (particleMetropolisHastings) and the SV model (particleMetropolisHastingsSVModel). 

Adapting the code for another model
--------------
See the discussion in *README.MD* in the directory *r/*.