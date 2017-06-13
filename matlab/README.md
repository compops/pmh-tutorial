# MATLAB code for PMH tutorial

This MATLAB code implements the Kalman filter (KF), particle filter (PF) and particle Metropolis-Hastings (PMH) algorithm for two different dynamical models: a linear Gaussian state-space (LGSS) model and a stochastic volatilty (SV) model. Note that the Kalman filter can only be employed for the first of these two models. The details of the code is described in the tutorial paper available at: < http://arxiv.org/pdf/1511.01707 >.

Note that the MATLAB code in this folder covers the basic implementations in the paper. See the R code in r/ for all the implementations and to recreate the results in the tutorial.

Requirements
--------------
The code is written and tested for MATLAB 2016b and makes use of the statistics toolbox and the Quandl package. See < https://github.com/quandl/Matlab > for more installation and to download the toolbox. Note that urlread2 is required by the Quandl toolbox and should be installed as detailed in the README file of the Quandl toolbox.

Main script files
--------------
These are the main script files that implement the various algorithms discussed in the tutorial.

**example1_lgss.m** State estimation in a LGSS model using the KM and a fully-adapted PF (faPF). The code is discussed in Section 3.1 and the results are presented in Section 3.2 as Figure 4 and Table 1.

**example2_lgss.m** Parameter estimation of one parameter in the LGSS model using PMH with the faPF as the likelihood estimator. The code is discussed in Section 4.1 and the results are presented in Section 4.2 as Figure 5.

**example3_sv.m** Parameter estimation of three parameters in the SV model using PMH with the bootstrap PF as the likelihood estimator. The code is discussed in Section 5.1 and the results are presented in Section 5.2 as Figure 6. The code takes about an hour to run.

Supporting files
--------------
**generateData.m** Implements data generation for the LGSS model.
**kalmanFilter.m** Implements the Kalman filter for the LGSS model.
**particleFilter.m** Implements the faPF for the LGSS model.
**particleFilterSVmodel.m** Implements the bPF for the SV model.
**particleMetropolisHastings.m** Implements the PMH algorithm for the LGSS model.
**particleMetropolisHastingsSVmodel.m** Implements the PMH algorithm for the SV model.

Adapting the code for another model
--------------
See the discussion in *README.MD* in the directory *r/*.