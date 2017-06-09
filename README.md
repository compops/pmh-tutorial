# pmh-tutorial

This code was downloaded from < https://github.com/compops/pmh-tutorial > or from < http://work.johandahlin.com/ > and contains the code used to produce the results in the tutorial

J. Dahlin and T. B. Schön, **Getting started with particle Metropolis-Hastings for inference in nonlinear models**. Pre-print, arXiv:1511:01707, 2017. 

The papers are available as a preprint from < http://arxiv.org/pdf/1511.01707 > and < http://work.johandahlin.com/ >.

Requirements
--------------
The code is written and tested for R 3.2.2, Matlab R2014a and Python 2.7.6 with some additional libraries/packages (see below).

The implementation in R makes use of the packages Quandl and mvtnorm. They can be installed by the command 
``` R
install.packages(c("mvtnorm", "Quandl")) 
```
The implementation in Python makes use of NumPy 1.9.2, SciPy 0.15.1, Matplotlib 1.4.3, Pandas 0.13.1 and Quandl 2.8.9. On Ubuntu, these packages can be installed/upgraded using
``` bash
sudo pip install --upgrade package-name
```
For more information about the Quandl library, see < https://www.quandl.com/tools/python >.

The implementation in Matlab makes use of the statistics toolbox and the Quandl package. See < https://github.com/quandl/Matlab > for more installation and to download the toolbox. Note that urlread2 is required by the Quandl toolbox and should be installed as detailed in the README file of the Quandl toolbox.

Included files (folders matlab, python and r)
--------------
**example1-lgss.[R,py,m]** Implements the numerical illustration in Section 3.2 of state estimation in a linear Gaussian state space (LGSS) model using the fully-adapted particle filter (faPF). The output is the filtered state estimated as presented in Figure 3.

**example2-lgss.[R,py,m]** Implements the numerical illustration in Section 3.4 of parameter estimation of the parameter phi in the LGSS model using particle Metropolis-Hastings (PMH) with the faPF as the likelihood estimator. The output is the estimated parameter posterior as presented in Figure 4.

**example3-sv.[R,py,m]** Implements the numerical illustration in Section 4 of parameter estimation of the three parameters in the stochastic volatility (SV) model using particle Metropolis-Hastings (PMH) with the bootstrap particle filter (bPF) as the likelihood estimator. The output is the estimated parameter posterior as presented in Figure 5. The code takes some time (hours to execute).

**example4-sv.[R,py,m]** Implements the numerical illustration in Section 5.3.1, which makes use of the same setup as in Section 4.1 but with a tuned proposal distribution. The output is the estimated ACF and IACT as presented in Figure 7. The code takes some time (hours to execute).

**example5-sv.[R]** Implements the numerical illustration in Section 5.3.2, which makes use of the same setup as in Section 4.1 but with a tuned proposal distribution. The output is the estimated ACF and IACT as presented in Figure 7. The code takes some time (hours to execute).

**example*.[RData,mat,spdata]** A saved copy of the workspace after running thr corresponding code. Can be used to directly recreate the plots in the tutorial and to conduct additional analysis.

Supporting files (folders matlab, python and r)
--------------
**stateEstimationHelper.[py,R]**
Implementes the data generation for the LGSS model (generateData), the faPF for the LGSS model (sm), the Kalman filter for the LGSS model (kf) and the bPF for the SV model (sm_sv). In Matlab, these functions are defined in four seperate m-files with the corresponding file names.

**parameterEstimationHelper.[py,R]**
Implementes the PMH algorithm for the LGSS model (pmh), the SV model (pmh_sv) and the reparameterised SV model (pmh_sv_reparametrised). In Matlab, these functions are defined in two seperate m-files with the corresponding file names.

Included files (folder r-package)
--------------
The files for the R package pmhtutorial. These should be virtually the same as the files in the folder r but packaged as an R package. The folder r is maintained to keep the code for the three languages as similar as possible. However, the r package is simple to download and use directly.

Included files (folder matlab-skeleton)
--------------
Skeleton code files for MATLAB to help step-by-step implementation during courses and seminars. 


