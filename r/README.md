# R code for PMH tutorial

This R code implements the Kalman filter (KF), particle filter (PF) and particle Metropolis-Hastings (PMH) algorithm for two different dynamical models: a linear Gaussian state-space (LGSS) model and a stochastic volatilty (SV) model. Note that the Kalman filter can only be employed for the first of these two models. The details of the code is described in the tutorial paper available at: < http://arxiv.org/pdf/1511.01707 >.

Requirements
--------------
The code is written and tested for R 3.2.2 and makes use of the packages Quandl and mvtnorm. These can be installed in R by executing the command:
``` R
install.packages(c("mvtnorm", "Quandl")) 
```

Main script files
--------------
These are the main script files that implement the various algorithms discussed in the tutorial.

**example1-lgss.R** State estimation in a LGSS model using the KM and a fully-adapted PF (faPF). The code is discussed in Section 3.1 and the results are presented in Section 3.2 as Figure 4 and Table 1.

**example2-lgss.R** Parameter estimation of one parameter in the LGSS model using PMH with the faPF as the likelihood estimator. The code is discussed in Section 4.1 and the results are presented in Section 4.2 as Figure 5.

**example3-sv.R** Parameter estimation of three parameters in the SV model using PMH with the bootstrap PF as the likelihood estimator. The code is discussed in Section 5.1 and the results are presented in Section 5.2 as Figure 6. The code takes about an hour to run.

**example4-sv.R** Modified version of the code in *example3-sv.R* to make use of a better tailored parameter proposal. The details are discussed in Section 6.3.2 and the results are presented in the same section as Figures 7 and 8. Note that the only difference in the code is that the variable stepSize is changed.

**example5-sv.R** Modified version of the code in *example3-sv.R* to make use of another parameterisation of the model and a better tailored parameter proposal. The details are discussed in Section 6.3.3 and the results are presented in the same section. Note that the differences in the code is the use of another implemenation of PMH ant that the variable stepSize is changed.


Additional script files for creating plots for tutorial (codeForPaper/)
--------------
These are some additional files to recreate some extra results discussed in the tutorial.

**example1-lgss-plotData.R** Some sample code for generate data and recreate the plot of the data presented as Figure 3.

**example2-lgss-varyingT.R** An extended version of *example2-lgss.R* and makes several runs while changing the number of observations. The results are presented in Section 3.2 as Table 1.

**example4-sv-plotProposals.R** Some (ugly) code to plot the estimates of the posterior distribution and the proposal distribution using the output from a run of *example3-sv.R*. This code generates Figure 7 in Section 6.3.2.


Supporting files (helpers/)
--------------
**stateEstimation.R**
Implementes data generation for the LGSS model (generateData), the faPF for the LGSS model (particleFilter), the Kalman filter for the LGSS model (kalmanFilter) and the bPF for the SV model (paticleFilterSVmodel).

**parameterEstimation.R**
Implementes the PMH algorithm for the LGSS model (particleMetropolisHastings), the SV model (particleMetropolisHastingsSVModel) and the reparameterised SV model (particleMetropolisHastingsSVModelReparameterised). Moreover, a script (makePlotsParticleMetropolisHastingsSVModel) is included to generate the Figures presented in the paper using the output of the PMH algorithm.


Saved results (savedWorkspaces/ and figures/)
--------------
**savedWorkspaces/** Saved copies of the workspace after running the corresponding code. These outputs are used to generate all the results in the aper. Can be used to directly recreate the plots in the tutorial by setting the flags loadSavedWorkspace and savePlotToFile to TRUE.