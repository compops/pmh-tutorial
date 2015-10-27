##############################################################################
#
# Example of particle Metropolis-Hastings
# in a linear Gaussian state space model
#
# (c) 2015 Johan Dahlin
# johan.dahlin (at) liu.se
#
##############################################################################

# Import some libraries
import matplotlib.pylab          as plt
import numpy                     as np
import stateEstimationHelper     as helpState
import parameterEstimationHelper as helpParam

# Set the random seed
np.random.seed( 10 );

##############################################################################
# Define the model
##############################################################################

# Here, we use the following model
#
# x[tt+1] = phi   * x[tt] + sigmav * v[tt]
# y[tt]   = x[tt]         + sigmae * e[tt]
#
# where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

# Set the parameters of the model (par[0],par[1],par[2])
par = np.zeros(3)
par[0] = 0.75
par[1] = 1.00
par[2] = 1.00

# Set the number of time steps to simulate
T      = 250;

# Set the initial state
x0     = 0;

##############################################################################
# Generate data
##############################################################################

(x, y) = helpState.generateData(par, T, x0)

##############################################################################
# Parameter estimation using PMH
##############################################################################

# The inital guess of the parameter
initPar  = 0.50;

# No. particles in the particle filter ( choose nPart ~ T )
nPart    = 500;

# The length of the burn-in and the no. iterations of PMH ( nBurnIn < nRuns )
nBurnIn  = 1000;
nRuns    = 5000;

# The standard deviation in the random walk proposal
stepSize = 0.10;

# Run the PMH algorithm
res = helpParam.pmh(y,initPar,par,nPart,T,x0,helpState.pf,nRuns,stepSize)


##############################################################################
# Plot the results
##############################################################################

# Plot the parameter posterior estimate
# Solid black line indicate posterior mean
plt.subplot(2,1,1);
plt.hist(res[nBurnIn:nRuns], np.floor(np.sqrt(nRuns-nBurnIn)), normed=1, facecolor='#7570B3')
plt.xlabel("par[0]"); plt.ylabel("posterior density estimate")
plt.axvline( np.mean(res[nBurnIn:nRuns]), linewidth=2, color = 'k' )

# Plot the trace of the Markov chain after burn-in
# Solid black line indicate posterior mean
plt.subplot(2,1,2);
plt.plot(np.arange(nBurnIn,nRuns,1),res[nBurnIn:nRuns],color = '#E7298A')
plt.xlabel("iteration"); plt.ylabel("par[0]")
plt.axhline( np.mean(res[nBurnIn:nRuns]), linewidth=2, color = 'k' )

##############################################################################
# End of file
##############################################################################
