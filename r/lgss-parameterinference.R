##############################################################################
#
# Example of particle Metropolis-Hastings 
# in a linear Gaussian state space model
#
# (c) 2015 Johan Dahlin
# johan.dahlin (at) liu.se
#
##############################################################################

# Import helper
source("stateEstimationHelper.R")
source("parameterEstimationHelper.R")

# Set the random seed to replicate results in tutorial
set.seed( 10 )

##############################################################################
# Define the model
##############################################################################

# Here, we use the following model
#
# x[tt+1] = phi   * x[tt] + sigmav * v[tt]
# y[tt]   = x[tt]         + sigmae * e[tt]
#
# where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

# Set the parameters of the model
phi     = 0.75
sigmav  = 1.00
sigmae  = 1.00

# Set the number of time steps to simulate
T      = 250;

# Set the initial state
x0     = 0;


##############################################################################
# Generate data
##############################################################################

data = generateData(phi,sigmav,sigmae,T,x0)

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
res = pmh(data$y,initPar,sigmae,sigmav,nPart,T,x0,nRuns,stepSize)

##############################################################################
# Plot the results
##############################################################################

# Export plot to file
#cairo_pdf("lgss-parameter.pdf", height = 6, width = 8)

# Plot the parameter posterior estimate
# Solid black line indicate posterior mean
layout(matrix(1:2, 2, 1, byrow = TRUE)); par(mar=c(4,5,0,0))
hist(res[nBurnIn:nRuns], breaks=floor(sqrt(nRuns-nBurnIn)), col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25),border=NA,xlab=expression(phi),ylab="posterior estimate",main="",xlim=c(0.6,1),ylim=c(0,12),freq=F)
lines(density(res[nBurnIn:nRuns],kernel="e"),lwd=2,col="#7570B3")
abline(v=mean(res[nBurnIn:nRuns]),lwd=1,lty="dotted")

# Plot the trace of the Markov chain after burn-in
# Solid black line indicate posterior mean
plot(seq(nBurnIn,nRuns,1),res[nBurnIn:nRuns],col = '#7570B3', type="l",xlab="iteration",ylab=expression(phi),bty="n",ylim=c(0.6,1))
abline(h=mean(res[nBurnIn:nRuns]),lwd=1,lty="dotted")

#dev.off()

##############################################################################
# End of file
##############################################################################
