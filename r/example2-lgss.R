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
sigmae  = 0.10

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

# No. particles in the particle filter
nPart    = 100;

# The length of the burn-in and the no. iterations of PMH ( nBurnIn < nIter )
nBurnIn  = 1000;
nIter    = 5000;

# The standard deviation in the random walk proposal
stepSize = 0.10;

# Run the PMH algorithm
res = pmh(data$y,initPar,sigmav,sigmae,nPart,T,x0,nIter,stepSize)

##############################################################################
# Plot the results
##############################################################################

# Export plot to file
#cairo_pdf("lgss-parameter.pdf", height = 9, width = 8)

# Plot the parameter posterior estimate
# Solid black line indicate posterior mean
layout(matrix(1:3, 3, 1, byrow = TRUE)); par(mar=c(4,5,0,0))
hist(res[nBurnIn:nIter], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25),border=NA,xlab=expression(phi),ylab="posterior estimate",main="",xlim=c(0.5,0.8),ylim=c(0,12),freq=F)
lines(density(res[nBurnIn:nIter],kernel="e",from=0.5,to=0.8),lwd=2,col="#7570B3")
lines(seq(0.5,1.0,0.01),dnorm(seq(0.5,1.0,0.01),0,1),lwd=1,col="darkgrey")
abline(v=mean(res[nBurnIn:nIter]),lwd=1,lty="dotted")

# Plot the trace of the Markov chain after burn-in
# Solid black line indicate posterior mean
plot(seq(nBurnIn,nIter,1),res[nBurnIn:nIter],col = '#7570B3', type="l",xlab="iteration",ylab=expression(phi),bty="n",ylim=c(0.5,0.8))
grid=seq(nBurnIn,nIter,1)
polygon(c(grid,rev(grid)),c(res[nBurnIn:nIter],rep(0,length(res[nBurnIn:nIter]))),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))
abline(h=mean(res[nBurnIn:nIter]),lwd=1,lty="dotted")

# Plot the ACF of the Markov chain
foo=acf( res[nBurnIn:nIter], plot=F)
plot(foo$lag,foo$acf,col = '#7570B3', type="l",xlab="iteration",ylab="ACF",bty="n",lwd=2)
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

#dev.off()

mean( res[nBurnIn:nIter] )

##############################################################################
# End of file
##############################################################################
