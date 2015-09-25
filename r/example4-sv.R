##############################################################################
#
# Example of particle Metropolis-Hastings 
# in a stochastic volatility model
#
#
# (c) 2015 Johan Dahlin
# johan.dahlin (at) liu.se
#
##############################################################################

# Import helpers
library("mvtnorm")

source("stateEstimationHelper.R")
source("parameterEstimationHelper.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

##############################################################################
# Define the model
##############################################################################

# Here, we use the following model
#
# x[tt+1] = phi  * x[tt] + sigma   * v[tt]
# y[tt]   = beta * exp( xt[tt]/2 ) * e[tt]
#
# where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

# Set the number of time steps to simulate
T      = 500;

##############################################################################
# Load data
##############################################################################
y = as.numeric( read.table("omxs30data.csv")$V1 )

##############################################################################
# Parameter estimation using PMH
##############################################################################

# The inital guess of the parameter
initPar  = c( 0, 0.9, 0.2);

# No. particles in the particle filter ( choose nPart ~ T )
nPart    = 500;

# The length of the burn-in and the no. iterations of PMH ( nBurnIn < nIter )
nBurnIn  = 2500;
nIter    = 7500;

# The standard deviation in the random walk proposal
stepSize = matrix( c( 0.0540146363,  0.0002821836, -0.0005373284, 
                      0.0002821836,  0.0003919984, -0.0008584207,
                     -0.0005373284, -0.0008584207,  0.0028724014), ncol=3, nrow=3);
stepSize = 0.8 * stepSize;

# Run the PMH algorithm
res = pmh_sv(y,initPar,nPart,T,nIter,stepSize)

##############################################################################
# Plot the results
##############################################################################

resTh = res$thhat[nBurnIn:nIter,];

# Find the MAP estimate of the parameters
foo1  = density(resTh[,1],kernel="e",from=-1,to=1)
foo2  = density(resTh[,2],kernel="e",to=1)
foo3  = density(resTh[,3],kernel="e",from=0)
thhat = colMeans(resTh)

# Export plot to file
# cairo_pdf("sv-parameter-2.pdf", height = 3, width = 8)

# Plot the ACF for each dimension of the Markov chain
layout(matrix(1:3, 1, 3, byrow = TRUE)); 
par(mar=c(4,5,0,0))
grid=seq(nBurnIn,nBurnIn+1000-1,1);

# Mu
foo1=acf( resTh[,1], plot=F, lag.max=100)
plot(foo1$lag,foo1$acf,col = '#7570B3', type="l",xlab="iteration",ylab=expression("ACF of "* mu),bty="n",lwd=2,ylim=c(-0.2,1))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

# Phi
foo2=acf( resTh[,2], plot=F, lag.max=100)
plot(foo2$lag,foo2$acf,col = '#E7298A', type="l",xlab="iteration",ylab=expression("ACF of "* phi),bty="n",lwd=2,ylim=c(-0.2,1))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

# Sigma[v]
foo3=acf( resTh[,3], plot=F, lag.max=100)
plot(foo3$lag,foo3$acf,col = '#66A61E', type="l",xlab="iteration",ylab=expression("ACF of "* sigma[v]),bty="n",lwd=2,ylim=c(-0.2,1))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

# dev.off()

# Compute an estimate of the IACT using the first 100 ACF coefficients
iact = 1 + 2 * c( sum(foo1$acf), sum(foo2$acf), sum(foo3$acf) )

##############################################################################
# End of file
##############################################################################
