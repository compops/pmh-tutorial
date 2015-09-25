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
stepSize = diag( c( 0.10, 0.01, 0.05 )^2 );

# Run the PMH algorithm
res = pmh_sv(y,initPar,nPart,T,nIter,stepSize)

##############################################################################
# Plot the results
##############################################################################

resTh = res$thhat[nBurnIn:nIter,];
resXh = res$xhat[nBurnIn:nIter,];

# Find the MAP estimate of the parameters
foo1  = density(resTh[,1],kernel="e",from=-1,to=1)
foo2  = density(resTh[,2],kernel="e",to=1)
foo3  = density(resTh[,3],kernel="e",from=0)
thhat = colMeans(resTh)

# Estimate the log-volatility
xhat = colMeans(resXh)

# Export plot to file
#cairo_pdf("sv-parameter.pdf", height = 10, width = 8)

# Plot the parameter posterior estimate, solid black line indicate posterior mean
# Plot the trace of the Markov chain after burn-in, solid black line indicate posterior mean
layout(matrix(c(1,1,1,2,3,4,5,6,7,8,9,10), 4, 3, byrow = TRUE)); 
par(mar=c(4,5,0,4.5))

plot(y,col="#1B9E77",lwd=1.5,type="l",xlab="time",ylab="log-returns",bty="n")
par(new=TRUE)
plot(xhat[-1],lwd=1.5,col="#D95F02",type="l",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(-1,1.5))
axis(4)
mtext("log-volatility",side=4,line=3,cex=0.75)

par(mar=c(4,5,0,0))

grid=seq(nBurnIn,nBurnIn+1000-1,1);

# Mu
hist(resTh[,1], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25),border=NA,xlab=expression(mu),ylab="posterior estimate",main="",xlim=c(-1,1),freq=F)
lines(seq(-1,1,0.01),dnorm(seq(-1,1,0.01),0,1),col="darkgrey")
lines(foo1,lwd=2,col="#7570B3"); abline(v=thhat[1],lwd=1,lty="dotted");

plot(grid,resTh[1:1000,1],col='#7570B3', type="l",xlab="iteration",ylab=expression(mu),bty="n",ylim=c(-1,1))
polygon(c(grid,rev(grid)),c(resTh[1:1000,1],rep(-1,length(resTh[1:1000,1]))),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))
abline(h=thhat[1],lwd=1,lty="dotted")

foo1=acf( resTh[,1], plot=F, lag.max=100)
plot(foo1$lag,foo1$acf,col = '#7570B3', type="l",xlab="iteration",ylab=expression("ACF of "* mu),bty="n",lwd=2,ylim=c(0,1))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

# Phi
hist(resTh[,2], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25),border=NA,xlab=expression(phi),ylab="posterior estimate",main="",xlim=c(0.88,1.0),freq=F)
lines(seq(0.88,1,0.01),dnorm(seq(0.88,1,0.01),0.95,0.05),col="darkgrey")
lines(foo2,lwd=2,col="#E7298A"); abline(v=thhat[2],lwd=1,lty="dotted");

plot(grid,resTh[1:1000,2],col='#E7298A', type="l",xlab="iteration",ylab=expression(phi),bty="n",ylim=c(0.9,1.1))
polygon(c(grid,rev(grid)),c(resTh[1:1000,2],rep(0,length(resTh[1:1000,2]))),border=NA,col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25))
abline(h=thhat[2],lwd=1,lty="dotted")

foo2=acf( resTh[,2], plot=F, lag.max=100)
plot(foo2$lag,foo2$acf,col = '#E7298A', type="l",xlab="iteration",ylab=expression("ACF of "* phi),bty="n",lwd=2,ylim=c(0,1))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

# Sigma[v]
hist(resTh[,3], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25),border=NA,xlab=expression(sigma[v]),ylab="posterior estimate",main="",xlim=c(0.0,0.4),freq=F)
lines(seq(0,0.4,0.01),dgamma(seq(0,0.4,0.01),2,10),col="darkgrey")
lines(foo3,lwd=2,col="#66A61E"); abline(v=thhat[3],lwd=1,lty="dotted");

plot(grid,resTh[1:1000,3],col='#66A61E', type="l",xlab="iteration",ylab=expression(sigma[v]),bty="n",ylim=c(0.0,0.4))
polygon(c(grid,rev(grid)),c(resTh[1:1000,3],rep(0,length(resTh[1:1000,3]))),border=NA,col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25))
abline(h=thhat[3],lwd=1,lty="dotted")

foo3=acf( resTh[,3], plot=F, lag.max=100)
plot(foo3$lag,foo3$acf,col = '#66A61E', type="l",xlab="iteration",ylab=expression("ACF of "* sigma[v]),bty="n",lwd=2,ylim=c(0,1))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

#dev.off()

# Compute an estimate of the IACT using the first 100 ACF coefficients
iact = 1 + 2 * c( sum(foo1$acf), sum(foo2$acf), sum(foo3$acf) )

##############################################################################
# End of file
##############################################################################
