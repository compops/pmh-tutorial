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

# Run the PMH algorithm
res1 = pmh(data$y,initPar,sigmav,sigmae,nPart,T,x0,nIter,stepSize = 0.01)
res2 = pmh(data$y,initPar,sigmav,sigmae,nPart,T,x0,nIter,stepSize = 0.10)
res3 = pmh(data$y,initPar,sigmav,sigmae,nPart,T,x0,nIter,stepSize = 0.50)

##############################################################################
# Plot the results
##############################################################################

# Export plot to file
#cairo_pdf("example2-lgss.pdf", height = 10, width = 8)

layout(matrix(1:9, 3, 3, byrow = TRUE)); par(mar=c(4,5,0,0))

# Plot the parameter posterior estimate
hist(res1[nBurnIn:nIter], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25),border=NA,xlab=expression(phi),ylab="posterior estimate",main="",xlim=c(0.5,0.8),ylim=c(0,12),freq=F)
lines(density(res1[nBurnIn:nIter],kernel="e",from=0.5,to=0.8),lwd=2,col="#7570B3")
lines(seq(0.5,1.0,0.01),dnorm(seq(0.5,1.0,0.01),0,1),lwd=1,col="darkgrey")
abline(v=mean(res1[nBurnIn:nIter]),lwd=1,lty="dotted")

hist(res2[nBurnIn:nIter], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25),border=NA,xlab=expression(phi),ylab="posterior estimate",main="",xlim=c(0.5,0.8),ylim=c(0,12),freq=F)
lines(density(res2[nBurnIn:nIter],kernel="e",from=0.5,to=0.8),lwd=2,col="#E7298A")
lines(seq(0.5,1.0,0.01),dnorm(seq(0.5,1.0,0.01),0,1),lwd=1,col="darkgrey")
abline(v=mean(res2[nBurnIn:nIter]),lwd=1,lty="dotted")

hist(res3[nBurnIn:nIter], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25),border=NA,xlab=expression(phi),ylab="posterior estimate",main="",xlim=c(0.5,0.8),ylim=c(0,12),freq=F)
lines(density(res3[nBurnIn:nIter],kernel="e",from=0.5,to=0.8),lwd=2,col="#66A61E")
lines(seq(0.5,1.0,0.01),dnorm(seq(0.5,1.0,0.01),0,1),lwd=1,col="darkgrey")
abline(v=mean(res3[nBurnIn:nIter]),lwd=1,lty="dotted")

# Plot the trace of the Markov chain during 1000 iterations after the burn-in
plot(seq(nBurnIn,nBurnIn+1000,1),res1[nBurnIn:(nBurnIn+1000)],col = '#7570B3', type="l",xlab="iteration",ylab=expression(phi),bty="n",ylim=c(0.5,0.8))
grid=seq(nBurnIn,nBurnIn+1000,1)
polygon(c(grid,rev(grid)),c(res1[nBurnIn:(nBurnIn+1000)],rep(0,length(res1[nBurnIn:(nBurnIn+1000)]))),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))
abline(h=mean(res1[nBurnIn:nIter]),lwd=1,lty="dotted")

plot(seq(nBurnIn,nBurnIn+1000,1),res2[nBurnIn:(nBurnIn+1000)],col = '#E7298A', type="l",xlab="iteration",ylab=expression(phi),bty="n",ylim=c(0.5,0.8))
grid=seq(nBurnIn,nBurnIn+1000,1)
polygon(c(grid,rev(grid)),c(res2[nBurnIn:(nBurnIn+1000)],rep(0,length(res2[nBurnIn:(nBurnIn+1000)]))),border=NA,col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25))
abline(h=mean(res2[nBurnIn:nIter]),lwd=1,lty="dotted")

plot(seq(nBurnIn,nBurnIn+1000,1),res3[nBurnIn:(nBurnIn+1000)],col = '#66A61E', type="l",xlab="iteration",ylab=expression(phi),bty="n",ylim=c(0.5,0.8))
grid=seq(nBurnIn,nBurnIn+1000,1)
polygon(c(grid,rev(grid)),c(res3[nBurnIn:(nBurnIn+1000)],rep(0,length(res3[nBurnIn:(nBurnIn+1000)]))),border=NA,col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25))
abline(h=mean(res3[nBurnIn:nIter]),lwd=1,lty="dotted")

# Plot the ACF of the Markov chain
foo=acf( res1[nBurnIn:nIter], plot=F, lag.max=60)
plot(foo$lag,foo$acf,col = '#7570B3', type="l",xlab="iteration",ylab="ACF",bty="n",lwd=2,ylim=c(-0.2,1))
polygon(c(foo$lag,rev(foo$lag)),c(foo$acf,rep(0,length(foo$acf))),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

foo=acf( res2[nBurnIn:nIter], plot=F, lag.max=60)
plot(foo$lag,foo$acf,col = '#E7298A', type="l",xlab="iteration",ylab="ACF",bty="n",lwd=2,ylim=c(-0.2,1))
polygon(c(foo$lag,rev(foo$lag)),c(foo$acf,rep(0,length(foo$acf))),border=NA,col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")


foo=acf( res3[nBurnIn:nIter], plot=F, lag.max=60)
plot(foo$lag,foo$acf,col = '#66A61E', type="l",xlab="iteration",ylab="ACF",bty="n",lwd=2,ylim=c(-0.2,1))
polygon(c(foo$lag,rev(foo$lag)),c(foo$acf,rep(0,length(foo$acf))),border=NA,col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

#dev.off()

# Estimate the parameter posterior mean
mean( res1[nBurnIn:nIter] )
mean( res2[nBurnIn:nIter] )
mean( res3[nBurnIn:nIter] )

##############################################################################
# End of file
##############################################################################
