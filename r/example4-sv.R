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
y = as.numeric( read.table("omxs30data.csv")[,1] )

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
stepSize = matrix( c( 0.0726444646,  0.0009686992, -0.0010722330, 
                      0.0009686992,  0.0003498988, -0.0006930118,
                     -0.0010722330, -0.0006930118,  0.0023934189), ncol=3, nrow=3);
stepSize = 0.8 * stepSize;

# Run the PMH algorithm
res = pmh_sv(y,initPar,nPart,T,nIter,stepSize)

##############################################################################
# Plot the results
##############################################################################

# Extract the states after burn-in
resTh = res$thhat[ nBurnIn:nIter, ];
resXh = res$xhat[ nBurnIn:nIter, ];

# Estimate the KDE of the marginal posteriors
kde1  = density( resTh[ , 1], kernel = "e", from = -1, to = 1 );
kde2  = density( resTh[ , 2], kernel = "e", to = 1 );
kde3  = density( resTh[ , 3], kernel = "e", from = 0 );

# Estimate the posterior mean and the corresponding standard deviation
thhat   = colMeans( resTh );
thhatSD = apply( resTh, 2, sd );

# Estimate the log-volatility and the corresponding standad deviation
xhat    = colMeans( resXh );
xhatSD  = apply( resXh, 2, sd );

# Export plot to file
#cairo_pdf("example4-sv.pdf", height = 10, width = 8)

# Plot the parameter posterior estimate, solid black line indicate posterior mean
# Plot the trace of the Markov chain after burn-in, solid black line indicate posterior mean
layout( matrix( c(1,1,1,2,2,2,3,4,5,6,7,8,9,10,11), 5, 3, byrow = TRUE ) ); 
par( mar = c(4,5,0,0) );

# Grid for plotting the data and log-volatility
gridy = seq( 1, length( y ), 1 );

plot( y, col = "#1B9E77", lwd = 1, type = "l", xlab = "time", ylab = "log-returns", bty = "n", ylim = c(-5,5) )
polygon( c( gridy, rev( gridy ) ),c( y, rep( -5, length( y ) ) ), border = NA, col = rgb( t( col2rgb( "#1B9E77" ) ) / 256, alpha = 0.25) )

plot( xhat[-1], col = "#D95F02", lwd = 1.5, type = "l", xlab = "time", ylab = "log-volatility estimate", bty = "n", ylim = c(-2,2) )
polygon( c( gridy, rev( gridy ) ), c( xhat[-1] - 1.96 * xhatSD[-1], rev( xhat[-1] + 1.96 * xhatSD[-1] ) ), border = NA,col = rgb( t( col2rgb( "#D95F02" ) ) /256, alpha = 0.25))

nPlot = 1500;
grid  = seq( nBurnIn, nBurnIn+nPlot-1, 1 );

# Mu
hist( resTh[,1], breaks = floor( sqrt( nIter - nBurnIn ) ), col = rgb( t( col2rgb( "#7570B3" ) ) / 256, alpha = 0.25 ),border=NA,xlab=expression(mu),ylab="posterior estimate",main="",xlim=c(-1,1),freq=F)
lines( seq(-1,1,0.01), dnorm(seq(-1,1,0.01),0,1),col="darkgrey")
lines( kde1, lwd=2,col="#7570B3"); abline(v=thhat[1],lwd=1,lty="dotted");

plot(grid,resTh[1:nPlot,1],col='#7570B3', type="l",xlab="iteration",ylab=expression(mu),bty="n",ylim=c(-1,1),xlim=c(2500,4000))
polygon(c(grid,rev(grid)),c(resTh[1:nPlot,1],rep(-1,length(resTh[1:nPlot,1]))),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))
abline(h=thhat[1],lwd=1,lty="dotted")

foo1=acf( resTh[,1], plot=F, lag.max=100)
plot(foo1$lag,foo1$acf,col = '#7570B3', type="l",xlab="iteration",ylab=expression("ACF of "* mu),bty="n",lwd=2,ylim=c(-0.5,1))
polygon(c(foo1$lag,rev(foo1$lag)),c(foo1$acf,rep(0,length(foo1$acf))),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

# Phi
hist(resTh[,2], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25),border=NA,xlab=expression(phi),ylab="posterior estimate",main="",xlim=c(0.88,1.0),freq=F)
lines(seq(0.88,1,0.01),dnorm(seq(0.88,1,0.01),0.95,0.05),col="darkgrey")
lines(kde2,lwd=2,col="#E7298A"); abline(v=thhat[2],lwd=1,lty="dotted");
plot(grid,resTh[1:nPlot,2],col='#E7298A', type="l",xlab="iteration",ylab=expression(phi),bty="n",ylim=c(0.9,1.1),xlim=c(2500,4000))

polygon(c(grid,rev(grid)),c(resTh[1:nPlot,2],rep(0.9,length(resTh[1:nPlot,2]))),border=NA,col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25))
abline(h=thhat[2],lwd=1,lty="dotted")

foo2=acf( resTh[,2], plot=F, lag.max=100)
plot(foo2$lag,foo2$acf,col = '#E7298A', type="l",xlab="iteration",ylab=expression("ACF of "* phi),bty="n",lwd=2,ylim=c(-0.5,1))
polygon(c(foo2$lag,rev(foo2$lag)),c(foo2$acf,rep(0,length(foo2$acf))),border=NA,col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

# Sigma[v]
hist(resTh[,3], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25),border=NA,xlab=expression(sigma[v]),ylab="posterior estimate",main="",xlim=c(0.0,0.4),freq=F)
lines(seq(0,0.4,0.01),dgamma(seq(0,0.4,0.01),2,10),col="darkgrey")
lines(kde3,lwd=2,col="#66A61E"); abline(v=thhat[3],lwd=1,lty="dotted");

plot(grid,resTh[1:nPlot,3],col='#66A61E', type="l",xlab="iteration",ylab=expression(sigma[v]),bty="n",ylim=c(0.0,0.4),xlim=c(2500,4000))
polygon(c(grid,rev(grid)),c(resTh[1:nPlot,3],rep(0,length(resTh[1:nPlot,3]))),border=NA,col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25))
abline(h=thhat[3],lwd=1,lty="dotted")

foo3=acf( resTh[,3], plot=F, lag.max=100)
plot(foo3$lag,foo3$acf,col = '#66A61E', type="l",xlab="iteration",ylab=expression("ACF of "* sigma[v]),bty="n",lwd=2,ylim=c(-0.5,1))
polygon(c(foo3$lag,rev(foo3$lag)),c(foo3$acf,rep(0,length(foo3$acf))),border=NA,col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25))
abline(h=1.96/sqrt(nIter-nBurnIn),lty="dotted")
abline(h=-1.96/sqrt(nIter-nBurnIn),lty="dotted")

#dev.off()

# Compute an estimate of the IACT using the first 100 ACF coefficients
iact = 1 + 2 * c( sum(foo1$acf), sum(foo2$acf), sum(foo3$acf) )

# Estimate the covariance of the posterior to tune the proposal
estCov = var( resTh )

##############################################################################
# End of file
##############################################################################
