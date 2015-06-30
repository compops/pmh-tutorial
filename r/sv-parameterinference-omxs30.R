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
initPar  = c( 0.99, 0.13, 0.90);

# No. particles in the particle filter ( choose nPart ~ T )
nPart    = 500;

# The length of the burn-in and the no. iterations of PMH ( nBurnIn < nRuns )
nBurnIn  = 2500;
nRuns    = 7500;

# The standard deviation in the random walk proposal
stepSize = c( 0.01, 0.05, 0.05 );

# Run the PMH algorithm
res = pmh_sv(y,initPar,nPart,T,nRuns,stepSize)

##############################################################################
# Plot the results
##############################################################################

# Find the MAP estimate of the parameters
foo1 = density(res[nBurnIn:nRuns,1],kernel="e")
foo2 = density(res[nBurnIn:nRuns,2],kernel="e")
foo3 = density(res[nBurnIn:nRuns,3],kernel="e")
thMAP1 = mean(res[nBurnIn:nRuns,1])
thMAP2 = mean(res[nBurnIn:nRuns,2])
thMAP3 = mean(res[nBurnIn:nRuns,3])

# Estimate the log-volatility
xhat <- sm_sv(y,thMAP1,thMAP2,thMAP3,nPart,T)$xh

# Export plot to file
#cairo_pdf("sv-parameter.pdf", height = 10, width = 8)

# Plot the parameter posterior estimate, solid black line indicate posterior mean
# Plot the trace of the Markov chain after burn-in, solid black line indicate posterior mean
layout(matrix(c(1,1,2,3,4,5,6,7), 4, 2, byrow = TRUE)); 
par(mar=c(4,5,0,4.5))

plot(y,col="#1B9E77",lwd=1.5,type="l",xlab="time",ylab="log-returns",bty="n")
par(new=TRUE)
plot(xhat[-1],lwd=1.5,col="#D95F02",type="l",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(-1,1.5))
axis(4)
mtext("log-volatility",side=4,line=3,cex=0.75)

par(mar=c(4,5,0,0))

# Phi
hist(res[nBurnIn:nRuns,1], breaks=floor(sqrt(nRuns-nBurnIn)), col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25),border=NA,xlab=expression(phi),ylab="posterior estimate",main="",xlim=c(0.85,1),freq=F)
lines(foo1,lwd=2,col="#7570B3"); abline(v=thMAP1,lwd=1,lty="dotted");

plot(seq(nBurnIn,nRuns,1),res[nBurnIn:nRuns,1],col='#7570B3', type="l",xlab="iteration",ylab=expression(phi),bty="n",ylim=c(0.85,1))
abline(h=thMAP1,lwd=1,lty="dotted")

# Sigma_v
hist(res[nBurnIn:nRuns,2], breaks=floor(sqrt(nRuns-nBurnIn)), col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25),border=NA,xlab=expression(sigma[v]),ylab="posterior estimate",main="",xlim=c(0,0.5),freq=F)
lines(foo2,lwd=2,col="#E7298A"); abline(v=thMAP2,lwd=1,lty="dotted");

plot(seq(nBurnIn,nRuns,1),res[nBurnIn:nRuns,2],col='#E7298A', type="l",xlab="iteration",ylab=expression(sigma[v]),bty="n",ylim=c(0,0.5))
abline(h=thMAP2,lwd=1,lty="dotted")

# Beta
hist(res[nBurnIn:nRuns,3], breaks=floor(sqrt(nRuns-nBurnIn)), col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25),border=NA,xlab=expression(beta),ylab="posterior estimate",main="",xlim=c(0.6,1.4),freq=F)
lines(foo3,lwd=2,col="#66A61E"); abline(v=thMAP3,lwd=1,lty="dotted");

plot(seq(nBurnIn,nRuns,1),res[nBurnIn:nRuns,3],col='#66A61E', type="l",xlab="iteration",ylab=expression(beta),bty="n",ylim=c(0.6,1.4))
abline(h=thMAP3,lwd=1,lty="dotted")

#dev.off()

##############################################################################
# End of file
##############################################################################
