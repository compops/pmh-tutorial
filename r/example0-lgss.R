##############################################################################
#
# Example of particle filtering 
# in a linear Gaussian state space model
#
# (c) 2015 Johan Dahlin
# johan.dahlin (at) liu.se
#
##############################################################################

# Import helper
source("stateEstimationHelper.R")

# Set the random seed to replicate results in tutorial
set.seed(10)

##############################################################################
# Define the model
##############################################################################

# Here, we use the following model
#
# x[tt+1] = phi   * x[tt] + sigmaV * v[tt]
# y[tt]   = x[tt]         + sigmaE * e[tt]
#
# where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

# Set the parameters of the model
sigmaV  = 0.10
sigmaE  = 1.00

# Set the number of time steps to simulate
T      = 250;

# Set the initial state
x0     = 0;

##############################################################################
# Generate data
##############################################################################

data1 = generateData(0.37,sigmaV,sigmaE,T,x0)
data2 = generateData(0.67,sigmaV,sigmaE,T,x0)
data3 = generateData(0.98,sigmaV,sigmaE,T,x0)

# Export plot to file
cairo_pdf("lgss-data.pdf", height = 8, width = 8)

grid = seq(1,T)

# Plot the latent state and observations
layout(matrix(1:9, 3, 3, byrow = TRUE)); par(mar=c(4,5,0,0))

plot(data1$x,col="#1B9E77",lwd=1,type="l",xlab="time",ylab=expression("latent state " * x[t] ),bty="n",ylim=c(-1,1.5))
polygon(c(grid,rev(grid)),c(data1$x,rep(-1,T)),border=NA,col=rgb(t(col2rgb("#1B9E77"))/256,alpha=0.25))

plot(data1$y,col="#1B9E77",lwd=1,type="l",xlab="time",ylab=expression("observation " * y[t] ),bty="n",ylim=c(-4,4))
polygon(c(grid,rev(grid)),c(data1$y,rep(-4.0,T)),border=NA,col=rgb(t(col2rgb("#1B9E77"))/256,alpha=0.25))

foo = acf(data1$y, plot=F, lag.max=25);
plot(foo$lag,foo$acf,col="#1B9E77",lwd=1.5,type="l",xlab="time",ylab=expression("ACF of "*y[t]),bty="n",ylim=c(-0.15,1),xlim=c(0,25))
abline(h=-1.96/sqrt(T),lty="dotted"); abline(h=1.96/sqrt(T),lty="dotted");

plot(data2$x,col="#D95F02",lwd=1,type="l",xlab="time",ylab=expression("latent state " * x[t] ),bty="n",ylim=c(-1,1.5))
polygon(c(grid,rev(grid)),c(data2$x,rep(-1,T)),border=NA,col=rgb(t(col2rgb("#D95F02"))/256,alpha=0.25))

plot(data2$y,col="#D95F02",lwd=1,type="l",xlab="time",ylab=expression("observation " * y[t] ),bty="n",ylim=c(-4,4))
polygon(c(grid,rev(grid)),c(data2$y,rep(-4.0,T)),border=NA,col=rgb(t(col2rgb("#D95F02"))/256,alpha=0.25))

foo = acf(data21$y, plot=F, lag.max=25);
plot(foo$lag,foo$acf,col="#D95F02",lwd=1.5,type="l",xlab="time",ylab=expression("ACF of "*y[t]),bty="n",ylim=c(-0.15,1),xlim=c(0,25))
abline(h=-1.96/sqrt(T),lty="dotted"); abline(h=1.96/sqrt(T),lty="dotted");

plot(data3$x,col="#7570B3",lwd=1,type="l",xlab="time",ylab=expression("latent state " * x[t] ),bty="n",ylim=c(-1,1.5))
polygon(c(grid,rev(grid)),c(data3$x,rep(-1,T)),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))

plot(data3$y,col="#7570B3",lwd=1,type="l",xlab="time",ylab=expression("observation " * y[t] ),bty="n",ylim=c(-4,4))
polygon(c(grid,rev(grid)),c(data3$y,rep(-4.0,T)),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))

foo = acf(data3$y, plot=F, lag.max=25);
plot(foo$lag,foo$acf,col="#7570B3",lwd=1.5,type="l",xlab="time",ylab=expression("ACF of "*y[t]),bty="n",ylim=c(-0.15,1),xlim=c(0,25))
abline(h=-1.96/sqrt(T),lty="dotted"); abline(h=1.96/sqrt(T),lty="dotted");

dev.off()

###################################################################################
# End of file
###################################################################################

