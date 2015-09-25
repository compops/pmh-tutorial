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
x    = data$x
y    = data$y

# Export plot to file
#cairo_pdf("lgss-state.pdf", height = 8, width = 8)

# Plot the latent state and observations
layout(matrix(1:3, 3, 1, byrow = TRUE));
par(mar=c(4,5,0,0))
grid = seq(0,T);

plot(y,col="#1B9E77",lwd=1.5,type="l",xlab="time",ylab="observation",bty="n",ylim=c(-6,6))
polygon(c(grid,rev(grid)),c(y,rep(-7,T)),border=NA,col=rgb(t(col2rgb("#1B9E77"))/256,alpha=0.25))

plot(x,col="#D95F02",lwd=1.5,type="l",xlab="time",ylab="latent state",bty="n",ylim=c(-6,6))
polygon(c(grid,rev(grid)),c(x,rep(-7,T)),border=NA,col=rgb(t(col2rgb("#D95F02"))/256,alpha=0.25))

###################################################################################
# State estimation using the particle filter
###################################################################################

grid = seq(0,T-1);

# Using N = 20 particles and plot the estimate of the latent state
N      <- 20;
outPF  <- sm(y,phi,sigmav,sigmae,N,T,x0)
plot(grid,outPF$xh,col="#7570B3",lwd=1.5,type="l",xlab="time",ylab="state estimate",bty="n",ylim=c(-6,6))
polygon(c(grid,rev(grid)),c(outPF$xh,rep(-7,T)),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))

###################################################################################
# State estimation using the Kalman filter
###################################################################################

outKF  <- kf(y,phi,sigmav,sigmae,x0,0.01);
lines( grid, outKF$xh[-(T+1)], col="black", lty="dashed",lwd=1.5)

#dev.off()

# Compute bias
log( mean( abs(sm(y,phi,sigmav,sigmae,10,T,x0)$xh-outKF$xh[-(T+1)]) ) )
log( mean( abs(sm(y,phi,sigmav,sigmae,20,T,x0)$xh-outKF$xh[-(T+1)]) ) )
log( mean( abs(sm(y,phi,sigmav,sigmae,50,T,x0)$xh-outKF$xh[-(T+1)]) ) )
log( mean( abs(sm(y,phi,sigmav,sigmae,100,T,x0)$xh-outKF$xh[-(T+1)]) ) )
log( mean( abs(sm(y,phi,sigmav,sigmae,200,T,x0)$xh-outKF$xh[-(T+1)]) ) )

# Compute MSE
log( mean( (sm(y,phi,sigmav,sigmae,10,T,x0)$xh-outKF$xh[-(T+1)])^2 ) )
log( mean( (sm(y,phi,sigmav,sigmae,20,T,x0)$xh-outKF$xh[-(T+1)])^2 ) )
log( mean( (sm(y,phi,sigmav,sigmae,50,T,x0)$xh-outKF$xh[-(T+1)])^2 ) )
log( mean( (sm(y,phi,sigmav,sigmae,100,T,x0)$xh-outKF$xh[-(T+1)])^2 ) )
log( mean( (sm(y,phi,sigmav,sigmae,200,T,x0)$xh-outKF$xh[-(T+1)])^2 ) )

###################################################################################
# End of file
###################################################################################