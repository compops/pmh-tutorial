##############################################################################
#
# Example of fully-adapted particle filtering 
# in a linear Gaussian state space model
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
cairo_pdf("example1-lgss.pdf", height = 10, width = 8)

# Plot the latent state and observations
layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE));
par(mar=c(4,5,0,0))
grid = seq(0,T);

plot(grid,y,col="#1B9E77",lwd=1.5,type="l",xlab="time",ylab="observation",bty="n",ylim=c(-6,6))
polygon(c(grid,rev(grid)),c(y,rep(-7,T+1)),border=NA,col=rgb(t(col2rgb("#1B9E77"))/256,alpha=0.25))

###################################################################################
# State estimation using the particle filter and Kalman filter
###################################################################################

# Using N = 20 particles and plot the estimate of the latent state
N      <- 20;
outPF  <- sm(y,phi,sigmav,sigmae,N,T,x0)
outKF  <- kf(y,phi,sigmav,sigmae,x0,0.01);
diff   <- outPF$xh-outKF$xh[-(T+1)];

grid = seq(0,T-1);
plot(grid,diff,col="#7570B3",lwd=1.5,type="l",xlab="time",ylab="difference in state estimate",bty="n",ylim=c(-0.1,0.1))
polygon(c(grid,rev(grid)),c(diff,rep(-0.10,T)),border=NA,col=rgb(t(col2rgb("#7570B3"))/256,alpha=0.25))

# Compute bias and MSE
logBiasMSE = matrix( 0, nrow=7, ncol=2);
gridN      = c( 10, 20, 50, 100, 200, 500, 1000);

for ( ii in 1:length(gridN) ) {
  logBiasMSE[ ii, 1 ] = log( mean( abs(sm(y,phi,sigmav,sigmae,gridN[ii],T,x0)$xh-outKF$xh[-(T+1)]) ) )  
  logBiasMSE[ ii, 2 ] = log( mean( (sm(y,phi,sigmav,sigmae,gridN[ii],T,x0)$xh-outKF$xh[-(T+1)])^2 ) )
}

# Plot the bias and MSE for comparison
plot(   gridN,logBiasMSE[ , 1 ],col="#E7298A",lwd=1.5,type="l",xlab="no. particles (N)",ylab="log-bias",bty="n",ylim=c(-7,-3))
points( gridN,logBiasMSE[ , 1 ],col="#E7298A",pch=19)
polygon(c(gridN,rev(gridN)),c(logBiasMSE[ , 1 ],rep(-7,length(gridN))),border=NA,col=rgb(t(col2rgb("#E7298A"))/256,alpha=0.25))

plot(   gridN,logBiasMSE[ , 2 ],col="#66A61E",lwd=1.5,type="l",xlab="no. particles (N)",ylab="log-MSE",bty="n",ylim=c(-12,-6))
points( gridN,logBiasMSE[ , 2 ],col="#66A61E",pch=19)
polygon(c(gridN,rev(gridN)),c(logBiasMSE[ , 2 ],rep(-12,length(gridN))),border=NA,col=rgb(t(col2rgb("#66A61E"))/256,alpha=0.25))

dev.off()

###################################################################################
# End of file
###################################################################################
