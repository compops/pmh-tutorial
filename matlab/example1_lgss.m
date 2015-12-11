%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of fully-adapted particle filtering 
% in a linear Gaussian state space model
%
% Copyright (C) 2015 Johan Dahlin < johan.dahlin (at) liu.se >
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here, we use the following model
%
% x[tt+1] = phi  * x[tt] + sigmav * v[tt]
% y[tt]   = x[tt]        + sigmae * e[tt]
%
% where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

% Set the parameters of the model
phi    = 0.75;
sigmav = 1.00;
sigmae = 0.10;

% Set the number of time steps to simulate
T      = 250;

% Set the initial state
x0     = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,y] = generateData(phi,sigmav,sigmae,T,x0);

% Plot the measurements and latent states
subplot(3,1,1); plot(y(2:(T+1)),'LineWidth',1.5,'Color',[27 158 119]/256); 
xlabel('time'); ylabel('measurement');

subplot(3,1,2); plot(x(2:(T+1)),'LineWidth',1.5,'Color',[217 95 2]/256); 
xlabel('time'); ylabel('latent state');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State estimation using the particle filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using N = 20 particles and plot the estimate of the latent state
N      = 20;
outPF  = sm(y,phi,sigmav,sigmae,N,T,x0);

subplot(3,1,3); plot(1:(T-1),outPF(2:T),'LineWidth',1.5,'Color',[117 112 179]/256);
xlabel('time'); ylabel('state estimate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State estimation using the Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outKF  = kf(y,phi,sigmav,sigmae,T,x0,0.01);
hold on;
plot(1:(T-1),outKF(2:T),'k.','LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
