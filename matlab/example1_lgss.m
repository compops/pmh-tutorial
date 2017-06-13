%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of state estimation in a LGSS model 
% using particle filters and Kalman filters
%
% Copyright (C) 2017 Johan Dahlin < liu (at) johandahlin.com.nospam >
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

% Set random seed
rng(0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here, we use the following model
%
% x[t + 1] = phi * x[t] + sigmav * v[t]
% y[t] = x[t] + sigmae * e[t]
%
% where v[t] ~ N(0, 1) and e[t] ~ N(0, 1)

% Set the parameters of the model
phi = 0.75;
sigmav = 1.00;
sigmae = 0.10;
theta = [phi sigmav sigmae];

% Set the number of time steps to simulate
T = 250;

% Set the initial state
initialState = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x, y] = generateData(theta, T, initialState);

% Plot the measurements and latent states
subplot(3,1,1); 
plot(y(2:(T + 1)), 'LineWidth', 1.5, 'Color', [27 158 119] / 256); 
xlabel('time'); 
ylabel('measurement');

subplot(3,1,2); 
plot(x(2:(T + 1)), 'LineWidth', 1.5, 'Color', [217 95 2] / 256); 
xlabel('time'); 
ylabel('latent state');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using N = 20 particles and plot the estimate of the latent state
outPF  = particleFilter(y, theta, 20, initialState);

% Kalman filter
outKF  = kalmanFilter(y, theta, initialState, 0.01);

% Plot the difference
subplot(3,1,3); 
difference = outPF(2:T) - outKF(2:T);
plot(1:(T - 1), difference, 'LineWidth', 1.5, 'Color', [117 112 179] / 256);
xlabel('time'); 
ylabel('difference in state estimate');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
