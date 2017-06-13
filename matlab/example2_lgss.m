%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of particle Metropolis-Hastings in a LGSS model.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using PMH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The inital guess of the parameter
initialPhi = 0.50;

% No. particles in the particle filter (choose noParticles ~ T)
noParticles = 100;

% The length of the burn-in and the no. iterations of PMH 
% (noBurnInIterations < noIterations)
noBurnInIterations = 1000;
noIterations = 5000;

% The standard deviation in the random walk proposal
stepSize = 0.10;

% Run the PMH algorithm
res = particleMetropolisHastings(y, initialPhi, [sigmav sigmae],...
                                 noParticles, initialState,...
                                 noIterations, stepSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noBins = floor(sqrt(noIterations - noBurnInIterations));
grid = noBurnInIterations:noIterations;
resPhi = res(noBurnInIterations:noIterations);

% Plot the parameter posterior estimate
% Solid black line indicate posterior mean
subplot(3, 1, 1);
hist(resPhi, noBins);
xlabel('phi'); 
ylabel('posterior density estimate');

h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [117 112 179] / 256, 'EdgeColor', 'w');

hold on;
plot([1 1] * mean(resPhi), [0 200], 'LineWidth', 3);
hold off;

% Plot the trace of the Markov chain after burn-in
% Solid black line indicate posterior mean
subplot(3, 1, 2);
plot(grid, resPhi, 'Color', [117 112 179] / 256, 'LineWidth', 1);
xlabel('iteration'); 
ylabel('phi');

hold on;
plot([grid(1) grid(end)], [1 1] * mean(resPhi), 'k', 'LineWidth', 3);
hold off;

% Plot ACF of the Markov chain after burn-in
subplot(3, 1, 3);
[acf, lags] = xcorr(resPhi - mean(resPhi), 100, 'coeff');
stem(lags(101:200), acf(101:200), 'Color', [117 112 179] / 256, 'LineWidth', 2);
xlabel('lag'); 
ylabel('ACF of phi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
