%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of particle Metropolis-Hastings in a stochastic volatility model
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
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = Quandl.get('NASDAQOMX/OMXS30',...
               'start_date', '2012-01-02',...
               'end_date', '2014-01-02',...
               'type', 'data'); 

y = 100 * diff(log(flipud(d(:, 2))));
T = length(y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using PMH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The inital guess of the parameter
initialTheta  = [0 0.9 0.2];

% No. particles in the particle filter (choose noParticles ~ T)
noParticles = 500;

% The length of the burn-in and the no. iterations of PMH 
% (noBurnInIterations < noIterations)
noBurnInIterations = 2500;
noIterations = 7500;

% The standard deviation in the random walk proposal
stepSize = diag([0.10 0.01 0.05].^2);

% Run the PMH algorithm
[thhat, xhhat] = particleMetropolisHastingsSVmodel(y,...
                                                   initialTheta,...
                                                   noParticles,...
                                                   noIterations,...
                                                   stepSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid = noBurnInIterations:noIterations;
noBins = floor(sqrt(noIterations - noBurnInIterations));
resXh = xhhat(grid, 2:(T + 1));
resTh = thhat(grid, :);

% Plot the log-returns
subplot(5, 3, [1 2 3]);
plot(y, 'LineWidth', 1, 'Color', [27 158 119] / 256)
xlabel('time'); 
ylabel('log-return');

% Plot the log-volatility
subplot(5, 3, [4 5 6]);
plot(mean(resXh, 1), 'LineWidth', 1, 'Color', [217 95 2] / 256)
xlabel('time'); 
ylabel('log-volatility estimate');

% Plot the parameter posterior estimate
% Plot the trace of the Markov chain after burn-in
% Solid black line indicate posterior mean

%--------------------------------------------------------------------------
% Mu
%--------------------------------------------------------------------------
% Histogram of marginal parameter posterior
subplot(5, 3, 7);
hist(resTh(:, 1), noBins);
xlabel('mu'); 
ylabel('posterior density estimate');

h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [117 112 179] / 256, 'EdgeColor', 'w');
hold on; 
plot([1 1] * mean(resTh(:, 1)), [0 500], 'k'); 
hold off;

% Trace plot
subplot(5, 3, 8);
plot(grid, resTh(:, 1), 'Color', [117 112 179] / 256);
hold on; 
plot([grid(1) grid(end)], [1 1] * mean(resTh(:, 1)), 'k'); 
hold off;
xlabel('iteration'); 
ylabel('trace of mu');

% Plot ACF of the Markov chain after burn-in
subplot(5, 3, 9);
[acf, lags] = xcorr(resTh(:, 1) - mean(resTh(:, 1)), 100, 'coeff');
stem(lags(101:200), acf(101:200), 'Color', [117 112 179] / 256, 'LineWidth', 2);
xlabel('lag'); 
ylabel('ACF of mu');


%--------------------------------------------------------------------------
% Phi
%--------------------------------------------------------------------------
% Histogram of marginal parameter posterior
subplot(5, 3, 10);
hist(resTh(:, 2), noBins);
xlabel('phi'); 
ylabel('posterior density estimate');

h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [231 41 138] / 256, 'EdgeColor', 'w');
hold on; 
plot([1 1] * mean(resTh(:, 2)), [0 500], 'k'); 
hold off;

% Trace plot
subplot(5, 3, 11);
plot(grid, resTh(:, 2), 'Color', [231 41 138] / 256);
xlabel('iteration'); 
ylabel('trace of phi');
hold on; 
plot([grid(1) grid(end)],[1 1] * mean(resTh(:, 2)), 'k'); 
hold off;

% Plot ACF of the Markov chain after burn-in
subplot(5, 3, 12);
[acf, lags] = xcorr(resTh(:, 2) - mean(resTh(:, 2)), 100, 'coeff');
stem(lags(101:200), acf(101:200), 'Color', [231 41 138] / 256, 'LineWidth', 2);
xlabel('lag'); 
ylabel('ACF of phi');


%--------------------------------------------------------------------------
% SigmaV
%--------------------------------------------------------------------------
% Histogram of marginal parameter posterior
subplot(5, 3, 13);
hist(resTh(:, 3), noBins);
xlabel('sigmav'); 
ylabel('posterior density estimate');

h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [102 166 30] / 256, 'EdgeColor', 'w');
hold on; 
plot([1 1] * mean(resTh(:, 3)), [0 500], 'k'); 
hold off;

% Trace plot
subplot(5, 3, 14);
plot(grid, resTh(:, 3), 'Color', [102 166 30] / 256);
hold on; 
plot([grid(1) grid(end)],[1 1] * mean(resTh(:, 3)), 'k'); 
hold off;
xlabel('iteration'); 
ylabel('trace of sigmav');

% Plot ACF of the Markov chain after burn-in
subplot(5, 3, 15);
[acf, lags] = xcorr(resTh(:, 3) - mean(resTh(:, 3)), 100, 'coeff');
stem(lags(101:200), acf(101:200), 'Color', [102 166 30] / 256, 'LineWidth', 2);
xlabel('lag'); 
ylabel('ACF of sigmav');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
