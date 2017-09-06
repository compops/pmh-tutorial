%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of particle Metropolis-Hastings in a stochastic volatility model
%
% Johan Dahlin <liu (at) johandahlin.com.nospam>
% Documentation at https://github.com/compops/pmh-tutorial
% Published under GNU General Public License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random seed
rng(0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = Quandl.get('NASDAQOMX/OMXS30', 'start_date', '2012-01-02', 'end_date', '2014-01-02', 'type', 'data'); 
logReturns = 100 * diff(log(flipud(data(:, 2))));
noObservations = length(logReturns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialTheta  = [0 0.9 0.2];
noParticles = 500;              % Use noParticles ~ noObservations
noBurnInIterations = 2500;
noIterations = 7500;
stepSize = diag([0.10 0.01 0.05].^2);

[parameterTrace, logVolatilityEstimate] = particleMetropolisHastingsSVmodel(logReturns, initialTheta, noParticles, noIterations, stepSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid = noBurnInIterations:noIterations;
noBins = floor(sqrt(noIterations - noBurnInIterations));
logVolatilityEstimate = logVolatilityEstimate(grid, 2:(noObservations + 1));
parameterTrace = parameterTrace(grid, :);

% Plot the log-returns
subplot(5, 3, [1 2 3]);
plot(logReturns, 'LineWidth', 1, 'Color', [27 158 119] / 256)
xlabel('time'); 
ylabel('log-return');

% Plot the log-volatility
subplot(5, 3, [4 5 6]);
plot(mean(logVolatilityEstimate, 1), 'LineWidth', 1, 'Color', [217 95 2] / 256)
xlabel('time'); 
ylabel('log-volatility estimate');

% Histogram of marginal parameter posterior of mu
subplot(5, 3, 7);
hist(parameterTrace(:, 1), noBins);
xlabel('mu'); 
ylabel('posterior density estimate');

h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [117 112 179] / 256, 'EdgeColor', 'w');
hold on; 
plot([1 1] * mean(parameterTrace(:, 1)), [0 500], 'k'); 
hold off;

% Trace plot for mu
subplot(5, 3, 8);
plot(grid, parameterTrace(:, 1), 'Color', [117 112 179] / 256);
hold on; 
plot([grid(1) grid(end)], [1 1] * mean(parameterTrace(:, 1)), 'k'); 
hold off;
xlabel('iteration'); 
ylabel('trace of mu');

% Plot ACF of the Markov chain for mu after burn-in
subplot(5, 3, 9);
[acf, lags] = xcorr(parameterTrace(:, 1) - mean(parameterTrace(:, 1)), 100, 'coeff');
stem(lags(101:200), acf(101:200), 'Color', [117 112 179] / 256, 'LineWidth', 2);
xlabel('lag'); 
ylabel('ACF of mu');

% Histogram of marginal parameter posterior of phi
subplot(5, 3, 10);
hist(parameterTrace(:, 2), noBins);
xlabel('phi'); 
ylabel('posterior density estimate');

h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [231 41 138] / 256, 'EdgeColor', 'w');
hold on; 
plot([1 1] * mean(parameterTrace(:, 2)), [0 500], 'k'); 
hold off;

% Trace plot for phi
subplot(5, 3, 11);
plot(grid, parameterTrace(:, 2), 'Color', [231 41 138] / 256);
xlabel('iteration'); 
ylabel('trace of phi');
hold on; 
plot([grid(1) grid(end)],[1 1] * mean(parameterTrace(:, 2)), 'k'); 
hold off;

% Plot ACF of the Markov chain for phi after burn-in
subplot(5, 3, 12);
[acf, lags] = xcorr(parameterTrace(:, 2) - mean(parameterTrace(:, 2)), 100, 'coeff');
stem(lags(101:200), acf(101:200), 'Color', [231 41 138] / 256, 'LineWidth', 2);
xlabel('lag'); 
ylabel('ACF of phi');

% Histogram of marginal parameter posterior of sigma_v
subplot(5, 3, 13);
hist(parameterTrace(:, 3), noBins);
xlabel('sigmav'); 
ylabel('posterior density estimate');

h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [102 166 30] / 256, 'EdgeColor', 'w');
hold on; 
plot([1 1] * mean(parameterTrace(:, 3)), [0 500], 'k'); 
hold off;

% Trace plot of sigma_v
subplot(5, 3, 14);
plot(grid, parameterTrace(:, 3), 'Color', [102 166 30] / 256);
hold on; 
plot([grid(1) grid(end)],[1 1] * mean(parameterTrace(:, 3)), 'k'); 
hold off;
xlabel('iteration'); 
ylabel('trace of sigmav');

% Plot ACF of the Markov chain of sigma_v after burn-in
subplot(5, 3, 15);
[acf, lags] = xcorr(parameterTrace(:, 3) - mean(parameterTrace(:, 3)), 100, 'coeff');
stem(lags(101:200), acf(101:200), 'Color', [102 166 30] / 256, 'LineWidth', 2);
xlabel('lag'); 
ylabel('ACF of sigmav');