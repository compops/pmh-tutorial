%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of particle Metropolis-Hastings in a LGSS model.
%
% Johan Dahlin <liu (at) johandahlin.com.nospam>
% Documentation at https://github.com/compops/pmh-tutorial
% Published under GNU General Public License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random seed
rng(0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the model and generate data
% x[t + 1] = phi * x[t] + sigmav * v[t],    v[t] ~ N(0, 1)
% y[t] = x[t] + sigmae * e[t],              e[t] ~ N(0, 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = 0.75;
sigmav = 1.00;
sigmae = 0.10;
parameters = [phi sigmav sigmae];
noObservations = 250;
initialState = 0;

[states, observations] = generateData(parameters, noObservations, initialState);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialPhi = 0.50;
noParticles = 100;          % Use noParticles ~ noObservations
noBurnInIterations = 1000;
noIterations = 5000;
stepSize = 0.10;

phiTrace = particleMetropolisHastings(observations, initialPhi, [sigmav sigmae], noParticles, initialState, noIterations, stepSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noBins = floor(sqrt(noIterations - noBurnInIterations));
grid = noBurnInIterations:noIterations;
phiTrace = phiTrace(noBurnInIterations:noIterations);

% Plot the parameter posterior estimate (solid black line = posterior mean)
subplot(3, 1, 1);
hist(phiTrace, noBins);
xlabel('phi'); 
ylabel('posterior density estimate');

h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [117 112 179] / 256, 'EdgeColor', 'w');

hold on;
plot([1 1] * mean(phiTrace), [0 200], 'LineWidth', 3);
hold off;

% Plot the trace of the Markov chain after burn-in  (solid black line = posterior mean)
subplot(3, 1, 2);
plot(grid, phiTrace, 'Color', [117 112 179] / 256, 'LineWidth', 1);
xlabel('iteration'); 
ylabel('phi');

hold on;
plot([grid(1) grid(end)], [1 1] * mean(phiTrace), 'k', 'LineWidth', 3);
hold off;

% Plot ACF of the Markov chain after burn-in
subplot(3, 1, 3);
[acf, lags] = xcorr(phiTrace - mean(phiTrace), 100, 'coeff');
stem(lags(101:200), acf(101:200), 'Color', [117 112 179] / 256, 'LineWidth', 2);
xlabel('lag'); 
ylabel('ACF of phi');