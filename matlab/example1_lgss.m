%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State estimation in a LGSS model using particle and Kalman filters
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

subplot(3,1,1); 
plot(observations(2:(noObservations + 1)), 'LineWidth', 1.5, 'Color', [27 158 119] / 256); 
xlabel('time'); 
ylabel('measurement');

subplot(3,1,2); 
plot(states(2:(noObservations + 1)), 'LineWidth', 1.5, 'Color', [217 95 2] / 256); 
xlabel('time'); 
ylabel('latent state');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Particle filter with N = 20 particles
stateEstPF = particleFilter(observations, parameters, 20, initialState);

% Kalman filter
stateEstKF = kalmanFilter(observations, parameters, initialState, 0.01);

subplot(3,1,3); 
difference = stateEstPF(2:noObservations) - stateEstKF(2:noObservations);
plot(1:(noObservations - 1), difference, 'LineWidth', 1.5, 'Color', [117 112 179] / 256);
xlabel('time'); 
ylabel('difference in state estimate');