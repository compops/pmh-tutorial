%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of fully-adapted particle filtering 
% in a linear Gaussian state space model
%
% (c) 2015 Johan Dahlin
% johan.dahlin (at) liu.se
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
