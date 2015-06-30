%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of particle Metropolis-Hastings in a stochastic volatility model
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
% x[tt+1] = phi  * x[tt] + sigma   * v[tt]
% y[tt]   = beta * exp( xt[tt]/2 ) * e[tt]
%
% where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)

% Set the number of time steps to simulate
T      = 500;

% Set the initial state
x0     = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = load('omxs30data.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using PMH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The inital guess of the parameter
initPar  = [ 0.99 0.13 0.90 ];

% No. particles in the particle filter ( choose nPart ~ T )
nPart    = 500;

% The length of the burn-in and the no. iterations of PMH ( nBurnIn < nRuns )
nBurnIn  = 2500;
nRuns    = 7500;

% The standard deviation in the random walk proposal
stepSize = [ 0.01 0.05 0.05]

% Run the PMH algorithm
res = sv_pmh(y,initPar,nPart,T,nRuns,stepSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xhat = sv_sm(y,mean(res(nBurnIn:nRuns,1)),mean(res(nBurnIn:nRuns,2)),mean(res(nBurnIn:nRuns,3)),N,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the log-returns and the log-volatility
subplot(4,2,[1 2]);
plot(y,'LineWidth',2,'Color',[27 158 119]/256)
xlabel('time'); ylabel('log-return');
hold on; plot(xhat(2:(T+1)),'LineWidth',2,'Color',[217 95 2]/256)

% Plot the parameter posterior estimate
% Plot the trace of the Markov chain after burn-in
% Solid black line indicate posterior mean
subplot(4,2,3);
hist(res(nBurnIn:nRuns,1), floor(sqrt(nRuns-nBurnIn)) );
xlabel('phi'); ylabel('posterior density estimate');

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[117 112 179]/256,'EdgeColor','w');
hold on; plot([1 1] * mean(res(nBurnIn:nRuns,1)), [0 500], 'k'); hold off;

subplot(4,2,4);
plot(nBurnIn:nRuns,res(nBurnIn:nRuns,1), 'Color', [117 112 179]/256);
hold on; plot([nBurnIn nRuns],[1 1] * mean(res(nBurnIn:nRuns,1)), 'k'); hold off;
xlabel('iteration'); ylabel('trace of phi');

subplot(4,2,5);
hist(res(nBurnIn:nRuns,2), floor(sqrt(nRuns-nBurnIn)) );
xlabel('sigmav'); ylabel('posterior density estimate');

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[231 41 138]/256,'EdgeColor','w');
hold on; plot([1 1] * mean(res(nBurnIn:nRuns,2)), [0 500], 'k'); hold off;

subplot(4,2,6);
plot(nBurnIn:nRuns,res(nBurnIn:nRuns,2), 'Color', [231 41 138]/256);
xlabel('iteration'); ylabel('trace of sigmav');
hold on; plot([nBurnIn nRuns],[1 1] * mean(res(nBurnIn:nRuns,2)), 'k'); hold off;

subplot(4,2,7);
hist(res(nBurnIn:nRuns,3), floor(sqrt(nRuns-nBurnIn)) );
xlabel('beta'); ylabel('posterior density estimate');

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[102 166 30]/256,'EdgeColor','w');
hold on; plot([1 1] * mean(res(nBurnIn:nRuns,3)), [0 500], 'k'); hold off;

subplot(4,2,8);
plot(nBurnIn:nRuns,res(nBurnIn:nRuns,3), 'Color', [102 166 30]/256);
hold on; plot([nBurnIn nRuns],[1 1] * mean(res(nBurnIn:nRuns,3)), 'k'); hold off;
xlabel('iteration'); ylabel('trace of beta');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
