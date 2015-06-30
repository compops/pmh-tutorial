%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of PMH in a LGSS model
%
% Main script
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
sigmae = 1.00;

% Set the number of time steps to simulate
T      = 250;

% Set the initial state
x0     = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,y] = lgss_generateData(phi,sigmav,sigmae,T,x0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using PMH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The inital guess of the parameter
initPar  = 0.50;

% No. particles in the particle filter ( choose nPart ~ T )
nPart    = 500;

% The length of the burn-in and the no. iterations of PMH ( nBurnIn < nRuns )
nBurnIn  = 1000;
nRuns    = 5000;

% The standard deviation in the random walk proposal
stepSize = 0.10;

% Run the PMH algorithm
res = lgss_pmh(y,initPar,sigmae,sigmav,nPart,T,x0,nRuns,stepSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the parameter posterior estimate
% Solid black line indicate posterior mean
subplot(2,1,1);
hist(res(nBurnIn:nRuns), floor(sqrt(nRuns-nBurnIn)) );
xlabel('phi'); ylabel('posterior density estimate');

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[117 112 179]/256,'EdgeColor','w');

hold on;
plot([1 1] * mean(res(nBurnIn:nRuns)), [0 500], 'LineWidth', 1);
hold off;

% Plot the trace of the Markov chain after burn-in
% Solid black line indicate posterior mean
subplot(2,1,2);
plot(nBurnIn:nRuns,res(nBurnIn:nRuns), 'Color', [231 41 138]/256, 'LineWidth', 1);
xlabel('iteration'); ylabel('phi');

hold on;
plot([nBurnIn,nRuns], [1 1] * mean(res(nBurnIn:nRuns)), 'k', 'LineWidth', 1);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
