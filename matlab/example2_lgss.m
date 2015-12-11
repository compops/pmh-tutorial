%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of PMH in a LGSS model
%
% Main script
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
T      = 500;

% Set the initial state
x0     = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,y] = generateData(phi,sigmav,sigmae,T,x0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using PMH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The inital guess of the parameter
initPar  = 0.50;

% No. particles in the particle filter ( choose nPart ~ T )
nPart    = 100;

% The length of the burn-in and the no. iterations of PMH ( nBurnIn < nRuns )
nBurnIn  = 1000;
nRuns    = 5000;

% The standard deviation in the random walk proposal
stepSize = 0.10;

% Run the PMH algorithm
res = pmh(y,initPar,sigmav,sigmae,nPart,T,x0,nRuns,stepSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the parameter posterior estimate
% Solid black line indicate posterior mean
subplot(3,1,1);
hist(res(nBurnIn:nRuns), floor(sqrt(nRuns-nBurnIn)) );
xlabel('phi'); ylabel('posterior density estimate');

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[117 112 179]/256,'EdgeColor','w');

hold on;
plot([1 1] * mean(res(nBurnIn:nRuns)), [0 500], 'LineWidth', 1);
hold off;

% Plot the trace of the Markov chain after burn-in
% Solid black line indicate posterior mean
subplot(3,1,2);
plot(nBurnIn:nRuns,res(nBurnIn:nRuns), 'Color', [117 112 179]/256, 'LineWidth', 1);
xlabel('iteration'); ylabel('phi');

hold on;
plot([nBurnIn,nRuns], [1 1] * mean(res(nBurnIn:nRuns)), 'k', 'LineWidth', 1);
hold off;

% Plot ACF of the Markov chain after burn-in
subplot(3,1,3);
macf = acf( res(nBurnIn:nRuns), 100, 0 );
plot(2:101, macf, 'Color', [117 112 179]/256, 'LineWidth', 2);
xlabel('lag'); ylabel('ACF of phi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
