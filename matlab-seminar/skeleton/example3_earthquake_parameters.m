%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of Particle Metropolis-Hastings (PMH) for the Earthquake model
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
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('earthquake_data.mat');
T = 114;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using PMH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The inital guess of the parameter

%///////////////              ADD CODE             ///////////////
initPar  = ...

% No. particles in the particle filter

%///////////////              ADD CODE             ///////////////
nPart    = ...

% The length of the burn-in and the no. iterations of PMH 
% ( nBurnIn < nIter )
nBurnIn  = 250;
nIter    = 2000;

% The covariance matrix in the random walk proposal

%///////////////              ADD CODE             ///////////////
stepSize = ...

% Run the PMH algorithm
[th, xh] = pmh_earthquake( y, initPar, nPart, T, nIter, stepSize );

% Compute posterior means
%///////////////              ADD CODE             ///////////////

thhat = 
xhhat = 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);

% Plot the expected versus observed no. earthquakes
subplot(4,2,[1 2]);
plot( y )
xlabel( 'time' ); 
ylabel( 'no. earthquakes' );

hold on; 
plot( thhat(3) * exp(xhhat), 'r' ); 
hold off;

% Plot the parameter posterior estimate
% Plot the trace of the Markov chain after burn-in
subplot(4,2,3);
hist( th( nBurnIn:nIter, 1 ), floor( sqrt( nIter - nBurnIn ) ) );
xlabel('phi'); 
ylabel('posterior density estimate');

subplot(4,2,4);
plot( nBurnIn:nIter, th( nBurnIn:nIter, 1 ) );
xlabel('iteration'); ylabel('trace of phi');

subplot(4,2,5);
hist( th( nBurnIn:nIter, 2 ), floor( sqrt( nIter - nBurnIn ) ) );
xlabel('sigmav'); ylabel('posterior density estimate');

subplot(4,2,6);
plot( nBurnIn:nIter, th( nBurnIn:nIter, 2 ) );
xlabel('iteration'); ylabel('trace of sigmav');

subplot(4,2,7);
hist( th( nBurnIn:nIter, 3 ), floor( sqrt( nIter - nBurnIn ) ) );
xlabel('beta'); ylabel('posterior density estimate');

subplot(4,2,8);
plot( nBurnIn:nIter, th( nBurnIn:nIter, 3 ) );
xlabel('iteration'); ylabel('trace of beta');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%