%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of Metropolis-Hastings (MH) for the Earthquake data
% with IID model
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
initPar  = [ 0.5 0.5 15 ];

% The length of the burn-in and the no. iterations of MH algorithm 
% ( nBurnIn < nIterations )
nBurnIn     = 10000;
nIterations = 50000;

% The covariance matrix in the random walk proposal
Sigma = [ [ 121.3141    1.4331    3.1442];
          [   1.4331    0.4692    0.0443];
          [   3.1442    0.0443    0.3539] ];
Sigma = 0.8 * Sigma;

% Run the MH algorithm
th = mh( y, initPar, nIterations, Sigma );

% Compute the parameter posterior mean
thhat = mean( th, 1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);

% Plot the data and the fitted model
subplot(4,2,[1 2]);
hist( y, 25 )
xlabel( 'time' ); 
ylabel( 'no. earthquakes' );

hold on; 
grid = 5:0.1:45;
plot( grid, T .* exp( dt( grid, thhat ) ), 'r', 'LineWidth', 2 ); 
hold off;

% Plot the parameter posterior estimate
% Plot the trace of the Markov chain after burn-in
subplot(4,2,3);
hist( th( nBurnIn:nIterations, 1 ), floor( sqrt( nIterations - nBurnIn ) ) );
xlabel('nu'); 
ylabel('posterior density estimate');

subplot(4,2,4);
plot( nBurnIn:nIterations, th( nBurnIn:nIterations, 1 ) );
xlabel('iteration'); ylabel('trace of nu');

subplot(4,2,5);
hist( th( nBurnIn:nIterations, 2 ), floor( sqrt( nIterations - nBurnIn ) ) );
xlabel('mu'); ylabel('posterior density estimate');

subplot(4,2,6);
plot( nBurnIn:nIterations, th( nBurnIn:nIterations, 2 ) );
xlabel('iteration'); ylabel('trace of mu');

subplot(4,2,7);
hist( th( nBurnIn:nIterations, 3 ), floor( sqrt( nIterations - nBurnIn ) ) );
xlabel('sigma'); ylabel('posterior density estimate');

subplot(4,2,8);
plot( nBurnIn:nIterations, th( nBurnIn:nIterations, 3 ) );
xlabel('iteration'); ylabel('trace of sigma');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%