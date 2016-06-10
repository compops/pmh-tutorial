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
initPar  = mean(y);

% The length of the burn-in and the no. iterations of MH algorithm 
% ( nBurnIn < nIterations )
nBurnIn     = 10000;
nIterations = 50000;

% The covariance matrix in the random walk proposal
Sigma = 0.17;

% Run the MH algorithm
th = mh( y, initPar, nIterations, Sigma );

% Compute the parameter posterior mean
thhat = mean( th, 1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);

% Plot the data and the fitted model
subplot(2,2,[1 2]);
hist( y, 25 )
xlabel( 'time' ); 
ylabel( 'no. earthquakes' );

hold on; 
grid = 5:1:45;
plot( grid, T .* exp( dpoisson( grid, thhat ) ), 'r', 'LineWidth', 2 ); 
hold off;

% Plot the parameter posterior estimate
% Plot the trace of the Markov chain after burn-in
subplot(2,2,3);
hist( th( nBurnIn:nIterations, 1 ), floor( sqrt( nIterations - nBurnIn ) ) );
xlabel('theta'); 
ylabel('posterior density estimate');

subplot(2,2,4);
plot( nBurnIn:nIterations, th( nBurnIn:nIterations, 1 ) );
xlabel('iteration'); ylabel('trace of theta');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
