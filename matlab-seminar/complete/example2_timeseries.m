%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Particle Metropolis-Hastings (PMH) 
% for the Earthquake data with state space model.
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

% No. particles in the particle filter ( choose nPart ~ T )
nParticles    = 100;

% The length of the burn-in and the no. iterations of the PMH algorithm
% ( nBurnIn < nIter )
nBurnIn     = 500;
nIterations = 2000;

% The covariance matrix in the random walk proposal
Sigma = diag( [0.07 0.03 2].^2 ); 
Sigma = 0.8 * Sigma;

% Run the PMH algorithm
th = pmh( y, initPar, nParticles, nIterations, Sigma );

% Compute the parameter posterior mean
thhat = mean( th, 1 );

% Run a particle filter to estimate the latent state
xhhat = pf( y, thhat, nParticles );

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
hist( th( nBurnIn:nIterations, 1 ), floor( sqrt( nIterations - nBurnIn ) ) );
xlabel('phi'); 
ylabel('posterior density estimate');

subplot(4,2,4);
plot( nBurnIn:nIterations, th( nBurnIn:nIterations, 1 ) );
xlabel('iteration'); ylabel('trace of phi');

subplot(4,2,5);
hist( th( nBurnIn:nIterations, 2 ), floor( sqrt( nIterations - nBurnIn ) ) );
xlabel('sigmav'); ylabel('posterior density estimate');

subplot(4,2,6);
plot( nBurnIn:nIterations, th( nBurnIn:nIterations, 2 ) );
xlabel('iteration'); ylabel('trace of sigmav');

subplot(4,2,7);
hist( th( nBurnIn:nIterations, 3 ), floor( sqrt( nIterations - nBurnIn ) ) );
xlabel('beta'); ylabel('posterior density estimate');

subplot(4,2,8);
plot( nBurnIn:nIterations, th( nBurnIn:nIterations, 3 ) );
xlabel('iteration'); ylabel('trace of beta');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%