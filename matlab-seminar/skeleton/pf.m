%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Bootstrap particle filter 
% for the Earthquake data with state space model.
%
% Inputs:
% observations:             observations from the system for t=1,...,T.
%  
% theta:                    parameter of the model
%
% nParticles:               number of particles (N)
%
% x0:                       the initial state.
%
% Outputs:
% state_estimate:           estimate of the filtered state for each t.
%
% loglikelihood_estimate:   estimate of the log-likelihood at T-1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[ state_estimate, loglikelihood_estimate ] = pf( observations, theta, nParticles )

  %===========================================================
  % Initialise variables
  %===========================================================
  nObservations     = length( observations );
  
  state_estimate    = zeros( nObservations+1, 1 );
  particles         = zeros( nParticles, nObservations+1 );
  weights           = ones(  nParticles, nObservations+1 ) / nParticles;
  
  loglikelihood_estimate = 0;
  particles(:,1)         = 0;
  state_estimate(1)      = 0;
    
  %===========================================================
  % Run main loop
  %===========================================================
  for tt = 2:(nObservations+1)

    %=========================================================
    % Resample ( multinomial )
    %=========================================================
    idx = randsample( nParticles, nParticles, true, weights(:,tt-1) ); 
    
    %=========================================================
    % Propagate
    %=========================================================
    particles(:,tt) = theta(1) .* particles(idx,tt-1) + theta(2) .* normrnd( 0, 1, nParticles, 1 );
    
    %=========================================================
    % Compute weights
    %=========================================================
    weights(:,tt) = dpoisson( observations(tt-1), theta(3) * exp( particles(:,tt) ) );
    
    % Rescale log-weights and recover weight
    wmax    = max( weights(:,tt) );
    weights(:,tt) = exp( weights(:,tt) - wmax );
    
    % Estimate the log-likelihood
    loglikelihood_estimate = loglikelihood_estimate + wmax + log(sum( weights(:,tt) )) - log(nParticles);
    
    % Normalize the weights
    weights(:,tt) = weights(:,tt) / sum( weights(:,tt) );
    
    % Compute state estimate
    state_estimate(tt) = sum( weights(:,tt) .* particles(:,tt) );
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper for computing the logarithm of the Poisson pmf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = dpoisson(x, lambda)
    out = -log( factorial(x) ) + x * log( lambda ) - lambda;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
