%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Bootstrap particle filter for the SV model
%
% (c) 2015 Johan Dahlin
% johan.dahlin (at) liu.se
%
%
% Inputs:
% y:                   observations from the system for t=1,...,T.
%  
% mu, phi, sigmav:     the mean and persistence of the state and the 
%                      standard deviations of the state innovations.
%
% nPart:               number of particles (N)
%
% T and x0:            the no. observations and initial state.
%
% Outputs:
% xhatf:               vector with T elements
%                      estimates of the filtered state
%                      for each t=0,1,...,T-1.
%
% llp:                 estimate of the log-likelihood at T-1
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[xhatf,llp] = sm_sv(y,mu,phi,sigmav,N,T)

  % Initalise variables
  xhatf = zeros( T+1, 1);
  p     = zeros( N, T+1);
  w     = ones(  N, T+1) / N;
  llp   = 0;
  
  p(:,1)   = mu + sigmav/sqrt(1-phi^2) * normrnd( 0, 1, N, 1 ) ;
  xhatf(1) = mean(p(:,1));

  %===========================================================
  % Run main loop
  %===========================================================
  for tt = 2:(T+1)

    %=========================================================
    % Resample ( multinomial )
    %=========================================================
    idx = randsample( N, N, true, w(:,tt-1) ); 
    idx = idx( randperm( N ) );
    
    %=========================================================
    % Propagate
    %=========================================================
    p(:,tt) = mu + phi * ( p(idx,tt-1) - mu) + sigmav * normrnd( 0, 1, N, 1 ) ;
    
    %=========================================================
    % Compute weights
    %=========================================================
    w(:,tt) = dnorm( y(tt-1), 0, exp(p(:,tt)/2) );
    
    % Rescale log-weights and recover weight
    wmax    = max( w(:,tt) );
    w(:,tt) = exp( w(:,tt) - wmax );
    
    % Estimate the log-likelihood
    llp    = llp + wmax + log(sum( w(:,tt) )) - log(N);
    
    % Normalize the weights
    w(:,tt) = w(:,tt) / sum( w(:,tt) );
    
    % Estimate the state
    xhatf(tt) = sum( w(:,tt) .* p(:,tt) );
    
  end
  
end

% Helper for computing the logarithm of the Gaussian density
% N(x;mu,sigma^2)

function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log( sigma.^2 ) - 0.5 ./ sigma.^2 .* ( x - mu ).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%