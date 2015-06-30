%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of particle filtering in a SV model
%
% Bootstrap particle filter
%
% (c) 2015 Johan Dahlin
% johan.dahlin (at) liu.se
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[xhatf,llp] = sv_sm(y,phi,sigma,beta,N,T)

  % Initalise variables
  xhatf = zeros( T+1, 1);
  p     = zeros( N, T+1);
  w     = ones(  N, T+1) / N;
  llp   = 0;
  
  p(:,1)   = sigma/sqrt(1-phi^2) * normrnd( 0, 1, N, 1 ) ;
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
    p(:,tt) = phi * p(idx,tt-1) + sigma * normrnd( 0, 1, N, 1 ) ;
    
    %=========================================================
    % Compute weights
    %=========================================================
    w(:,tt) = dnorm( y(tt-1), 0, beta*exp(p(:,tt)/2) );
    
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