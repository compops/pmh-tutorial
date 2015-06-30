%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of particle filtering in a LGSS model
%
% Bootstrap particle filter
%
% (c) 2015 Johan Dahlin
% johan.dahlin (at) liu.se
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[xhatf,llp] = lgss_sm(y,phi,sigmav,sigmae,nPart,T,x0)

  % Initalise variables
  xhatf = zeros( T+1, 1);
  p     = zeros( nPart, T+1);
  w     = ones(  nPart, T+1) / nPart;
  llp   = 0;
  
  xhatf(1) = x0;
  p(:,1)   = x0;
    
  %===========================================================
  % Run main loop
  %===========================================================
  for tt = 2:(T+1)

    %=========================================================
    % Resample ( multinomial )
    %=========================================================
    idx = randsample( nPart, nPart, true, w(:,tt-1) ); 
    idx = idx( randperm( nPart ) );
    
    %=========================================================
    % Propagate
    %=========================================================
    p(:,tt) = phi * p(idx,tt-1) + sigmav * normrnd( 0, 1, nPart, 1 ) ;
    
    %=========================================================
    % Compute weights
    %=========================================================
    w(:,tt) = dnorm( y(tt-1), p(:,tt), sigmae );
    
    % Rescale log-weights and recover weight
    wmax    = max( w(:,tt) );
    w(:,tt) = exp( w(:,tt) - wmax );
    
    % Estimate the log-likelihood
    llp    = llp + wmax + log(sum( w(:,tt) )) - log(nPart);
    
    % nPartormalize the weights
    w(:,tt) = w(:,tt) / sum( w(:,tt) );
    
    % Estimate the state
    xhatf(tt) = sum( w(:,tt) .* p(:,tt) );
    
  end
  
end

% Helper for computing the logarithm of the Gaussian density
% nPart(x;mu,sigma^2)

function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log( sigma.^2 ) - 0.5 ./ sigma.^2 .* ( x - mu ).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%