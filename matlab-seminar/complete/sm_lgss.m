%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Particle filter for the linear Gaussian state space model
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
%
% Inputs:
% y:                   observations from the system for t=1,...,T.
%  
% par:                 parameter of the model
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

function[xhatf,llp] = sm_lgss( y, par, nPart, T, x0 )

  %===========================================================
  % Initialise variables
  %===========================================================
  xhatf = zeros( T, 1 );
  p     = zeros( nPart, T );
  w     = ones(  nPart, T) / nPart;
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
    
    %=========================================================
    % Propagate
    %=========================================================
    p(:,tt) = par(1) .* p(idx,tt-1) + par(2) .* normrnd( 0, 1, nPart, 1 ) ;
    
    %=========================================================
    % Compute weights
    %=========================================================
    w(:,tt) = dnorm( y(tt-1), p(:,tt), par(3) );
    
    % Rescale log-weights and recover weight
    wmax    = max( w(:,tt) );
    w(:,tt) = exp( w(:,tt) - wmax );
    
    % Estimate the log-likelihood
    llp    = llp + wmax + log(sum( w(:,tt) )) - log(nPart);
    
    % Normalize the weights
    w(:,tt) = w(:,tt) / sum( w(:,tt) );
    
    % Estimate the state
    xhatf(tt) = sum( w(:,tt) .* p(:,tt) );
    
  end
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper for computing the logarithm of the Gaussian density
% N(x;mu,sigma^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log( sigma.^2 ) - 0.5 ./ sigma.^2 .* ( x - mu ).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
