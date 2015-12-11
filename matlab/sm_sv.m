%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Bootstrap particle filter for the SV model
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
  a     = zeros( N, T+1);
  p     = zeros( N, T+1);
  v     = ones(  N, T+1);
  w     = ones(  N, T+1) / N;
  llp   = 0;
  
  a(:,1)   = 1:N;
  p(:,1)   = mu + sigmav/sqrt(1-phi^2) * normrnd( 0, 1, N, 1 ) ;

  %===========================================================
  % Run main loop
  %===========================================================
  for tt = 2:(T+1)

    %=========================================================
    % Resample ( multinomial )
    %=========================================================
    idx = randsample( N, N, true, w(:,tt-1) ); 
    idx = idx( randperm( N ) );
    
    % Resample the ancestory linage
    a(:,1:tt-1) = a(idx,1:tt-1);
    
    % Add the most recent ancestors
    a(:,tt)     = idx;    
    
    %=========================================================
    % Propagate
    %=========================================================
    p(:,tt) = mu + phi * ( p(idx,tt-1) - mu) + sigmav * normrnd( 0, 1, N, 1 ) ;
    
    %=========================================================
    % Compute weights
    %=========================================================
    v(:,tt) = dnorm( y(tt-1), 0, exp(p(:,tt)/2) );
    
    % Rescale log-weights and recover weight
    vmax    = max( v(:,tt) );
    v(:,tt) = exp( v(:,tt) - vmax );
    
    % Estimate the log-likelihood
    llp    = llp + vmax + log(sum( v(:,tt) )) - log(N);
    
    % Normalize the weights
    w(:,tt) = v(:,tt) / sum( v(:,tt) );
    
  end
  
  % Sample the state estimate using the weights at tt=T
  nIdx  = randsample( N, 1, true, w(:,T) );
  indices = sub2ind(size(p), a(nIdx,:), 1:(T+1));
  xhatf = p( indices );
end

% Helper for computing the logarithm of the Gaussian density
% N(x;mu,sigma^2)

function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log( sigma.^2 ) - 0.5 ./ sigma.^2 .* ( x - mu ).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
