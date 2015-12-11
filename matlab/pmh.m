%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Particle Metropolis-Hastings (PMH) for the LGSS model
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
% initPar:             initial value for phi (persistence of the state)
%
% sigmav, sigmae:      the standard deviations of the state innovations 
%                      and observation noise.
%
% nPart:               number of particles (N)
%
% T and xo:            the no. observations and initial state.
%
% nIter and stepSize:  the number of iterations in PMH and the 
%                      standard deviation of the RW proposal.
%
% Outputs:
% th:                  K samples from the parameter posterior.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[th] = pmh(y,initPar,sigmav,sigmae,nPart,T,xo,nIter,stepSize)
  
  %=====================================================================
  % Initialise variables
  %=====================================================================
  th     = zeros(nIter,1);
  thp    = zeros(nIter,1);
  ll     = zeros(nIter,1);
  llp    = zeros(nIter,1);
  accept = zeros(nIter,1);
  
  % Set the initial parameter and estimate the initial log-likelihood
  th(1)     = initPar;  
  [~,ll(1)] = sm(y,th(1),sigmav,sigmae,nPart,T,xo);
  
  %=====================================================================
  % Run main loop
  %=====================================================================
  for kk = 2:nIter
    
    % Propose a new parameter
    thp(kk) = th(kk-1) + stepSize * normrnd(0,1);
  
    % Estimate the log-likelihood (don't run if unstable system)
    if ( abs( thp(kk) ) < 1.0 )
      [~,llp(kk)] = sm(y,thp(kk),sigmav,sigmae,nPart,T,xo);
    end
    
    % Compute the acceptance probability
    aprob = exp( dnorm(thp(kk),0,1) - dnorm(th(kk-1),0,1) + llp(kk) - ll(kk-1) );

    % Generate uniform random variable in U[0,1]
    u = unifrnd(0,1);
    
    % Accept / reject step
    % Check if | phi | > 1.0, in that case always reject.
    if ( (u < aprob) && ( abs(thp(kk)) < 1.0 ) )
      
      % Accept the parameter
      th(kk)     = thp(kk);
      ll(kk)     = llp(kk);
      accept(kk) = 1.0;
      
    else
      % Reject the parameter
      th(kk)     = th(kk-1);
      ll(kk)     = ll(kk-1);
      accept(kk) = 0.0;
    end
    
    % Write out progress
    if ( rem(kk,100) == 0 )
      
      disp(['#####################################################################']);
      disp([' Iteration: ',num2str(kk),' of : ', num2str(nIter) ,' completed.']);
      disp([' Current state of the Markov chain: ', num2str(th(kk),2) ]);
      disp([' Proposed next state of the Markov chain: ', num2str(thp(kk),2) ]);
      disp([' Current posterior mean: ', num2str( mean( th(1:kk) ) ,2 ) ]);
      disp([' Current acceptance rate: ', num2str( mean( accept(1:kk) ),2 ) ]);
      disp(['#####################################################################']);
    end
  end
  
end

% Helper for computing the logarithm of the Gaussian density
% N(x;mu,sigma^2)

function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log( sigma.^2 ) - 0.5 ./ sigma.^2 .* ( x - mu ).^2;
end
