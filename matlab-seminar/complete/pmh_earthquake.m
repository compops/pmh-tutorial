%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Particle Metropolis-Hastings (PMH) for the Earthquake model
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
% initPar:             initial values of the parameters
%
% nPart:               number of particles (N)
%
% T and xo:            the no. observations and initial state.
%
% nIter and stepSize:  the number of iterations in PMH and the 
%                      covariance matrix of the RW proposal.
%
% Outputs:
% th:                  K samples from the parameter posterior.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[th, xh] = pmh_earthquake( y, initPar, nPart, T, nIter, stepSize )
  
 % Initalise variables
  th     = zeros( nIter, length(initPar) );
  thp    = zeros( nIter, length(initPar) );
  xh     = zeros( nIter, T+1 );  
  xhp    = zeros( nIter, T+1 );  
  ll     = zeros( nIter, 1 );
  llp    = zeros( nIter, 1 );
  accept = zeros( nIter, 1 );
  
  % Set the initial parameter
  th(1,:) = initPar;  
  
  % Run the initial particle filter to estimate the log-likelihood and
  % the latent state
  [ xh(1,:), ll(1) ] = sm_earthquake( y, th(1,:), nPart, T );
  
  %=====================================================================
  % Run main loop
  %=====================================================================
  for kk = 2:nIter
    
    % Propose a new parameter
    thp(kk,:) = mvnrnd( th(kk-1,:), stepSize );
    
    % Estimate the log-likelihood (don't run if unstable system)
    if ( abs( thp(kk,1) ) < 1.0 )
      [ xhp(kk,:), llp(kk) ] = sm_earthquake( y, thp(kk,:), nPart, T );
    end
    
    % Compute the acceptance probability
    aprob = exp( llp(kk) - ll(kk-1) );
    
    % Generate uniform random variable in U[0,1]
    u = unifrnd(0,1);
    
    % Accept / reject step
    % Check if | phi | > 1.0, in that case always reject.
    if ( (u < aprob) && ( abs( thp(kk,1) ) < 1.0 ) )
      
      % Accept the parameter
      th(kk,:)   = thp(kk,:);
      xh(kk,:)   = xhp(kk,:);
      ll(kk)     = llp(kk);
      accept(kk) = 1.0;
      
    else
      % Reject the parameter
      th(kk,:)   = th(kk-1,:);
      xh(kk,:)   = xh(kk-1,:);
      ll(kk)     = ll(kk-1);
      accept(kk) = 0.0;
    end
    
    % Write out progress
    if ( rem(kk,100) == 0 )
      
      disp(['#####################################################################################']);
      disp([' Iteration: ',num2str(kk),' of : ', num2str(nIter) ,' completed.']);
      disp([' Current state of the Markov chain:       ', num2str(th(kk,:),3) ]);
      disp([' Proposed next state of the Markov chain: ', num2str(thp(kk,:),3) ]);
      disp([' Current posterior mean:                  ', num2str( mean( th(1:kk),2 ) ,3 ) ]);
      disp([' Current acceptance rate:                 ', num2str( mean( accept(1:kk) ),3 ) ]);
      disp(['#####################################################################################']);
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%