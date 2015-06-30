%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of PMH in a SV model
%
% Particle Metropolis-Hastings sampler
%
% (c) 2015 Johan Dahlin
% johan.dahlin (at) liu.se
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[th] = sv_pmh(y,initPar,nPart,T,nIter,stepSize)
  
 % Initalise variables
  th     = zeros(nIter,3);
  thp    = zeros(nIter,3);
  ll     = zeros(nIter,1);
  llp    = zeros(nIter,1);
  accept = zeros(nIter,1);
  
  % Set the initial parameter and estimate the initial log-likelihood
  th(1,:)   = initPar;  
  [~,ll(1)] = sv_sm(y,th(1,1),th(1,2),th(1,3),nPart,T);
  
  %=====================================================================
  % Run main loop
  %=====================================================================
  for kk = 2:nIter
    
    % Propose a new parameter
    thp(kk,:) = th(kk-1,:) + stepSize .* normrnd(0,1,1,3);
  
    % Estimate the log-likelihood (don't run if unstable system)
    if ( ( abs( thp(kk,1) ) < 1.0 ) && ( thp(kk,2) > 0.0 ) && ( thp(kk,3) > 0.0 ) )
      [~,llp(kk)] = sv_sm(y,thp(kk,1),thp(kk,2),thp(kk,3),nPart,T);
    end
    
    % Compute the acceptance probability
    aprob = exp( llp(kk) - ll(kk-1) );
    
    % Generate uniform random variable in U[0,1]
    u = unifrnd(0,1);
    
    % Accept / reject step
    % Check if | phi | > 1.0, in that case always reject.
    if ( (u < aprob) && ( abs(thp(kk)) < 1.0 ) && ( thp(kk,2) > 0.0 ) && ( thp(kk,3) > 0.0 ))
      
      % Accept the parameter
      th(kk,:)   = thp(kk,:);
      ll(kk)     = llp(kk);
      accept(kk) = 1.0;
      
    else
      % Reject the parameter
      th(kk,:)   = th(kk-1,:);
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