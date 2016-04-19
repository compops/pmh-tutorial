%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Metropolis-Hastings (MH) 
% for the Earthquake data with IID model.
%
% Inputs:
% observations:     observations from the system for t=1,...,T.
%
% initialTheta:     initial values of the parameters.
%
% nIterations:      the number of iterations in the PMH algorithm.
%
% Sigma:            the covariance matrix of the RW proposal.
%
% Outputs:
% theta:             nIterations samples from the parameter posterior.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[ theta ] = mh( observations, initialTheta, nIterations, Sigma )
  
 % Initalise variables
  theta         = zeros( nIterations, length(initialTheta) ); 
  loglikelihood = zeros( nIterations, 1 );
  accept        = zeros( nIterations, 1 );
  
  % Set the initial parameter
  theta(1,:) = initialTheta;  
  
  % Compute the initial log-likelihood
  loglikelihood(1) = sum( dt(observations, theta(1,:)) );
  
  %=====================================================================
  % Run main loop
  %=====================================================================
  for kk = 2:nIterations
    
    % Propose a new parameter
    theta_proposed = mvnrnd( theta(kk-1,:), Sigma );
    
    % Estimate the log-likelihood if DOF and sigma is positive
    if ( ( theta_proposed(1) > 0.0 ) && ( theta_proposed(3) > 0.0 ) )
        loglikelihood_proposed = sum( dt(observations, theta_proposed) );
    else
        loglikelihood_proposed = -inf;
    end
    
    % Compute the acceptance probability
    aprob = exp( loglikelihood_proposed - loglikelihood(kk-1) );
    
    % Generate uniform random variable in U[0,1]
    u = unifrnd(0,1);
    
    % Accept / reject step
    if ( u < aprob )
      
      % Accept the parameter
      theta(kk,:)        = theta_proposed;
      loglikelihood(kk)  = loglikelihood_proposed;
      accept(kk)         = 1.0;
      
    else
      % Reject the parameter
      theta(kk,:)        = theta(kk-1,:);
      loglikelihood(kk)  = loglikelihood(kk-1);
      accept(kk)         = 0.0;
    end
    
    % Write out progress
    if ( rem(kk,5000) == 0 )
      
      disp(['#####################################################################################']);
      disp([' Iteration: ',num2str(kk),' of : ', num2str(nIterations) ,' completed.']);
      disp([' Current state of the Markov chain:       ', num2str(theta(kk,:),3) ]);
      disp([' Proposed next state of the Markov chain: ', num2str(theta_proposed,3) ]);
      disp([' Current posterior mean:                  ', num2str( mean( theta(1:kk,:),1 ) ,3 ) ]);
      disp([' Current acceptance rate:                 ', num2str( mean( accept(1:kk) ),3 ) ]);
      disp(['#####################################################################################']);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%