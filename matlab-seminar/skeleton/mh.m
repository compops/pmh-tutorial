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
  
 % Initialise variables
  theta         = zeros( nIterations, length(initialTheta) ); 
  logposterior  = zeros( nIterations, 1 );
  accept        = zeros( nIterations, 1 );
  
  % Set the initial parameter
  theta(1,:) = initialTheta;  
  
  % Compute the initial log-posterior
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%% ADD CODE HERE %%%%%%%%%%%%%%%%%%%%%
  % Hint: dpossion and observations                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  logposterior(1) = ...
  
  %=====================================================================
  % Run main loop
  %=====================================================================
  for kk = 2:nIterations
    
    % Propose a new parameter
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% ADD CODE HERE %%%%%%%%%%%%%%%%%%%%%
    % Hint: mvnrnd and Sigma                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta_proposed = ...
    
    % Compute the log-posterior if the intensity is positive
    if ( theta_proposed(1) > 0.0 )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% ADD CODE HERE %%%%%%%%%%%%%%%%%%%%%
        % Hint: dpossion and observations                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        logposterior_proposed = ...
          
    else
        logposterior_proposed = -inf;
    end
    
    % Compute the acceptance probability
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% ADD CODE HERE %%%%%%%%%%%%%%%%%%%%%
    % Hint: exp, logposterior_proposed and logposterior
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aprob = ...
    
    % Generate uniform random variable in U[0,1]
    u = unifrnd(0,1);
    
    % Accept / reject step
    if ( u < aprob )
      
      % Accept the parameter
      theta(kk,:)        = theta_proposed;
      logposterior(kk)   = logposterior_proposed;
      accept(kk)         = 1.0;
      
    else
      % Reject the parameter
      theta(kk,:)        = theta(kk-1,:);
      logposterior(kk)   = logposterior(kk-1);
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