%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle Metropolis-Hastings (PMH) for the SV model
%
% Johan Dahlin <liu (at) johandahlin.com.nospam>
% Documentation at https://github.com/compops/pmh-tutorial
% Published under GNU General Public License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[theta, xHatFiltered] = particleMetropolisHastingsSVmodel(observations, initialParameters, noParticles, noIterations, stepSize)
  noObservations = length(observations);
  
  theta = zeros(noIterations, 3);
  thetaProposed = zeros(noIterations, 3);
  xHatFiltered = zeros(noIterations, noObservations + 1);  
  xHatFilteredProposed = zeros(noIterations, noObservations + 1);  
  logLikelihood = zeros(noIterations, 1);
  logLikelihoodProposed = zeros(noIterations, 1);
  proposedThetaAccepted = zeros(noIterations, 1);
  
  % Set the initial parameter and estimate the initial log-likelihood
  theta(1, :) = initialParameters;  
  [xHatFiltered(1, :), logLikelihood(1)] = particleFilterSVmodel(observations, theta(1, :), noParticles);
  
  for k = 2:noIterations  
    % Propose a new parameter
    thetaProposed(k, :) = mvnrnd(theta(k-1, :), stepSize);
  
    % Estimate the log-likelihood (don't run if unstable system)
    if (abs(thetaProposed(k, 2)) < 1.0) && (thetaProposed(k, 3) > 0.0)
      [xHatFilteredProposed(k, :), logLikelihoodProposed(k)] = particleFilterSVmodel(observations, thetaProposed(k, :), noParticles);
    end
    
    % Compute the acceptance probability (reject if unstable)
    prior = dnorm(thetaProposed(k, 1), 0, 1);
    prior = prior - dnorm(theta(k - 1, 1), 0, 1);
    prior = prior + dnorm(thetaProposed(k, 2), 0.95, 0.05);
    prior = prior - dnorm(theta(k - 1, 2), 0.95, 0.05);
    prior = prior + dgamma(thetaProposed(k, 3), 2, 10);
    prior = prior - dgamma(theta(k - 1, 3), 2, 10);
    likelihoodDifference = logLikelihoodProposed(k) - logLikelihood(k - 1);
    acceptProbability = exp(prior + likelihoodDifference);  
    acceptProbability = acceptProbability * (abs(thetaProposed(k, 2)) < 1.0); 
    acceptProbability = acceptProbability * (thetaProposed(k, 3) > 0.0);
        
    % Accept / reject step
    uniformRandomVariable = unifrnd(0, 1);
    if (uniformRandomVariable < acceptProbability)
      % Accept the parameter
      theta(k, :) = thetaProposed(k, :);
      xHatFiltered(k, :) = xHatFilteredProposed(k, :);
      logLikelihood(k) = logLikelihoodProposed(k);
      proposedThetaAccepted(k) = 1.0;
    else
      % Reject the parameter
      theta(k, :) = theta(k - 1, :);
      xHatFiltered(k, :) = xHatFiltered(k - 1, :);
      logLikelihood(k) = logLikelihood(k - 1);
      proposedThetaAccepted(k) = 0.0;
    end
    
    % Write out progress
    if ( rem(k, 100) == 0 )
      disp(['#####################################################################################']);
      disp([' Iteration: ',num2str(k),' of : ', num2str(noIterations) ,' completed.']);
      disp([' Current state of the Markov chain:       ', num2str(theta(k, :), 3)]);
      disp([' Proposed next state of the Markov chain: ', num2str(thetaProposed(k, :), 3)]);
      disp([' Current posterior mean:                  ', num2str(mean(theta(1:k, :), 1), 3)]);
      disp([' Current acceptance rate:                 ', num2str(mean(proposedThetaAccepted(1:k)), 3)]);
      disp(['#####################################################################################']);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper for computing the logarithm of N(x; mu, sigma^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log(sigma.^2) - 0.5 ./ sigma.^2 .* (x - mu).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper for computing the logarithm of Gamma(x; a, b) with mean a/b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = dgamma(x, a, b)
    out = a * log(b) - gammaln(a) + (a-1) * log(x) - b * x;
end
