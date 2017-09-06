%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle Metropolis-Hastings (PMH) for the LGSS model
%
% Johan Dahlin <liu (at) johandahlin.com.nospam>
% Documentation at https://github.com/compops/pmh-tutorial
% Published under GNU General Public License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[phi] = particleMetropolisHastings(observations, initialPhi, parameters, noParticles, initialState, noIterations, stepSize)

  sigmav = parameters(1);
  sigmae = parameters(2);
  
  phi = zeros(noIterations, 1);
  phiProposed = zeros(noIterations, 1);
  logLikelihood = zeros(noIterations, 1);
  logLikelihoodProposed = zeros(noIterations, 1);
  proposedPhiAccepted = zeros(noIterations, 1);
  
  % Set the initial parameter and estimate the initial log-likelihood
  phi(1) = initialPhi;  
  parameters = [phi(1) sigmav sigmae];
  [~, logLikelihood(1)] = particleFilter(observations, parameters, noParticles, initialState);
  
  for k = 2:noIterations  
    % Propose a new parameter
    phiProposed(k) = phi(k-1) + stepSize * normrnd(0, 1);
  
    % Estimate the log-likelihood (don't run if unstable system)
    if (abs(phiProposed(k)) < 1.0)
      thetaProposed = [phiProposed(k), sigmav, sigmae];
      [~, logLikelihoodProposed(k)] = particleFilter(observations, thetaProposed, noParticles, initialState);
    end
    
    % Compute the acceptance probability (reject if unstable system)
    prior = dnorm(phiProposed(k), 0, 1) - dnorm(phi(k - 1), 0, 1);
    likelihoodDifference = logLikelihoodProposed(k) - logLikelihood(k - 1);
    acceptProbability = exp(prior + likelihoodDifference);
    acceptProbability = acceptProbability * (abs(phiProposed(k)) < 1.0);
    
    % Accept / reject step
    uniformRandomVariable = unifrnd(0, 1);
    if (uniformRandomVariable < acceptProbability)
      % Accept the parameter
      phi(k) = phiProposed(k);
      logLikelihood(k) = logLikelihoodProposed(k);
      proposedPhiAccepted(k) = 1.0;
    else
      % Reject the parameter
      phi(k) = phi(k - 1);
      logLikelihood(k) = logLikelihood(k - 1);
      proposedPhiAccepted(k) = 0.0;
    end
    
    % Write out progress
    if ( rem(k, 100) == 0 )  
      disp(['#####################################################################']);
      disp([' Iteration: ',num2str(k),' of : ', num2str(noIterations) ,' completed.']);
      disp([' Current state of the Markov chain: ', num2str(phi(k), 2)]);
      disp([' Proposed next state of the Markov chain: ', num2str(phiProposed(k), 2)]);
      disp([' Current posterior mean: ', num2str(mean(phi(1:k)), 2)]);
      disp([' Current acceptance rate: ', num2str(mean(proposedPhiAccepted(1:k)), 2)]);
      disp(['#####################################################################']);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper for computing the logarithm of N(x; mu, sigma^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log(sigma.^2) - 0.5 ./ sigma.^2 .* (x - mu).^2;
end
