%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Particle Metropolis-Hastings (PMH) for the SV model
%
% Copyright (C) 2017 Johan Dahlin < liu (at) johandahlin.com.nospam >
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[theta, xHatFiltered] = particleMetropolisHastingsSVmodel(y, initialTheta, noParticles, noIterations, stepSize)
  %
  % Particle Metropolis-Hastings (PMH) for the SV model
  %
  % Inputs:
  % y:                   observations from the system for t=1,...,T.
  %
  % initialTheta:        initial value for the parameters (mu, phi, sigmav)
  %
  % noParticles:         number of particles (N)
  %
  %
  % noIterations:        the number of iterations in PMH
  %
  % stepSize:            the standard deviation of the RW proposal.
  %
  % Outputs:
  % theta:               K samples from the parameter posterior of theta.
  %
  
  %=====================================================================
  % Initialise variables
  %=====================================================================
  T = length(y);
  
  theta = zeros(noIterations, 3);
  thetaProposed = zeros(noIterations, 3);
  xHatFiltered = zeros(noIterations, T + 1);  
  xHatFilteredProposed = zeros(noIterations, T + 1);  
  logLikelihood = zeros(noIterations, 1);
  logLikelihoodProposed = zeros(noIterations, 1);
  proposedThetaAccepted = zeros(noIterations, 1);
  
  % Set the initial parameter and estimate the initial log-likelihood
  theta(1, :) = initialTheta;  
  [xHatFiltered(1, :), logLikelihood(1)] = particleFilterSVmodel(y, theta(1, :), noParticles);
  
  %=====================================================================
  % Run main loop
  %=====================================================================
  for k = 2:noIterations
    
    % Propose a new parameter
    thetaProposed(k, :) = mvnrnd(theta(k-1, :), stepSize);
  
    % Estimate the log-likelihood (don't run if unstable system)
    if (abs(thetaProposed(k, 2)) < 1.0) && (thetaProposed(k, 3) > 0.0)
      [xHatFilteredProposed(k, :), logLikelihoodProposed(k)] = particleFilterSVmodel(y, thetaProposed(k, :), noParticles);
    end
    
    % Compute the acceptance probability
    prior = dnorm(thetaProposed(k, 1), 0, 1);
    prior = prior - dnorm(theta(k - 1, 1), 0, 1);
    prior = prior + dnorm(thetaProposed(k, 2), 0.95, 0.05);
    prior = prior - dnorm(theta(k - 1, 2), 0.95, 0.05);
    prior = prior + dgamma(thetaProposed(k, 3), 2, 10);
    prior = prior - dgamma(theta(k - 1, 3), 2, 10);
    
    likelihoodDifference = logLikelihoodProposed(k) - logLikelihood(k - 1);
    acceptProbability = exp(prior + likelihoodDifference);  
    
    % Set acceptance probabilty to zero if the system is unstable
    acceptProbability = acceptProbability * (abs(thetaProposed(k, 2)) < 1.0); 
    acceptProbability = acceptProbability * (thetaProposed(k, 3) > 0.0);
    
    % Generate uniform random variable in U[0,1]
    uniformRandomVariable = unifrnd(0, 1);
    
    % Accept / reject step
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper for computing the logarithm of the Gaussian density
% N(x; mu, sigma^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log(sigma.^2) - 0.5 ./ sigma.^2 .* (x - mu).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper for computing the logarithm of the Gamma density
% Gamma(x; a, b) with mean a/b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = dgamma(x, a, b)
    out = a * log(b) - gammaln(a) + (a-1) * log(x) - b * x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
