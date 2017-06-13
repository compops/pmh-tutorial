%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Particle Metropolis-Hastings (PMH) for the LGSS model
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

function[phi] = particleMetropolisHastings(y, initialPhi, theta, noParticles, initialState, noIterations, stepSize)
  %
  % Particle Metropolis-Hastings (PMH) for the LGSS model
  %
  % Inputs:
  % y:                   observations from the system for t=1,...,T.
  %
  % initialPhi:          initial value for phi (persistence of the state)
  %
  % sigmav, sigmae:      the standard deviations of the state innovations 
  %                      and observation noise.
  %
  % noParticles:         number of particles (N)
  %
  % initialState:        the initial state.
  %
  % noIterations:        the number of iterations in PMH
  %
  % stepSize:            the standard deviation of the RW proposal.
  %
  % Outputs:
  % phi:                  K samples from the parameter posterior of phi.
  %
  
  %=====================================================================
  % Initialise variables
  %=====================================================================
  sigmav = theta(1);
  sigmae = theta(2);
  
  phi = zeros(noIterations, 1);
  phiProposed = zeros(noIterations, 1);
  logLikelihood = zeros(noIterations, 1);
  logLikelihoodProposed = zeros(noIterations, 1);
  proposedPhiAccepted = zeros(noIterations, 1);
  
  % Set the initial parameter and estimate the initial log-likelihood
  phi(1) = initialPhi;  
  theta = [phi(1) sigmav sigmae];
  [~, logLikelihood(1)] = particleFilter(y, theta,...
                                         noParticles, initialState);
  
  %=====================================================================
  % Run main loop
  %=====================================================================
  for k = 2:noIterations
    
    % Propose a new parameter
    phiProposed(k) = phi(k-1) + stepSize * normrnd(0, 1);
  
    % Estimate the log-likelihood (don't run if unstable system)
    if (abs(phiProposed(k)) < 1.0)
      thetaProposed = [phiProposed(k), sigmav, sigmae];
      [~, logLikelihoodProposed(k)] = particleFilter(y,...
                                                     thetaProposed,... 
                                                     noParticles,... 
                                                     initialState);
    end
    
    % Compute the acceptance probability
    prior = dnorm(phiProposed(k), 0, 1) - dnorm(phi(k - 1), 0, 1);
    likelihoodDifference = logLikelihoodProposed(k) - logLikelihood(k - 1);
    acceptProbability = exp(prior + likelihoodDifference);
    
    % Set acceptance probabilty to zero if the system is unstable
    acceptProbability = acceptProbability * (abs(phiProposed(k)) < 1.0);

    % Generate uniform random variable in U[0,1]
    uniformRandomVariable = unifrnd(0, 1);
    
    % Accept / reject step
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

% Helper for computing the logarithm of the Gaussian density
% N(x;mu,sigma^2.
function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log(sigma.^2) - 0.5 ./ sigma.^2 .* (x - mu).^2;
end
