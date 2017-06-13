%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fully-adapted particle filter for the linear Gaussian SSM
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
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[xHatFiltered, logLikelihood] = particleFilter(y, theta, noParticles, initialState)
  %
  % Fully-adapted particle filter for the linear Gaussian SSM
  %
  % Inputs:
  % y:                   observations from the system for t=1,...,T.
  %  
  % phi, sigmav, sigmae: the persistence of the state and the 
  %                      standard deviations of the state innovations and 
  %                      observation noise.
  %
  % noParticles:         number of particles (N)
  %
  % initialState:        initial state.
  %
  % Outputs:
  % xHatFiltered:        vector with T elements
  %                      estimates of the filtered state
  %                      for each t=0,1,...,T-1.
  %
  % logLikelihood:       estimate of the log-likelihood at T-1

  %===========================================================
  % Initialise variables
  %===========================================================
  T = length(y) - 1;
  phi = theta(1);
  sigmav = theta(2);
  sigmae = theta(3);  
  
  particles = zeros(noParticles, T + 1);
  ancestorIndices = zeros(noParticles, T + 1);
  weights = ones(noParticles, T + 1);
  normalisedWeights = ones(noParticles, T + 1) / noParticles;
  xHatFiltered = zeros(T + 1, 1);
  
  logLikelihood = 0;
  ancestorIndices(:, 1)= 1:noParticles;  
  xHatFiltered(1) = initialState;
  particles(:, 1) = initialState;
  
  %===========================================================
  % Run main loop
  %===========================================================
  for t = 2:T

    %=========================================================
    % Resample ( multinomial )
    %=========================================================
    newAncestors = randsample(noParticles, noParticles, true, normalisedWeights(:,t - 1)); 

    % Resample the ancestory linage
    ancestorIndices(:, 1:(t - 1)) = ancestorIndices(newAncestors, 1:(t - 1));
    
    % Add the most recent ancestors
    ancestorIndices(:, t) = newAncestors;
    
    %=========================================================
    % Propagate
    %=========================================================
    part1 = ( sigmav^(-2) + sigmae^(-2) )^(-1);
    part2 = sigmae^(-2) .* y(t);
    part2 = part2 + sigmav^(-2) .* phi .* particles(newAncestors, t - 1);
    particles(:, t) = part1 .* part2 + sqrt(part1) .* normrnd(0, 1, noParticles, 1);
    
    %=========================================================
    % Compute weights
    %=========================================================
    weights(:, t) = dnorm(y(t + 1), phi .* particles(:, t), sqrt(sigmae^2 + sigmav^2));
    
    % Rescale log-weights and recover weight
    maxWeight = max(weights(:, t));
    weights(:, t) = exp(weights(:, t) - maxWeight);

    % Normalize the weights
    sumWeights = sum(weights(:, t));
    normalisedWeights(:, t) = weights(:, t) / sumWeights;    
    
    % Estimate the log-likelihood
    predictiveLikelihood = maxWeight + log(sumWeights) - log(noParticles);
    logLikelihood  = logLikelihood + predictiveLikelihood;
    
    % Estimate the state
    xHatFiltered(t) = mean(particles(:,t));
  end
end

% Helper for computing the logarithm of the Gaussian density
% N(x; mu, sigma^2)

function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log(sigma.^2) - 0.5 ./ sigma.^2 .* (x - mu).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
