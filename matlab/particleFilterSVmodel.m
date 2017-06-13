%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Bootstrap particle filter for the SV model
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

function[xHatFiltered, logLikelihood] = particleFilterSVmodel(y, theta, noParticles)
  %
  % Bootstrap particle filter for the SV model
  %
  % Inputs:
  % y:                   observations from the system for t=1,...,T.
  %  
  % mu, phi, sigmav:     the mean and persistence of the state and the 
  %                      standard deviations of the state innovations.
  %
  % noParticles:         number of particles (N)
  %
  % Outputs:
  % xHatFiltered:        sample from filtered state distribution
  %
  % logLikelihood:       estimate of the log-likelihood at T-1

  %===========================================================
  % Initialise variables
  %===========================================================
  T = length(y);
  mu = theta(1);
  phi = theta(2);
  sigmav = theta(3);
  
  particles = zeros(noParticles, T + 1);
  ancestorIndices = zeros(noParticles, T + 1);
  weights = ones(noParticles, T + 1);
  normalisedWeights = ones(noParticles, T + 1) / noParticles;
  
  logLikelihood = 0;
  ancestorIndices(:, 1)= 1:noParticles;
  particles(:, 1) = mu + sigmav / sqrt(1 - phi^2) * normrnd(0, 1, noParticles, 1);
  xHatFiltered(1) = mean(particles(:, 1));
  
  %===========================================================
  % Run main loop
  %===========================================================
  for t = 2:(T + 1)
    
    %=========================================================
    % Resample ( multinomial )
    %=========================================================
    newAncestors = randsample(noParticles, noParticles, true, normalisedWeights(:, t - 1)); 

    % Resample the ancestory linage
    ancestorIndices(:, 1:(t - 1)) = ancestorIndices(newAncestors, 1:(t - 1));
    
    % Add the most recent ancestors
    ancestorIndices(:, t) = newAncestors;
    
    %=========================================================
    % Propagate
    %=========================================================
    part1 = mu + phi * (particles(newAncestors, t - 1) - mu);
    part2 = sigmav * normrnd(0, 1, noParticles, 1);
    particles(:, t) = part1 + part2;
    
    %=========================================================
    % Compute weights
    %=========================================================
    weights(:, t) = dnorm(y(t - 1), 0, exp(particles(:, t) / 2));
    
    % Rescale log-weights and recover weight
    maxWeight = max(weights(:, t));
    weights(:, t) = exp(weights(:, t) - maxWeight);

    % Normalize the weights
    sumWeights = sum(weights(:, t));
    normalisedWeights(:, t) = weights(:, t) / sumWeights;    
    
    % Estimate the log-likelihood
    predictiveLikelihood = maxWeight + log(sumWeights) - log(noParticles);
    logLikelihood  = logLikelihood + predictiveLikelihood;
  end
  
  % Sample the state estimate using the weights at tt=T
  xHatFiltered = zeros(1, T + 1);
  ancestorIndex  = randsample(noParticles, 1, true, normalisedWeights(:, T));
  
  for t = 2:(T + 1)
    xHatFiltered(t) = particles(ancestorIndices(ancestorIndex, t), t);
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
