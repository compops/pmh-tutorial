%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap particle filter for the SV model
%
% Johan Dahlin <liu (at) johandahlin.com.nospam>
% Documentation at https://github.com/compops/pmh-tutorial
% Published under GNU General Public License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[xHatFiltered, logLikelihood] = particleFilterSVmodel(observations, parameters, noParticles)

  noObservations = length(observations);
  mu = parameters(1);
  phi = parameters(2);
  sigmav = parameters(3);
  
  particles = zeros(noParticles, noObservations + 1);
  ancestorIndices = zeros(noParticles, noObservations + 1);
  weights = ones(noParticles, noObservations + 1);
  normalisedWeights = ones(noParticles, noObservations + 1) / noParticles;
  
  logLikelihood = 0;
  ancestorIndices(:, 1)= 1:noParticles;
  particles(:, 1) = mu + sigmav / sqrt(1 - phi^2) * normrnd(0, 1, noParticles, 1);
  xHatFiltered(1) = mean(particles(:, 1));
  
  for t = 2:(noObservations + 1)    
    % Resample (multinomial)
    newAncestors = randsample(noParticles, noParticles, true, normalisedWeights(:, t - 1)); 
    ancestorIndices(:, 1:(t - 1)) = ancestorIndices(newAncestors, 1:(t - 1));
    ancestorIndices(:, t) = newAncestors;
    
    % Propagate
    part1 = mu + phi * (particles(newAncestors, t - 1) - mu);
    part2 = sigmav * normrnd(0, 1, noParticles, 1);
    particles(:, t) = part1 + part2;
    
    % Compute weights
    weights(:, t) = dnorm(observations(t - 1), 0, exp(particles(:, t) / 2));
    
    maxWeight = max(weights(:, t));
    weights(:, t) = exp(weights(:, t) - maxWeight);
    sumWeights = sum(weights(:, t));
    normalisedWeights(:, t) = weights(:, t) / sumWeights;    
    
    % Estimate the log-likelihood
    predictiveLikelihood = maxWeight + log(sumWeights) - log(noParticles);
    logLikelihood  = logLikelihood + predictiveLikelihood;
  end
  
  % Sample the state estimate using the weights at t = T
  xHatFiltered = zeros(1, noObservations + 1);
  ancestorIndex  = randsample(noParticles, 1, true, normalisedWeights(:, noObservations));
  
  for t = 2:(noObservations + 1)
    xHatFiltered(t) = particles(ancestorIndices(ancestorIndex, t), t);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper for computing the logarithm of N(x; mu, sigma^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log(sigma.^2) - 0.5 ./ sigma.^2 .* (x - mu).^2;
end