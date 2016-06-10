%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Particle Metropolis-Hastings (PMH) 
% for the Earthquake data with state space model.
%
% Inputs:
% observations:     observations from the system for t=1,...,T.
%
% initialTheta:     initial values of the parameters.
%
% nParticles:       number of particles (N).
%
% nIterations:      the number of iterations in the PMH algorithm.
%
% Sigma:            the covariance matrix of the RW proposal.
%
% Outputs:
% theta:            nIterations samples from the parameter posterior.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[ theta ] = pmh( observations, initialTheta, nParticles, nIterations, Sigma )
  
  % Hint: modify mh.m by replacing dpoisson with pf to compute the
  % likelihood for the acceptance probability
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
