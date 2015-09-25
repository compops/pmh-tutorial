%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates data from the LGSS model with parameters (phi,sigmav,sigmae)
%
% (c) 2015 Johan Dahlin
% johan.dahlin (at) liu.se
%
%
% Inputs:
% phi, sigmav, sigmae: the persistence of the state and the 
%                      standard deviations of the state innovations and 
%                      observation noise.
%
% T and xo:            the no. observations and initial state.
%
% Outputs:
% x,y:                 the latent state and observations
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[x,y] = lgss_generateData(phi,sigmav,sigmae,T,xo)

  % Pre-allocate vectors for the state (x) and observations (y)
  x    = zeros(T+1,1);
  y    = zeros(T+1,1);
  
  % Set the initial state
  x(1) = xo;
  
  % Simulate the system for each time step
  for tt = 2:(T+1)
    x(tt) = phi  * x(tt-1) + sigmav  * normrnd(0,1);
    y(tt) = x(tt)          + sigmae  * normrnd(0,1);
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%