%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of particle filtering in a LGSS model
%
% Data generation
%
% (c) 2015 Johan Dahlin
% johan.dahlin (at) liu.se
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[x,y] = lgss_generateData(phi,sigmav,sigmae,T,xo)

  % Pre-allocate vectors for log-volatility/state (x) 
  % and log-returns/observations (y)
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