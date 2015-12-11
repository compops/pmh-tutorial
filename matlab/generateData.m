%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates data from the LGSS model with parameters (phi,sigmav,sigmae)
%
% Copyright (C) 2015 Johan Dahlin < johan.dahlin (at) liu.se >
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
