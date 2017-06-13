%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates data from the LGSS model
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

function[x, y] = generateData(theta, T, initialState)
  %
  % Generates data from the LGSS model
  %
  % Inputs:
  % theta:               the persistence of the state and the 
  % phi, sigmav, sigmae  standard deviations of the state innovations and 
  %                      observation noise.
  %
  % initialState:        the initial state.
  %
  % Outputs:
  % x, y:                the latent state and observations
  %

  % Pre-allocate vectors for the state (x) and observations (y)
  x = zeros(T+1, 1);
  y = zeros(T+1, 1);
  
  % Set the initial state
  x(1) = initialState;
  
  % Extract parameters
  phi = theta(1);
  sigmav = theta(2);
  sigmae = theta(3);
  
  % Simulate the system for each time step
  for t = 2:(T + 1)
    x(t) = phi * x(t-1) + sigmav * normrnd(0, 1);
    y(t) = x(t) + sigmae * normrnd(0, 1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
