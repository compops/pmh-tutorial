%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kalman filtering for LGSS model
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

function xHatFiltered = kalmanFilter(y, theta, initialState, initialStateCovariance)

  T = length(y);

  % Initalise variables
  xHatFiltered = initialState * ones( T, 1);
  xHatPredicted = initialState * ones( T, 1);
  predictiveCovariance = initialStateCovariance;
  
  A = theta(1);
  C = 1;
  Q = theta(2)^2;
  R = theta(3)^2;
      
  %===========================================================
  % Run main loop
  %===========================================================
  for t = 1:T

    % Compute Kalman Gain
    S = C * predictiveCovariance * C + R;
    kalmanGain = predictiveCovariance * C / S;
    
    % Compute state estimate
    xHatFiltered(t)   = xHatPredicted(t) + kalmanGain * ( y(t) - C * xHatPredicted(t) );
    xHatPredicted(t+1) = A * xHatFiltered(t); 
    
    % Update covariance
    filteredCovariance = predictiveCovariance - kalmanGain * S * kalmanGain;
    predictiveCovariance = A * filteredCovariance * A + Q;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
