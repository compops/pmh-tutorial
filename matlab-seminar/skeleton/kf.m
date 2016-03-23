%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of Kalman filtering for the linear Gaussian state space model
%
% Bootstrap particle filter
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xhatf = kf(y,par,T,x0,P0)

  % Initalise variables
  xhatf = x0 * ones( T+1, 1 );
  xhatp = x0 * ones( T+1, 1 );
  Pp    = P0;
  
  % Set the parameters for the Kalman filter
  A = par(1);
  C = 1;
  Q = par(2)^2;
  R = par(3)^2;
      
  % Run main loop
    for tt = 2:(T+1)

    % Compute Kalman Gain
    S = C * Pp * C + R;
    K = Pp * C / S;
    
    % Compute state estimate
    xhatf(tt)   = xhatp(tt) + K * ( y(tt-1) - C * xhatp(tt) );
    xhatp(tt+1) = A * xhatf(tt); 
    
    % Update covariance
    Pf = Pp - K * S * K;
    Pp = A * Pf * A + Q;
    
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
