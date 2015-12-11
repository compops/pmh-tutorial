%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of particle filtering in a LGSS model
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

function xhatf = lgss_kf(y,phi,sigmav,sigmae,T,x0,P0)

  % Initalise variables
  xhatf = x0 * ones( T, 1);
  xhatp = x0 * ones( T, 1);
  Pp    = P0;
  
  A = phi;
  C = 1;
  Q = sigmav^2;
  R = sigmae^2;
      
  %===========================================================
  % Run main loop
  %===========================================================
  for tt = 1:T

    % Compute Kalman Gain
    S = C * Pp * C + R;
    K = Pp * C / S;
    
    % Compute state estimate
    xhatf(tt)   = xhatp(tt) + K * ( y(tt) - C * xhatp(tt) );
    xhatp(tt+1) = A * xhatf(tt); 
    
    % Update covariance
    Pf = Pp - K * S * K;
    Pp = A * Pf * A + Q;
    
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
