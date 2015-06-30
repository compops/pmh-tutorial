%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of particle filtering in a LGSS model
%
% Bootstrap particle filter
%
% (c) 2015 Johan Dahlin
% johan.dahlin (at) liu.se
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