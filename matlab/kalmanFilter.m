% Kalman filtering for LGSS model
function xHatFiltered = kalmanFilter(observations, parameters, initialState, initialStateCov)

  noObservations = length(observations);
  A = parameters(1);
  C = 1;
  Q = parameters(2)^2;
  R = parameters(3)^2;
  
  xHatFiltered = initialState * ones( noObservations, 1);
  xHatPredicted = initialState * ones( noObservations, 1);
  predictiveCovariance = initialStateCov;
  
  % Main loop
  for t = 1:noObservations

    % Correction step
    S = C * predictiveCovariance * C + R;
    kalmanGain = predictiveCovariance * C / S;
    filteredCovariance = predictiveCovariance - kalmanGain * S * kalmanGain;
    xHatFiltered(t)   = xHatPredicted(t) + kalmanGain * ( observations(t) - C * xHatPredicted(t) );
    
    % Prediction step
    xHatPredicted(t+1) = A * xHatFiltered(t); 
    predictiveCovariance = A * filteredCovariance * A + Q;
  end
end