%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates data from the LGSS model
%
% Johan Dahlin <liu (at) johandahlin.com.nospam>
% Documentation at https://github.com/compops/pmh-tutorial
% Published under GNU General Public License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[states, observations] = generateData(parameters, noObservations, initialState)
  states = zeros(noObservations+1, 1);
  observations = zeros(noObservations+1, 1);
  
  states(1) = initialState;
  phi = parameters(1);
  sigmav = parameters(2);
  sigmae = parameters(3);
  
  for t = 2:(noObservations + 1)
    states(t) = phi * states(t-1) + sigmav * normrnd(0, 1);
    observations(t) = states(t) + sigmae * normrnd(0, 1);
  end
end