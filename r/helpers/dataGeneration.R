##############################################################################
# Generating data from a LGSS model
# (c) Johan Dahlin 2017 under MIT license <liu@johandahlin.com.nospam>
##############################################################################
generateData <- function(theta, noObservations, initialState) 
{
  phi <- theta[1] 
  sigmav <- theta[2]
  sigmae <- theta[3]
  
  state <- matrix(0, nrow = noObservations + 1, ncol = 1)
  observation <- matrix(0, nrow = noObservations + 1, ncol = 1)
  
  state[1] <- initialState
  observation[1] <- NA
  
  for (t in 2:(noObservations + 1)) {
    state[t] <- phi * state[t - 1] + sigmav * rnorm(1)
    observation[t] <- state[t] + sigmae * rnorm(1)
  }
  
  list(x = state, y = observation)
}