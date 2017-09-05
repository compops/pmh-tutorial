from __future__ import print_function, division
import numpy as np
from numpy.random import randn

# Generates data from the LGSS model
def generateData(theta, noObservations, initialState):

    phi = theta[0]
    sigmav = theta[1]
    sigmae = theta[2]

    state = np.zeros(noObservations + 1)
    observation = np.zeros(noObservations)
    state[0] = initialState

    for t in range(1, noObservations):
        state[t] = phi * state[t - 1] + sigmav * randn()
        observation[t] = state[t] + sigmae * randn()

    return(state, observation)
