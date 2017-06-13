##############################################################################
#
# Generating data from a LGSS model
#
# Copyright (C) 2017 Johan Dahlin < liu (at) johandahlin.com.nospam >
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
##############################################################################

from __future__ import print_function, division

import numpy as np
from numpy.random import randn


##############################################################################
##############################################################################
#
# Generates data from the LGSS model with parameters (phi, sigmav, sigmae)
#
# Inputs:
# theta:               the persistence of the state theta[0] and the
#                      standard deviations of the state noise theta[1] and
#                      observation noise theta[2].
#
# T and initialState:  the no. observations and initial state.
#
# Outputs:
# x, y:                the latent state and observations
#
#
##############################################################################
##############################################################################

def generateData(theta, T, initialState):
    # Pre-allocate vectors for state (x)
    # and measurements (y)
    x = np.zeros(T + 1)
    y = np.zeros(T)

    # Set the initial state
    x[0] = initialState

    # Simulate for each time step
    for t in range(1, T):
        x[t] = theta[0] * x[t - 1] + theta[1] * randn()
        y[t] = x[t] + theta[2] * randn()

    return(x, y)
