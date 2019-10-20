#############################################################################
## Define the function to calculate energy biases                          ##
#############################################################################

import numpy as np
from scipy import constants

def mbar_biases(x, k, centers, T):

    # x: reaction coordinate
    # k: force constanta
    # centers: coordinate center in each umbrella window
    # T: simulation temperature
    # return: 
    #   reduced bias energy of each umbrella window

    R    = (1/constants.calorie)/1000.0 * constants.R
    dim1 = x.shape[0]
    dim2 = centers.shape[0]
    x = x.reshape(dim1,1)
    T = float(T)

    x_array = x * np.ones([dim1,dim2])
    c_array = centers * np.ones([dim1,dim2])
    energy = .5 * k * (x_array - c_array)**2
    return energy / (R * T)
