#############################################################################
## Integrate the PMF functions                                             ##
#############################################################################

import numpy as np
from scipy import constants

def integrate_pmf(dx, W, dW, W0, T):
    # dx: bin size
    # W : PMF value
    # dW: PMF error value
    # W0: reference state
    # T : temperature
    
    R    = (1/constants.calorie)/1000.0 * constants.R 
    beta = 1/(R*T)
    W    = R*T*W
    dW   = R*T*dW
    W0   = R*T*W0

    integral = np.sum(np.exp(-beta*(W-W0))*dx)
    sigma_sq = np.sum(((beta*np.exp(-beta*(W-W0))*dx)**2)*(dW**2))
    sigma    = sigma_sq**0.5
    G        = -R*T*np.log(integral)
    dG       = ((R*T)**2*sigma_sq/integral**2)**0.5

    return G, dG, integral, sigma


def integrate_pmf_potential(x, dx, W, dW, x0, k, T):
    # x : 
    # dx: bin size
    # W : PMF value
    # dW: PMF error value
    # x0: centers
    # k : force constant
    # T : temperature

    R    = (1/constants.calorie)/1000.0 * constants.R 
    beta = 1/(R*T)
    x    = np.array(x)
    W    = R*T*W
    dW   = R*T*dW
    beta = 1/(R*T)

    integral = np.sum(np.exp(-beta*W)*np.exp(-beta*(.5*k*(x-x0)**2))*dx)/np.sum(np.exp(-beta*W)*dx)
    p        = (np.exp(-beta*W)*dx)/np.sum(np.exp(-beta*W)*dx)
    sigma_sq = np.sum((beta*p*(integral-np.exp(-beta*0.5*k*(x-x0)**2)))**2*(dW**2))
    sigma    = sigma_sq**0.5
    G        = -R*T*np.log(integral)
    dG       = ((R*T)**2*sigma_sq/integral**2)**0.5

    return G, dG, integral, sigma
