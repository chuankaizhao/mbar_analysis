#############################################################################
## Define Constants                                                        ##
#############################################################################

from scipy import constants
import numpy as np

### constants 
R    = (1/constants.calorie)/1000.0 * constants.R      # gas constant in kcal/mol/K
T    = 300                                             # temperature
beta = 1 / (R*T)                                       # inverse temperature
num_bins  = 100                                        # number of bins for MBAR
stride    = 2                                          # stride for reading the colvars, stride is 1 if frequencies for saving history and colvar are the same
subsample_scheme = 'eq'                                # if the input trajectory is pre-equilibrated or not

### variables
n_windows = 45                                         # number of windows
x_min     = 38.                                        # minimum of reaction coordinate
x_max     = 60.                                        # maximum of reaction coordinate
centers   = np.linspace(x_min, x_max, num=n_windows)   # centers of umbrella windows
k         = 10.                                        # force constant
f_prefix  = 'PYL2-HAB1-APO-REUS-SEP'                   # prefix name of the output files
n_jobs    = 1                                          # number of simulation runs per umbrella window
cv_index  = 1                                          # the column index of collective variable in colvar file
subsample_scheme = 'eq'                                # if the input trajectory is pre-equilibrated or not
integtype = 'sep'                                      # 'sep': separation or 'rst': restraint 
