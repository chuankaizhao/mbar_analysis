#############################################################################
## Define the function to load colvar data and subsample                   ##
## according to indices                                                    ##
#############################################################################

import numpy as np

def loadcolvar(f_prefix, n_windows, n_jobs, stride, cv_index, uncorr_sample_idx):

    # f_prefix: prefix name of the output files
    # n_windows: number of umbrella windows
    # n_jobs: number of sequential simulations per umbrella window
    # uncorr_sample_idx: subsampled data indices
    # stride: stride for reading the colvars
    # cv_index: the column index of collective variable in colvar file
    # return: x, subsampled colvar data

    x = []
    for i in range(n_windows):
      data = []
      for j in range(n_jobs):
        temp = np.loadtxt(str(i) + '/' + f_prefix + ".job%i.%i.sort.colvars.traj" % (j,i))
        _, idx = np.unique(temp[:,0], return_index=True)
        idx    = np.array(np.sort(idx)[stride-1::stride])
        temp_corr = temp[idx, cv_index:cv_index+1]
        data.append(temp_corr)
      data_uncorr = np.vstack(data)[uncorr_sample_idx[i]]
      x.append(data_uncorr)

    return x
