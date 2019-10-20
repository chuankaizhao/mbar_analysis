#############################################################################
## Create bins for MBAR calculations                                       ##
#############################################################################

import numpy as np

def binning(x, num_bins):
    # x: input array, assuming 2D
    # num_bins: initial number of bins
    # return: 
    #   x_avg_bins, average value of each bin 
    #   num_bins, final number of bins
    #   bins, actual bins
    #   n_bins, the indices of the bins to which each value in input array belongs

    x_stack = np.sort(np.vstack(x), axis=0)
    min_bin = x_stack[0,0] - (x_stack[1,0] - x_stack[0,0])/2.0
    max_bin = x_stack[-1,0] + (x_stack[-1,0] - x_stack[-2,0])/2.0
    bins    = np.linspace(min_bin, max_bin, num=num_bins+1)
    n_bins  = np.digitize(np.vstack(x).flatten(), bins) - 1
    real_num_bins  = np.unique(n_bins).shape[0]

    # If there are empty bins, increase bin size in order to avoid empty bins
    while real_num_bins != num_bins:
      num_bins -= 1
      bins      = np.linspace(min_bin, max_bin, num=num_bins+1)
      n_bins    = np.digitize(np.vstack(x).flatten(),bins) - 1
      real_num_bins = np.unique(n_bins).shape[0]

    # Calculate the average value of each bin
    x_all      = np.vstack(x).flatten()
    x_avg_bins = []
    for i in range(num_bins):
      bin_avg = x_all[np.where(n_bins.flatten()==i)[0]].mean()
      x_avg_bins.append(bin_avg)

    return x_avg_bins, num_bins, bins, n_bins
