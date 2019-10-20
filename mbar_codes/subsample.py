#############################################################################
## Define Subsampling Function                                             ##
#############################################################################

import numpy  as np
import pymbar as pm

def subsample(f_prefix, n_windows, n_jobs, subsample_scheme, percent):

    # f_prefix: prefix name of the output files
    # n_windows: number of umbrella windows
    # n_jobs: number of sequential simulations per umbrella window
    # return:
    #   N_k, number of samples per window
    #   N_tot_samples, total number of samples in subsampled data
    #   uncorr_sample_idx, subsampled data indices

    uncorr_sample_idx = [] 
    N_k  = []
    N_tot_samples = 0

    for i in range(n_windows):
      N_k_window = 0
      data       = []
      for j in range(n_jobs):
        data.append(np.loadtxt(str(i) + '/' + f_prefix + ".job%i.%i.sort.history" % (j,i)))
      data = np.vstack(data)
      data = data[0:int(len(data)*percent)]
      N_tot_samples += len(data)

      ### Chose the subsample scheme
      ###   eq: if the input data is equilibrated trajectory
      ###   noteq: if the input data is not equilibrated
      if subsample_scheme == 'eq':
        sample = pm.timeseries.subsampleCorrelatedData(data[:,3])

      elif subsample_scheme == 'noteq':
        ### identify and discard the equilibration part
        [ t0, g, Neff_max ] = pm.timeseries.detectEquilibration(data[:,3])
        data = data[t0:]
        sample = pm.timeseries.subsampleCorrelatedData(data[:,3],g=g)
        sample = sample + t0
        
      sample = np.array(sample)
      uncorr_sample_idx.append(sample)
      N_k.append(len(sample))
    N_k = np.array(N_k)

    return N_k, N_tot_samples, uncorr_sample_idx
