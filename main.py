#############################################################################
## Main function                                                           ##
#############################################################################

import numpy as np
import pymbar as pm
import pickle

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import container
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size': '12', 'weight':'bold'})
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)
#plt.rc('xtick', labelsize=16)
#plt.rc('ytick', labelsize=16)

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('12')

from mbar_codes.inputs import *
from mbar_codes.loadcolvar import loadcolvar
from mbar_codes.subsample import subsample
from mbar_codes.biases import mbar_biases
from mbar_codes.binning import binning
from mbar_codes.integrate import *

##############################################################################
## Call the Functions and Perform MABR Calculation                          ##
##############################################################################

### Load the history and subsample the correlated history data
percent = 1.0     ### using the full data to compute final PMF
N_k, N_tot_samples, uncorr_sample_idx = subsample(f_prefix, n_windows, n_jobs, subsample_scheme, percent) 

### Load the colvars data and subsample according to the subsampled data indices
x = loadcolvar(f_prefix, n_windows, n_jobs, stride, cv_index, uncorr_sample_idx)

### Create a good number of bins for MBAR calculation
x_avg_bins, num_bins, bins, n_bins = binning(x, num_bins)

### Determine the reduced energy bias
u_kn = mbar_biases(np.vstack(x), k, centers, T).T

### Launch the MBAR calculations
mbar        = pm.MBAR(u_kn, N_k)
pmf_results = mbar.computePMF(u_kn[0].flatten()*0, n_bins, num_bins)

### Integrate the PMF to calculate ensemble average
if integtype == 'sep':
  r_index = np.argmin( np.abs( np.array(x_avg_bins) - x_max ) )
  b_site_index  = np.argmin(pmf_results[0])
  left, right   = b_site_index - 3, b_site_index + 4
  G, dG, integral, sigma = integrate_pmf( np.diff(bins)[0], pmf_results[0][left:right], pmf_results[1][left:right], pmf_results[0][r_index], T)
elif integtype == 'rst':
  x0 = [ centers[ np.argmin(np.abs(x_avg_bins[i] - centers)) ] for i in range(len(x_avg_bins)) ]
  G, dG, integral, sigma = integrate_pmf_potential(x_avg_bins, np.diff(bins)[0], pmf_results[0], pmf_results[1], x0, k, T)

### Prepare for output
results = {}
results['mbar']              = mbar
results['x']                 = np.vstack(x)
results['N_tot_samples']     = N_tot_samples
results['pmf_results']       = pmf_results
results['x_avg_bins']        = x_avg_bins
results['uncorr_sample_idx'] = uncorr_sample_idx
results['num_bins']          = num_bins
results['bins']              = bins
results['G']                 = G
results['dG']                = dG
results['integral']          = integral
results['sigma']             = sigma 

### Plot the PMF profile
ax = plt.subplot()
ms = 3
caps = 1.5
ew   = 0.8
opacity = 0.8
ax.plot(x_avg_bins, R*T*pmf_results[0], color='blue', linewidth=2)
ax.errorbar(x_avg_bins, R*T*pmf_results[0], yerr=R*T*pmf_results[1], ls='none', capsize=caps, color="red", elinewidth=ew, alpha=opacity)
ax.set_xlabel(r'Reaction Coordinate')
ax.set_ylabel(r'PMF (kcal$\cdot$mol$^{-1}$)')
plt.tight_layout()
plt.savefig('PMF-Final.png', figsize=(5,5), dpi=300, bbox_inches='tight')
plt.clf()

#save mbar, np.vstack(x), pmf_results, x_avg_bins, N_samples, uncorr_sample_idx, n_bins, bins as pickle file
pickle.dump(results, open(f_prefix + "-analysis.pkl", "wb"))

##############################################################################
## Check the sampling convergence and plot PMF                              ##
##############################################################################

ax   = plt.subplot2grid((2,2),(0,0))
cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=0, vmax=1)

for i in range(16):
  ### Load the history and subsample the correlated history data
  percent = 0.0625 * (i+1)
  color   = cmap(percent)  

  N_k, N_tot_samples, uncorr_sample_idx = subsample(f_prefix, n_windows, n_jobs, subsample_scheme, percent)

  ### Load the colvars data and subsample according to the subsampled data indices
  x    = loadcolvar(f_prefix, n_windows, n_jobs, stride, cv_index, uncorr_sample_idx)

  ### Create a good number of bins for MBAR calculation
  x_avg_bins, num_bins, bins, n_bins = binning(x, num_bins)

  ### Determine the reduced energy bias
  u_kn = mbar_biases(np.vstack(x), k, centers, T).T

  ### Launch the MBAR calculations
  mbar        = pm.MBAR(u_kn, N_k)
  pmf_results = mbar.computePMF(u_kn[0].flatten()*0, n_bins, num_bins)

  ax.plot(x_avg_bins, R*T*pmf_results[0], color=color, linewidth=2)

ax_inset = plt.axes([.2, .15, .2, .025])
cb1 = mpl.colorbar.ColorbarBase(ax_inset,cmap=cmap,norm=norm,orientation='horizontal')
cb1.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
cb1.set_ticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])
cb1.set_label('Percentage of data usage', fontweight="bold")
ax.set_xlabel(r'Separation $\bf{\it{r}}$ (${\AA}$)', fontweight="bold")
ax.set_ylabel(r"Separation PMF (kcal/mol)", fontweight="bold")
ax.set_ylim(0,50)
ax.set_yticks((0,10,20,30,40,50))
ax.minorticks_on()
plt.tight_layout()
plt.savefig('PMF-Convergence.png', figsize=(10,10), dpi=1000, bbox_inches='tight')
