# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>
# <codecell>
# Fan Zhang 2015/10/14

import sys
import os
import numpy as np
import scipy as sp
from scipy.stats import beta
import matplotlib.pyplot as plt
import matplotlib
import h5py
import pdb
import rvd27
import rvd3

# Insert the src/python/directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)

# Plot mu of true positive position (344, 84)
# Plot mu of false positive position (38, 263)
position = [263] 
dsample = 10

path_mcmc = "../2015-10-09_Run_rvd2_with_without_chi2_for_comparison_with_rvd3_synthetic_data/hdf5"
control_mcmc='%s/%d/Control.hdf5' %(path_mcmc, dsample)
case_mcmc='%s/%d/Case0_1.hdf5' %(path_mcmc, dsample)

path_var = "../2015-09-28_Run_rvd3_synthetic_data_set/hdf5_ELBOupdate_1"
control_var='%s/%d/Control.hdf5' %(path_var, dsample)
case_var='%s/%d/Case0_1.hdf5' %(path_var, dsample)


def main():
    ################### Read mu of MCMC (rvd2) ################################
    with h5py.File(control_mcmc, 'r') as f:
            muControl = f['mu'][...]
            locControl = f['loc'][...]
    with h5py.File(case_mcmc, 'r') as f:
            muCase = f['mu'][...]
            locCase = f['loc'][...]
    idx = []
    for pos in position:             
            idx.append(pos)   
            muControl1 = muControl[idx]
            muCase1 = muCase[idx]
            #N = 2000
            #(muZ,_,_) =rvd27.sample_post_diff(muCase1, muControl1, N) # sample Z

    ## plot histogram
    num_bins = 25
    for i in xrange(len(position)):
            fig = plt.figure(figsize=(12,8))

            ########### Plot mu of MCMC (rvd2) vs Variational (rvd3) ##################
            # normed=True, the integral of the histogram will sum to 1.
            plt.hist(muCase1[i,:].T, num_bins, normed=True, facecolor='r', alpha=0.5, label = 'Case (MCMC)')
            plt.hist(muControl1[i,:].T, num_bins, normed=True, facecolor='k', alpha=0.5, label = 'Control (MCMC)')

            ############# Plot mu of Variational (rvd3) ################################
            caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(case_var)
            casegam = caseq['gam']
            a = casegam[position, 0]
            b = casegam[position, 1]
            cov_case = int(np.median(caseN))
            x_case = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
            plt.plot(x_case, beta.pdf(x_case, a, b), 'r--', lw=4, alpha=1.0, label="Case (Variational)")
            r_case = beta.rvs(a, b, size=2000)
            plt.hist(r_case, num_bins, normed=True, histtype='stepfilled', alpha=0.2, facecolor='r')

            controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(control_var)
            controlgam = controlq['gam']    
            a = controlgam[position, 0]
            b = controlgam[position, 1]
            x_control = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
            plt.plot(x_control, beta.pdf(x_control, a, b), 'k--', lw=4, alpha=1.0, label='Control (Variational)')
            r_control = beta.rvs(a, b, size=2000)
            plt.hist(r_control, num_bins, normed=True, histtype='stepfilled', alpha=0.2, facecolor='k')

            plt.legend(loc='best', frameon=False)
            plt.xlabel('$\hat{\mu} = \mu-\mu_0$', fontsize = 20)
            plt.xticks(rotation=25) 
            plt.title('$\hat{\mu}$ at position %s when median depth is %d' %((position[i]+1), cov_case), fontsize = 18)
            plt.xticks(rotation=25) 
            plt.savefig('position_%s_%d_mcmc_vs_var.png' %((position[i]+1), cov_case))
            plt.tight_layout()

if __name__ == "__main__":
     main()