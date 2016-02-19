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
import rvd27
import rvd3

# Plot mu of true positive position (344, 84)
# Plot mu of false positive position (38, 263)
# Pick two positions for the RVD3 paper:
# 159, dsample=1000, Case0_3; 84, dsample=100, Case1_0
position = [84]
dsample = 100

path_mcmc = "../2015-10-09_Run_rvd2_with_without_chi2_for_comparison_with_rvd3_synthetic_data/hdf5"
control_mcmc='%s/%d/Control.hdf5' %(path_mcmc, dsample)
case_mcmc='%s/%d/Case1_0.hdf5' %(path_mcmc, dsample)

path_var = "../2015-09-28_Run_rvd3_synthetic_data_set/hdf5"
control_var='%s/%d/Control.hdf5' %(path_var, dsample)
case_var='%s/%d/Case1_0.hdf5' %(path_var, dsample)


def beta_mean(p):
    return p[0]*1.0/np.sum(p)
    
def get_a_b(p, casegam, center):
    par = casegam[p,:]
    a_est = casegam[p, 0]
    b_est = casegam[p, 1]
    mu_est = par[0][0]*1.0/np.sum(par)
    #print mu_est, center
    mu = mu_est + center
    #mu = mu_est+0.001
    M = a_est + b_est
    # get the exact alpha and beta 
    a = M*mu
    b = M*(1-mu)
    return a, b

def main():
    ################### Read mu of MCMC (rvd2) ################################
    with h5py.File(control_mcmc, 'r') as f:
            muControl = f['mu'][...]
            locControl = f['loc'][...]
            mu0Control = f['phi/mu0'][()]
    with h5py.File(case_mcmc, 'r') as f:
            muCase = f['mu'][...]
            locCase = f['loc'][...]
            mu0Case = f['phi/mu0'][()]
    idx = []
    for pos in position:             
            idx.append(pos)   
            muControl1 = muControl[idx]
            muCase1 = muCase[idx]
            #N = 2000
            #(muZ,_,_) =rvd27.sample_post_diff(muCase1, muControl1, N) # sample Z

    ## plot histogram
    num_bins = 30
    for i in xrange(len(position)):
            fig = plt.figure(figsize=(25,18))

            ########### Plot mu of MCMC (rvd2) vs Variational (rvd3) ##################
            # normed=True, the integral of the histogram will sum to 1.
            print 'MCMCmu0Case'
            print mu0Case
            print 'MCMCmu0Control'
            print mu0Control
            plt.hist(muCase1[i,:].T + mu0Case, num_bins, normed=True, facecolor='b', alpha=0.5, label = 'Case (MCMC)')
            plt.hist(muControl1[i,:].T + mu0Control, num_bins, normed=True, facecolor='g', alpha=0.5, label = 'Control (MCMC)')
            #plt.hist(muCase1[i,:].T+0.0007, num_bins, normed=True, facecolor='b', alpha=0.5, label = 'Case (MCMC)')
            #plt.hist(muControl1[i,:].T+0.0005, num_bins, normed=True, facecolor='g', alpha=0.5, label = 'Control (MCMC)')
			
            ############# Plot mu of Variational (rvd3) ################################
            caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(case_var)
            casegam = caseq['gam']
            print 'mu0case'
            print casephi['mu0']
            a, b = get_a_b(position, casegam, casephi['mu0'])
            cov_case = int(np.median(caseN))
            x_case = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
            plt.plot(x_case, beta.pdf(x_case, a, b), 'b--', lw=14, alpha=1.0, label="Case (Variational)")
            #r_case = beta.rvs(a, b, size=2000)
            #plt.hist(r_case, num_bins, normed=True, histtype='stepfilled', alpha=0.2, facecolor='r')

            controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(control_var)
            controlgam = controlq['gam']    
            print 'mu0control'
            print controlphi['mu0']
            a, b = get_a_b(position, controlgam, controlphi['mu0'])
            x_control = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
            plt.plot(x_control, beta.pdf(x_control, a, b), 'g--', lw=14, alpha=1.0, label='Control (Variational)')
            #r_control = beta.rvs(a, b, size=2000)
            #plt.hist(r_control, num_bins, normed=True, histtype='stepfilled', alpha=0.2, facecolor='k')

            plt.legend(loc='best', fontsize=45, frameon=False)
            #plt.xlim(0, 0.012)
            plt.xlabel('$\mu_{%s}$' %(position[i]+1), fontsize = 65)
            #plt.xticks(rotation=25) 
            #plt.yticks(rotation=25) 
            plt.setp(plt.gca().get_xticklabels(), fontsize=55)
            plt.setp(plt.gca().get_yticklabels(), fontsize=55)
            #plt.title('$\mu$ at position %s when median depth is %d' %((position[i]+1), cov_case), fontsize = 30)
            #plt.title('Median depth is %d' %(cov_case), fontsize = 50)
            plt.savefig('position_%s_%d_mcmc_vs_var_mu_fig1.png' %((position[i]+1), cov_case))
            plt.tight_layout()
            #plt.show()

if __name__ == "__main__":
     main()