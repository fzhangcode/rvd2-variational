# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import numpy as np
import rvd3
import scipy.stats as ss
from scipy.stats import beta
import matplotlib.pyplot as plt

def main():
    control='./hdf5/E1/c4-1_MTH1_100.hdf5' #c4-1_MTH1_top400_1986.hdf5
    case='./hdf5/E1/c4-47_MTH1_198605220525.hdf5'  #./hdf5/E1/c4-34_MTH1_top400_20160201.hdf5'

    pos = [339]  # 21, 56, 306, 339, 369
    for p in pos:
        #bayestest(case, control, p)
        plot(case, control, p)

def bayestest(caseHDF5Name, controlHDF5Name, position):
    alpha=0.05
    tau=0
    caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(caseHDF5Name)
    casedelta = caseq['delta']
    controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(controlHDF5Name)
    controldelta = controlq['delta']

    def beta_mean(p):
        return p[0]*1.0/np.sum(p)    

    def beta_var(p):
        s = np.sum(p)
        return p[0]*p[1]/(s**2*(s+1))

    mu = (beta_mean(casedelta[0][position,:]) - casephi['mu0'])- (beta_mean(controldelta[0][position,:])-controlphi['mu0'])
    sigma = beta_var(casedelta[0][position,:]) + beta_var(controldelta[0][position,:])
    z = (tau - mu)/sigma
    print (z)
    p = ss.norm.cdf(z)
    print (p[0])
    #bayescall.append(p[0]<alpha)


def plot(caseHDF5Name, controlHDF5Name, position):
    print ("Plot theta of position %d" %(position) )

    ###################### plot position for case ###################
    caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(caseHDF5Name)
    casedelta = caseq['delta']
    a = casedelta[0][position, 0]
    b = casedelta[0][position, 1]
    #print (a, b)
    
    fig, ax = plt.subplots()
    # display the pdf
    # ppf (percentage point function) is the inverse CDF.
    # median read depth of case file
    cov_case = int(np.median(caseN))
    x_case = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax.plot(x_case, beta.pdf(x_case, a, b), 'b-', lw=4, alpha=0.8, label="Case, Depth=%d" %cov_case)
    # generate random variables
    r_case = beta.rvs(a, b, size=1000)
    ax.hist(r_case, normed=True, histtype='stepfilled', alpha=0.2)
    ax.legend(loc='best', frameon=False)
    

    ###################### plot position for control ###################
    controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(controlHDF5Name)
    controldelta = controlq['delta']    
    a = controldelta[0][position, 0]
    b = controldelta[0][position, 1]
    #print (a, b)
    
    # display the pdf
    # ppf (percentage point function) is the inverse CDF.
    cov_control = int(np.median(controlN))
    x_control = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax.plot(x_control, beta.pdf(x_control, a, b), 'g-', lw=4, alpha=0.8, label='Control, Depth=%d' %cov_control)
    # generate random variables
    r_control = beta.rvs(a, b, size=1000)
    ax.hist(r_control, normed=True, histtype='stepfilled', alpha=0.2)
    ax.legend(loc='best', frameon=False)
    ax.set_title('$\\beta$ variational distribution of $\\theta$ at position %d' %(position))
    ax.set_xlabel('$\\theta$', fontsize = 20)
    ax.set_ylabel('PDF', fontsize = 18)
       
    # position_VAF_downsample
    #plt.savefig('%d.png' %(position))
    plt.show()

if __name__ == '__main__':
    main()
    
