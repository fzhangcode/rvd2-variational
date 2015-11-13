# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import numpy as np
import rvd3
import scipy.stats as ss
from scipy.stats import beta
import matplotlib.pyplot as plt

def main():
    control='../22015-09-28_Run_rvd3_synthetic_data_set/hdf5/10/Control.hdf5'
    case='../2015-09-28_Run_rvd3_synthetic_data_set/hdf5/10/Case100_0.hdf5'
    pos = [84,104,124,144,164,184,204,224,244,264,284,304,324,344]  
    for p in pos:
        #bayestest(case, control, p)
        plot(case, control, p)

def bayestest(caseHDF5Name, controlHDF5Name, position):
    alpha=0.05
    tau=0
    caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(caseHDF5Name)
    casegam = caseq['gam']
    controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(controlHDF5Name)
    controlgam = controlq['gam']

    def beta_mean(p):
        return p[0]*1.0/np.sum(p)    

    def beta_var(p):
        s = np.sum(p)
        return p[0]*p[1]/(s**2*(s+1))

    mu = (beta_mean(casegam[position,:]) - casephi['mu0'])- (beta_mean(controlgam[position,:])-controlphi['mu0'])
    sigma = beta_var(casegam[position,:]) + beta_var(controlgam[position,:])
    z = (tau - mu)/sigma
    print (z)
    p = ss.norm.cdf(z)
    print (p[0])
    #bayescall.append(p[0]<alpha)


def plot(caseHDF5Name, controlHDF5Name, position):
    print ("Plot mu of position %d" %(position+1) )

    ###################### plot position for case ###################
    caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(caseHDF5Name)
    casegam = caseq['gam']
    a = casegam[position, 0]
    b = casegam[position, 1]
    #print (a, b)
    
    fig, ax = plt.subplots()
    # display the pdf
    # ppf (percentage point function) is the inverse CDF.
    # median read depth of case file
    cov_case = int(np.median(caseN))
    x_case = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax.plot(x_case, beta.pdf(x_case, a, b), 'b-', lw=4, alpha=0.8, label="Case(VAF=1.0), Depth=%d" %cov_case)
    # generate random variables
    r_case = beta.rvs(a, b, size=1000)
    ax.hist(r_case, normed=True, histtype='stepfilled', alpha=0.2)
    ax.legend(loc='best', frameon=False)
    

    ###################### plot position for control ###################
    controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(controlHDF5Name)
    controlgam = controlq['gam']    
    a = controlgam[position, 0]
    b = controlgam[position, 1]
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
    ax.set_title('$\\beta$ variational distribution of $\mu$ at position %d' %(position+1))
    ax.set_xlabel('$\mu$', fontsize = 20)
    ax.set_ylabel('PDF', fontsize = 18)
       
    # position_VAF_downsample
    plt.savefig('%d_1.0_10.png' %(position+1))
    #plt.show()

if __name__ == '__main__':
    main()
    