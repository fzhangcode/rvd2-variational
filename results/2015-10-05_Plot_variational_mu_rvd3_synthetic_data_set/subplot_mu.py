# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import numpy as np
import rvd3
import scipy.stats as ss
from scipy.stats import beta
import matplotlib.pyplot as plt

path = "../2015-09-28_Run_rvd3_synthetic_data_set/hdf5"

def main():
    control_10000='%s/10000/Control.hdf5' %path
    control_1000='%s/1000/Control.hdf5' %path
    control_100='%s/100/Control.hdf5' %path
    control_10='%s/10/Control.hdf5' %path
    case_10000='%s/10000/Case1_0.hdf5' %path
    case_1000='%s/1000/Case1_0.hdf5' %path
    case_100='%s/100/Case1_0.hdf5' %path
    case_10='%s/10/Case1_0.hdf5' %path

    pos = [84] #[84,104,124,144,164,184,204,224,244,264,284,304,324,344]  
    for p in pos:
        plot(case_10000, case_1000, case_100, case_10, control_10000, control_1000, control_100, control_10, p)

def beta_mean(p):
    return p[0]*1.0/np.sum(p)
    
    
def get_a_b(p, casegam, center):
    par = casegam[p,:]
    a_est = casegam[p, 0]
    b_est = casegam[p, 1]
    mu_est = par[0]*1.0/np.sum(par)
    mu = mu_est + center
    M = a_est + b_est

    # get the exact alpha and beta 
    a = M*mu
    b = M*(1-mu)
    return a, b
    

def plot(case_10000, case_1000, case_100, case_10, control_10000, control_1000, control_100, control_10, position):
    print ("Plot mu of position %d" %(position+1) )
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    
    xstep = 0.002
    xmax = 0.0152
    ymax = 1200
    size = 15
    ################ Downsample = 10000 #################################################
    caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(case_10000)
    casegam = caseq['gam']
    a, b = get_a_b(position, casegam, casephi['mu0'])
    
    cov_case = int(np.median(caseN))
    x_case = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax1.plot(x_case, beta.pdf(x_case, a, b), 'b-', lw=5, alpha=0.8, label="Case")
    # generate random variables
    r_case = beta.rvs(a, b, size=1000)
    ax1.hist(r_case, normed=True, histtype='stepfilled', alpha=0.2)
    ax1.legend(loc='best', frameon=False)

    controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(control_10000)
    controlgam = controlq['gam']    
    a, b = get_a_b(position, controlgam, controlphi['mu0'])
        
    x_control = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax1.plot(x_control, beta.pdf(x_control, a, b), 'g-', lw=5, alpha=0.8, label='Control')
    # generate random variables
    r_control = beta.rvs(a, b, size=1000)
    ax1.hist(r_control, normed=True, histtype='stepfilled', alpha=0.2)
    ax1.legend(loc='best', frameon=False)
    ax1.set_title('Depth=%d' %cov_case, fontsize = size)
    xticks=np.arange(0, xmax, xstep)
    ax1.set_xticks(xticks)           
    ax1.set_ylim(0,ymax)
    
    print ('mu0^control:', controlphi['mu0'])
    print ('mu0^case', casephi['mu0'], '\n')
    
	 
    #################### Downsample = 1000 #############################################
    caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(case_1000)
    casegam = caseq['gam']
    a, b = get_a_b(position, casegam, casephi['mu0'])
    
    cov_case = int(np.median(caseN))
    x_case = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax2.plot(x_case, beta.pdf(x_case, a, b), 'b-', lw=5, alpha=0.8, label="Case")
    # generate random variables
    r_case = beta.rvs(a, b, size=1000)
    ax2.hist(r_case, normed=True, histtype='stepfilled', alpha=0.2)
    ax2.legend(loc='best', frameon=False)
    
    controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(control_1000)
    controlgam = controlq['gam']    
    a, b = get_a_b(position, controlgam, controlphi['mu0'])
    
    x_control = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax2.plot(x_control, beta.pdf(x_control, a, b), 'g-', lw=5, alpha=0.8, label='Control')
    # generate random variables
    r_control = beta.rvs(a, b, size=1000)
    ax2.hist(r_control, normed=True, histtype='stepfilled', alpha=0.2)
    ax2.legend(loc='best', frameon=False)
    ax2.set_title('Depth=%d' %cov_case, fontsize = size)
    xticks=np.arange(0, xmax, xstep)
    ax2.set_xticks(xticks)
    ax2.set_ylim(0,ymax)

    print ('mu0^control:', controlphi['mu0'])
    print ('mu0^case', casephi['mu0'], '\n')
    
    ################### Downsample = 100 #############################################
    caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(case_100)
    casegam = caseq['gam']
    a, b = get_a_b(position, casegam, casephi['mu0'])
    
    cov_case = int(np.median(caseN))
    x_case = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax3.plot(x_case, beta.pdf(x_case, a, b), 'b-', lw=5, alpha=0.8, label="Case")
    # generate random variables
    r_case = beta.rvs(a, b, size=1000)
    ax3.hist(r_case, normed=True, histtype='stepfilled', alpha=0.2)
    ax3.legend(loc='best', frameon=False)
    
    controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(control_100)
    controlgam = controlq['gam']    
    a, b = get_a_b(position, controlgam, controlphi['mu0'])
    
    x_control = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax3.plot(x_control, beta.pdf(x_control, a, b), 'g-', lw=5, alpha=0.8, label='Control')
    # generate random variables
    r_control = beta.rvs(a, b, size=1000)
    ax3.hist(r_control, normed=True, histtype='stepfilled', alpha=0.2)
    ax3.legend(loc='best', frameon=False)
    ax3.set_title('Depth=%d' %cov_case, fontsize = size)    
    xticks=np.arange(0, xmax, xstep)
    ax3.set_xticks(xticks)
    ax3.set_xlabel('$\mu$', fontsize = size)
    ax3.set_ylim(0,ymax)

    print ('mu0^control:', controlphi['mu0'])
    print ('mu0^case', casephi['mu0'], '\n')
    	
    #################### Downsample = 10 #############################################
    caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(case_10)
    casegam = caseq['gam']
    a, b = get_a_b(position, casegam, casephi['mu0'])
    
    cov_case = int(np.median(caseN))
    x_case = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax4.plot(x_case, beta.pdf(x_case, a, b), 'b-', lw=5, alpha=0.8, label="Case")
    # generate random variables
    r_case = beta.rvs(a, b, size=1000)
    ax4.hist(r_case, normed=True, histtype='stepfilled', alpha=0.2)
    ax4.legend(loc='best', frameon=False)

    controlR, controlN, controlphi, controlq, _, _ = rvd3.load_model(control_10)
    controlgam = controlq['gam']    
    a, b = get_a_b(position, controlgam, controlphi['mu0'])
    
    x_control = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax4.plot(x_control, beta.pdf(x_control, a, b), 'g-', lw=5, alpha=0.8, label='Control')
    # generate random variables
    r_control = beta.rvs(a, b, size=1000)
    ax4.hist(r_control, normed=True, histtype='stepfilled', alpha=0.2)
    ax4.legend(loc='best', frameon=False)
    ax4.set_title('Depth=%d' %cov_case, fontsize = size)	
    xticks=np.arange(0, xmax, xstep)
    ax4.set_xticks(xticks)
    ax4.set_xlabel('$\mu$', fontsize = size)
    ax4.set_ylim(0,ymax)

    print ('mu0^control:', controlphi['mu0'])
    print ('mu0^case', casephi['mu0'], '\n')

    plt.suptitle('$\\beta$ variational distribution of $\mu$ at position %d when VAF=1.0%% ' %(position+1), fontsize = size)
    # manually adjust the spacing of suptitle
    plt.subplots_adjust(top=0.9)
    #plt.tight_layout(fig, rect=[0, 0.03, 1, 0.95]) 
    #plt.show()
    
    fig = plt.gcf()
    fig.set_size_inches(12, 8)
    plt.savefig('mu_%d_VAF=1.0.png' %(position+1))
	
if __name__ == '__main__':
    main()