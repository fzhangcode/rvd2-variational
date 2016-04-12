# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import numpy as np
import rvd3
import scipy.stats as ss
from scipy.stats import beta
import matplotlib.pyplot as plt

def main():
    f1 = './c4-34_MTH1_M0.hdf5' 
    f2 = './c4-34_MTH1_M0devide10.hdf5'  
    f3 = './c4-34_MTH1_M0times10.hdf5'  
    f4 = './c4-34_MTH1_M0times100.hdf5'  
    #f5 = './c4-34_MTH1_M0times1000.hdf5'  
	
    pos = [339]  
    for p in pos:
        plot(f2, f1, f3, f4, p)

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



def plot(f2HDF5Name, f1HDF5Name, f3HDF5Name, f4HDF5Name, position):

    ###################### plot position for f2 ###################
    f2R, f2N, f2phi, f2q, loc, refb = rvd3.load_model(f2HDF5Name)
    f2gam= f2q['gam']
    a = f2gam[position, 0]
    b = f2gam[position, 1]
    #print (a, b)
    M0_f2 = f2phi['M0']
    
    fig, ax = plt.subplots()
    x_f2 = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax.plot(x_f2, beta.pdf(x_f2, a, b), 'g--', lw=8, alpha=0.8, label="M0=%.4f" %M0_f2)
    
    ###################### plot position for f1 ###################
    f1R, f1N, f1phi, f1q, _, _ = rvd3.load_model(f1HDF5Name)
    f1gam = f1q['gam']    
    a = f1gam[position, 0]
    b = f1gam[position, 1]
    M0_f1 = f1phi['M0']
    x_f1 = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax.plot(x_f1, beta.pdf(x_f1, a, b), 'b-', lw=8, alpha=0.8, label='M0=%.3f' %M0_f1)

    ###################### plot position for f3 ###################
    f3R, f3N, f3phi, f3q, _, _ = rvd3.load_model(f3HDF5Name)
    f3gam = f3q['gam']    
    a = f3gam[position, 0]
    b = f3gam[position, 1]
    M0_f3 = f3phi['M0']
    x_f3 = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax.plot(x_f3, beta.pdf(x_f3, a, b), 'm-.', lw=8, alpha=0.8, label='M0=%.2f' %M0_f3)	
	
    ###################### plot position for f4 ###################
    f4R, f4N, f4phi, f4q, _, _ = rvd3.load_model(f4HDF5Name)
    f4gam = f4q['gam']    
    a = f4gam[position, 0]
    b = f4gam[position, 1]
    M0_f4 = f4phi['M0']
    x_f4 = np.linspace(beta.ppf(0.001, a, b), beta.ppf(0.999, a, b), 100)
    ax.plot(x_f4, beta.pdf(x_f4, a, b), 'yo:', lw=8, alpha=0.8, label='M0=%.1f' %M0_f4)	
	

    legend = ax.legend(loc='upper left')
    for label in legend.get_texts():
        label.set_fontsize(38)	
	
    ax.set_xlabel('$\hat{\mu}_{1,014,740}$', fontsize = 38)
    ax.set_xlim([0.95, 1])    
    plt.setp(plt.gca().get_xticklabels(), fontsize=35)
    plt.setp(plt.gca().get_yticklabels(), fontsize=35)	

    plt.show()

if __name__ == '__main__':
    main()
    
