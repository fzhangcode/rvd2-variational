# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import numpy as np

import scipy.stats as ss
import scipy.optimize as so
from scipy.special import gammaln, psi, betaln
from scipy import linalg, integrate

#import pandas as pd

import multiprocessing as mp
from itertools import repeat

import h5py
import tempfile
import logging
import time
from datetime import datetime
from datetime import date
import warnings
import pdb
import re
from timeit import default_timer as timer

def main():
    log_level = logging.DEBUG # default logging level
    logging.basicConfig(level=log_level, format='%(levelname)s:%(module)s:%(message)s')

    ## Generate simulation data
    J = 10 # number of positions
    N = 3 # number of replicates
    
    n = np.empty([N,J], dtype=np.int64) # read depth (rep, location)
    n.fill(1000)
    
    phi = {'M0':100, 'mu0':0.1, 'M':[1000]*J}
    (r, theta, mu) = generate_sample(phi, N, J, n, seedint=19891129)

    ## model optimization  
    (phiHat, qHat) = ELBO_opt(r, n, seed = 19891129, pool = 40)

    ## save the parameters.
    save_model('case_model.hdf5', r, n, phiHat, qHat)


    # (r, theta, mu) = generate_sample(phi, N, J, n, seedint=19900906)

    # ## model optimization  
    # (phiHat, qHat) = ELBO_opt(r, n, seed = 19900906, pool = 40)
    # save_model('control_model.hdf5', r, n, phiHat, qHat)

    # test('case_model.hdf5','control_model.hdf5')


def test(caseHDF5Name, controlHDF5Name, alpha=0.05, tau=0, chi2=False, outputFile=None):
    # pdb.set_trace()
    caseR, caseN, casephi, caseq, loc, refb = load_model(caseHDF5Name)
    casegam = caseq['gam']
    controlR, controlN, controlphi, controlq, _, _ = load_model(controlHDF5Name)
    controlgam = controlq['gam']

    (N,J) = np.shape(caseR)[0:2]

    def beta_mean(p):
        return p[0]*1.0/np.sum(p)    

    def beta_var(p):
        s = np.sum(p)
        return p[0]*p[1]/(s**2*(s+1))
    
    # pdb.set_trace()
    bayescall = []    
    for j in xrange(J):
        mu = (beta_mean(casegam[j,:]) - casephi['mu0'])- (beta_mean(controlgam[j,:])-controlphi['mu0'])
        sigma = np.sqrt(beta_var(casegam[j,:]) + beta_var(controlgam[j,:]))
        z = (tau - mu)/sigma
        p = ss.norm.cdf(z)
        # pdb.set_trace()
        bayescall.append(p[0]<alpha)
    # logging.debug(call)
    # pdb.set_trace()

    ## combine the chi2 goodness of fit test
    if chi2:
        chi2call, chi2P = chi2combinetest(caseR, caseN, bayescall)
        call = np.logical_and(bayescall,chi2call)
    else:
        call = bayescall
    # pdb.set_trace()

    if outputFile is not None:

        vcfFilename = outputFile+'.vcf'

        write_dualvcf(vcfFilename, loc, call, refb, controlR, controlN, caseR, caseN)
        # output hdf5 file
        h5Filename = outputFile +'.hdf5'
        h5file = h5py.File(h5Filename, 'w')

        h5file.create_dataset('call', data=call)
        h5file.create_dataset('refb', data=refb)
        h5file.create_dataset('loc', data=loc, 
                              chunks=True, fletcher32=True, compression='gzip')               
        h5file.create_dataset('controlN', data=controlN, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('caseN', data=caseN, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('controlR', data=controlR, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('caseR', data=caseR, 
                              chunks=True, fletcher32=True, compression='gzip')
        if chi2:
            h5file.create_dataset('chi2call', data=chi2call, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('bayescall', data=bayescall)
        h5file.close()


    ## output the results
def write_dualvcf(outputFile, loc, call, refb, controlR=None, controlN=None, caseR=None, caseN=None):

    controlR = np.median(controlR,0)
    caseR = np.median(caseR,0)
    '''
        Write high confidence variant calls from somatic test when there are both control and case sample to VCF 4.2 file.
    '''
    J = len(loc)
    
    today=date.today()
    
    chrom = [x.split(':')[0][3:] for x in loc]
    pos = [int(x.split(':')[1]) for x in loc]
    
    vcfF = open(outputFile,'w')
    
    print("##fileformat=VCFv4.1", file=vcfF)
    print("##fileDate=%0.4d%0.2d%0.2d" % (today.year, today.month, today.day), file=vcfF)

    print("##source=rvd2", file=vcfF)

    print('##PosteriorTestSample= control-case-paired_sample.', file=vcfF)


    # print("##Posterior difference threshold = %0.2f" %tau, file=vcfF)
    # print("##Probability threshold alpha = %0.2f" %alpha, file=vcfF)
    
    # print("##Chi square test is included", file=vcfF)

    uniquechrom = set(chrom)
    uniquechrom = list(uniquechrom)

    for i in xrange(len(uniquechrom)):
        seq = [idx for idx, name in enumerate(chrom) if name==uniquechrom[i]]
        seqlen = len(seq)
        print("##contig=<ID=%(chr)s,length=%(seqlen)d>" %{'chr': uniquechrom[i],'seqlen': seqlen}, file=vcfF)

    
    print("##INFO=<ID=COAF,Number=1,Type=Float,Description=\"Control Allele Frequency\">", file=vcfF)
    print("##INFO=<ID=CAAF,Number=1,Type=Float,Description=\"Case Allele Frequency\">", file=vcfF)

    print("##FORMAT=<ID=AU,Number=1,Type=Integer,Description=\"Number of 'A' alleles\">", file=vcfF)
    print("##FORMAT=<ID=CU,Number=1,Type=Integer,Description=\"Number of 'C' alleles\">", file=vcfF)
    print("##FORMAT=<ID=GU,Number=1,Type=Integer,Description=\"Number of 'G' alleles\">", file=vcfF)
    print("##FORMAT=<ID=TU,Number=1,Type=Integer,Description=\"Number of 'T' alleles\">", file=vcfF)
    
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNormal\tCase", file=vcfF)

    for i in xrange(J):
        # pdb.set_trace()
        if call[i]:           
            # restore R
            actg = ['A','C','G','T']

            idx = actg.index(refb[i])
            caseR4 = np.zeros(4)
            controlR4 = np.zeros(4)
            caseR4[idx] = np.median(caseN[:,i])-np.sum(caseR[i,:])
            controlR4[idx] = np.median(controlN[:,i])-np.sum(controlR[i,:])
            for d in xrange(idx):
               caseR4[d] = caseR[i,d]
               controlR4[d] = controlR[i,d]
            for d in xrange(3-idx):
               caseR4[d+idx+1] = caseR[i,d+idx]
               controlR4[d+idx+1] = controlR[i,d+idx]

            print ("chr%s\t%d\t.\t%s\t.\t.\tPASS\t.\tAU:CU:GU:TU\t%d:%d:%d:%d\t%d:%d:%d:%d" \
                 % (chrom[i], pos[i], refb[i],\
                    controlR4[0], controlR4[1], controlR4[2], controlR4[3],\
                    caseR4[0], caseR4[1], caseR4[2], caseR4[3]), file=vcfF)
                
    vcfF.close()

    '''print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", file=vcfF)

    for i in xrange(J):
        # pdb.set_trace()
        if call[i]:           
            # restore R
            actg = ['A','C','G','T']

            idx = actg.index(refb[i])
            caseR4 = np.zeros(4)
            controlR4 = np.zeros(4)
            caseR4[idx] = np.median(caseN[:,i])-np.sum(caseR[i,:])
            controlR4[idx] = np.median(controlN[:,i])-np.sum(controlR[i,:])
            for d in xrange(idx):
                caseR4[d] = caseR[i,d]
                controlR4[d] = controlR[i,d]
            for d in xrange(3-idx):
                caseR4[d+idx+1] = caseR[i,d+idx]
                controlR4[d+idx+1] = controlR[i,d+idx]

            print ("chr%s\t%d\t.\t%s\t.\t.\tPASS\t.\tAU:CU:GU:TU" % (chrom[i], pos[i], refb[i]), file=vcfF)
                
    vcfF.close()'''

def chi2combinetest(R, N, bayescall = 1, pvalue = 0.05):

    nRep = R.shape[0]
    J = R.shape[1]
    chi2Prep = np.zeros((J,nRep))
    chi2P = np.zeros((J,1))
    for j in xrange(J):
            chi2Prep[j,:] = np.array([chi2test(R[i,j,:] ) for i in xrange(nRep)] )
            if np.any(np.isnan(chi2Prep[j,:])):
                chi2P[j] = np.nan
            else:
                chi2P[j] = 1-ss.chi2.cdf(-2*np.sum(np.log(chi2Prep[j,:] + np.finfo(float).eps)), 2*nRep) # combine p-values using Fisher's Method

    nbayescall = sum(bayescall)
    if nbayescall < 1:
        nbayescall = 1

    if np.median(N) > 500: #Benjamini-Hochberg method FWER control
        chi2call = chi2P < pvalue/nbayescall
    else:
        chi2call = chi2P < pvalue

    chi2call = chi2call.flatten()
    chi2P = chi2P.flatten()

    return  chi2call, chi2P


def chi2test(X, lamda=2.0/3, pvector=np.array([1.0/3]*3)):
    """ Do chi2 test to decide how well the error reads fits uniform multinomial distribution. P-value returned.
        lamda=1 Pearson's chi-square
        lamda=0 the log likelihood ratio statistic/ G^2
        lamda=-1/2 Freeman-Tukey's F^2
        lamda=-1  Neyman modified chi-square
        lamda=-2  modified G^2
    """
    X=np.array(X)

    nsum=np.sum(X)
    if nsum == 0: return np.nan # return NaN if there are no counts
    E=nsum*pvector


    if lamda==0 or lamda==-1:
        C=2.0*np.sum(X*np.log(X*1.0/E))
    else:
        C=2.0/(lamda*(lamda+1))*np.sum(X*((X*1.0/E)**lamda-1))
        
    df=len(pvector)-1
    #p=scipy.special.gammainc(C,df)
    # p=1-gammainc(df/2,C/2)
    p = 1 - ss.chi2.cdf(C, df) 
    return(p)
 

    
def generate_sample(phi, N=3, J=100, n=100, seedint=None):
    """Returns a sample with n reads, N replicates, and
    J locations. The parameters of the model are in the structure phi.
    """
    
    if seedint is not None: 
        np.random.seed(seedint)
    
    #TODO: test for size of n and make an array if a scalar
    
    # Draw J location-specific error rates from a Beta
    alpha0 = phi['M0']*phi['mu0']
    beta0 = phi['M0']*(1-phi['mu0'])
    mu = ss.beta.rvs(alpha0, beta0, size=J)
    
    # Draw sample error rate and error count
    theta=np.zeros((N,J))
    r = np.zeros((N,J))
    for j in xrange(0, J):
        alpha = mu[j]*phi['M'][j]
        beta = (1-mu[j])*phi['M'][j]
        theta[:,j] = ss.beta.rvs(alpha, beta, size=N)
        r[:,j] = ss.binom.rvs(n[:,j], theta[:,j])
    return r, theta, mu
    

## compute sufficient statistics
def EqlogTheta(delta):
    if delta[0] < np.finfo(float).eps:
        delta[0] += np.finfo(float).eps
    return psi(delta[0]) - psi(np.sum(delta))

def Eqlog1_Theta(delta):
    if delta[1] < np.finfo(float).eps:
        delta[1]+=np.finfo(float).eps
    return psi(delta[1]) - psi(np.sum(delta))

def EqMu(gam):
    return gam[0] / (np.sum(gam)) # eps?

    # with warnings.catch_warnings():
    #     warnings.filterwarnings('error')
    #     # pdb.set_trace()
    #     try:
    #         x = gam[0] / (np.sum(gam)) # eps?
    #     except RuntimeWarning: 
    #         # print 'Raised!'
    #         pdb.set_trace()


def EqlogMu(gam):
    if gam[0] < np.finfo(float).eps:
        gam[0] += np.finfo(float).eps
    return psi(gam[0]) - psi(np.sum(gam))    

def Eqlog1_Mu(gam):
	return psi(gam[1]) - psi(np.sum(gam))

def EqlogGamma(gam, M): 
    # Expectation of Beta function coefficient
    # logGamma = integrate.quad(kernel, np.finfo(float).eps, 1-np.finfo(float).eps, args=(gam, M), full_output=1)
    logGamma = integrate.quad(kernel, 1e-3, 1-1e-3, args=(gam, M), full_output=1)
    return logGamma[0]

def kernel(mu, gam, M): 
    return -ss.beta.pdf(mu, gam[0], gam[1])*betaln(mu*M, (1-mu)*M)
    
'''
def kernel(mu, gam, M): 
    a1 = mu*M
    a2 = (1-mu)*M
    # Use Normal distribution to approximate the Beta distribution when a1 and a2 are sufficiently large. 
    # beta (a1, a2)=normal(a1/(a1+a2), (a1*a2/((a1+a2)**2*(a1+a2+1)))**(0.5))
    #return ss.beta.pdf(mu, gam[0], gam[1]) * ss.norm.pdf(a1/(a1+a2), (a1*a2/((a1+a2)**2*(a1+a2+1)))**(0.5)) 
    return ss.norm.pdf(a1/(a1+a2), (a1*a2/((a1+a2)**2*(a1+a2+1)))**(0.5))     '''

## compute entropy
def BetaEntropy(x):
	# To compute EqlogQmu and EqlogQtheta
    return betaln(x[0], x[1]) - (x[0]-1) * psi(x[0]) - (x[1] - 1) * psi(x[1]) + (x[0] + x[1] -2) * psi(x[0] + x[1])

## compute ELBO
def ELBO(r, n, M, mu0, M0, delta, gam):
    if np.ndim(r) == 1: 
        N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: 
        N, J = np.shape(r)

    # Compute the expectations  
    try:
        Mu = np.array([EqMu(gam[j,:]) for j in xrange(J)])
    except TypeError:
        pdb.set_trace()
    # Mu = np.array([EqMu(gam[j,:]) for j in xrange(J)])
    logMu = np.array([EqlogMu(gam[j,:]) for j in xrange(J)])
    log1_Mu = np.array([Eqlog1_Mu(gam[j,:]) for j in xrange(J)])

    logTheta = np.zeros((N,J))
    log1_Theta = np.zeros((N,J))

    for j in xrange(J):
        for i in xrange(N):
            logTheta[i,j] = EqlogTheta(delta[i,j,:])
            log1_Theta[i,j] = Eqlog1_Theta(delta[i,j,:])

    # Eq[log p(r|theta, n)]
    EqlogPr = 0.0
    for j in xrange(J):
        for i in xrange(N):
            EqlogPr += -betaln(r[i,j] + 1, n[i,j] - r[i,j] +1)-np.log(n[i,j]+1) 
            EqlogPr += r[i,j]*logTheta[i,j] + (n[i,j] - r[i,j]) * log1_Theta[i,j] 

    # Eq[log p(theta|mu, M)]
    EqlogPtheta = 0.0
    for j in xrange(J):
        
        EqlogPtheta += N*EqlogGamma(gam[j,:], M[j])
        for i in xrange(N):
            EqlogPtheta += (M[j]* Mu[j]- 1)*logTheta[i,j] +\
            (M[j]*(1 - Mu[j]) - 1)*log1_Theta[i,j]
    # Eq[log p(mu; mu0, M0)]
    EqlogPmu = -J * betaln(mu0*M0, (1-mu0)*M0)
    for j in xrange(J):
        EqlogPmu += (M0*mu0-1)*logMu[j] + (M0*(1-mu0)-1)*log1_Mu[j]

    EqlogQtheta = 0.0
    for j in xrange(J):
        for i in xrange(N):
            EqlogQtheta -= BetaEntropy(delta[i,j,:])

    EqlogQmu = 0.0
    for j in xrange(J):
        EqlogQmu -= BetaEntropy(gam[j,:])

    return EqlogPr + EqlogPtheta + EqlogPmu - EqlogQtheta - EqlogQmu
    # return EqlogPr + EqlogPtheta + EqlogPmu

def ELBO_delta_ij(r, n, M, delta, gam):
    ## partial ELBO from replicate i position j
    ## ELBO used to optimize delta
    ## Commented out all items that don't depend on delta
  
    Mu = EqMu(gam)
    logTheta = EqlogTheta(delta)
    log1_Theta = Eqlog1_Theta(delta)

    # EqlogPr = -betaln(r+1, n-r+1) -np.log(n+1) + r*logTheta + (n - r)*log1_Theta
    # EqlogPtheta = EqlogGamma(gam, M) + (M*Mu - 1)*logTheta + (M*(1-Mu)-1)*log1_Theta
    # EqlogQtheta = BetaEntropy(delta)

    EqlogPr = r*logTheta + (n - r)*log1_Theta

    EqlogPtheta = (M*Mu - 1)*logTheta + (M*(1-Mu)-1)*log1_Theta

    EqlogQtheta = -BetaEntropy(delta)

    return EqlogPr + EqlogPtheta - EqlogQtheta
    # return EqlogPr + EqlogPtheta

def neg_ELBO_delta_ij(logdelta, gam, r, n, M):
    return -ELBO_delta_ij(r, n, M, np.exp(logdelta), gam)

def opt_delta_ij(args):
    # pdb.set_trace()
    r, n, M, delta, gam = args
    # pdb.set_trace()
    # limit delta to [0.001, 1000], np.log(delta) is [-6.9, 6.9]
    #bnds = [[-7, 7]]*2 
    # limit delta to [0.0001, 10000], np.log(delta) is [-10, 10]
    bnds = [[-10, 10]]*2 
    args=(gam, r, n, M)

    # logging.debug(bnds)
    # logging.debug(np.log(delta))
    logdelta = opt_par(neg_ELBO_delta_ij, np.log(delta), args, bnds, 'delta')
    delta = np.exp(logdelta)

    return delta

def opt_delta(r, n, M, delta, gam, pool = None):
    logging.debug("Optimizing delta")

    if np.ndim(r) == 1: N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: N, J = np.shape(r)

    # pdb.set_trace()
    st = time.time()
    if pool is not None:
        for i in xrange(N):
            # pdb.set_trace()
            args = zip (r[i,:], n[i,:], M, delta[i,:], gam)
            temp  = pool.map(opt_delta_ij, args)
            delta[i,:] = np.array(temp)
    else:
        logging.debug('Optimizing delta in single thread')
        for i in xrange(N):
            for j in xrange(J):
                logging.debug('Optimizing position %d of %d and replicate %d of %d' % (j,J,i,N))
                args = (r[i,j], n[i,j], M[j], delta[i,j,:], gam[j,:])
                delta[i,j,:] = opt_delta_ij(args)
                # logging.debug(delta[i,j
    logging.debug('Delta update elapsed time is %0.3f sec for %d samples %d replicates.' % (time.time() - st, J, N))
    return delta

def ELBO_gam_j( M, mu0, M0, delta, gam):
    ## partial ELBO depending on gam from each position j
    ## ELBO used to gam

    if np.ndim(delta) == 1: N = 1
    elif np.ndim(delta) == 2: N= np.shape(delta)[0]

    Mu = EqMu(gam)
    logMu = EqlogMu(gam)
    log1_Mu = Eqlog1_Mu(gam)

    logTheta = np.zeros((N,1))
    log1_Theta = np.zeros((N,1))

    for i in xrange(N):
        logTheta[i] = EqlogTheta(delta[i,:])
        log1_Theta[i] = Eqlog1_Theta(delta[i,:])

    EqlogPtheta = N*EqlogGamma(gam,M)
    for i in xrange(N):
        EqlogPtheta += (M*Mu-1) * logTheta[i] + (M*(1-Mu)-1)*log1_Theta[i] ## I had a typo here (M->Mu)

    EqlogPmu= -betaln(mu0*M0, (1-mu0)*M0)+ (M0*mu0-1)*logMu + (M0*(1-mu0)-1)*log1_Mu

    EqlogQmu = -BetaEntropy(gam)

    return  EqlogPtheta + EqlogPmu - EqlogQmu
    # return  EqlogPtheta + EqlogPmu

def neg_ELBO_gam_j(loggam, delta, M, mu0, M0):
    return -ELBO_gam_j(M, mu0, M0, delta, np.exp(loggam))

def opt_gam_j(args):
    M, mu0, M0, delta, gam = args
    # pdb.set_trace()
    # limit gam to [0.001, 1000], np.log(gam) is [-6.9, 6.9]
    #bnds = [[-7, 7]]*2 
    # limit gam to [0.0001, 10000], np.log(gam) is [-10, 10]
    bnds = [[-10, 10]]*2 
    args = (delta, M, mu0, M0)

    # def opt_par(func, x, args, bnds, parlabel):
    loggam = opt_par(neg_ELBO_gam_j, np.log(gam), args, bnds, 'gamma')
    gam = np.exp(loggam)
    # logging.debug(bnds)
    # logging.debug(loggam)
    return gam

def opt_gam(M, mu0, M0, delta, gam, pool = None):
    logging.debug("Optimizing gam")

    if np.ndim(gam) == 1: J=1
    elif np.ndim(gam) == 2: J=np.shape(gam)[0]

    st = time.time()

    if pool is not None:
        args = zip( M, repeat(mu0,J), repeat(M0,J), np.transpose(delta,axes=(1,0,2)),gam)
        gam = pool.map(opt_gam_j, args)
        gam = np.array(gam)

    else:
        for j in xrange(J):
            # pdb.set_trace()
            logging.debug("Optimizing gamma %d of %d" % (j, J))
            args = ( M[j], mu0, M0, delta[:,j,:], gam[j] )
            gam[j] = opt_gam_j(args)
            
    logging.debug('Gamma update elapsed time is %0.3f sec for %d samples.' % (time.time() - st, J))
    return gam

def ELBO_0(mu0, M0, gam):
    ## Items in ELBO depends on mu0 and M0
    ## FOr optimization of mu0 and M0

    J = gam.shape[0]
    
    logMu = np.array([EqlogMu(gam[j,:]) for j in xrange(J)])

    log1_Mu = np.array([Eqlog1_Mu(gam[j,:]) for j in xrange(J)])

    EqlogPmu = -J * betaln(mu0*M0, (1-mu0)*M0)
    for j in xrange(J):
        EqlogPmu += (M0*mu0-1)*logMu[j] + (M0*(1-mu0)-1)*log1_Mu[j]

    return EqlogPmu

def neg_ELBO_mu0(mu0, M0, gam):
    return -ELBO_0(mu0, M0, gam)

def opt_mu0(mu0, M0, gam):
    logging.debug("Optimizing mu0")
    #bnds = np.array([[0.01,0.99]])
    bnds = np.array([[0.0,1.0]])
    args=(M0, gam)
    mu0 = opt_par(neg_ELBO_mu0, mu0, args, bnds, 'mu0' )

    return mu0

def neg_ELBO_M0(logM0, mu0, gam):
    return -ELBO_0(mu0, np.exp(logM0), gam)
    
def opt_M0(mu0, M0, gam):
    logging.debug("Optimizing M0")
    
    #bnds = np.array([[-7,7]])    
    bnds = np.array([[-10,10]])    
    
    args = (mu0, gam)
    logM0 = opt_par(neg_ELBO_M0, np.log(M0), args, bnds, 'M0' )
    M0 = np.exp(logM0)
    return M0

def ELBO_M_j(M, delta, gam):
    ## partial ELBO depending on M from each position j 
    ## ELBO used to optimize M 

    if np.ndim(delta) == 1: N = 1
    elif np.ndim(delta) == 2: N= np.shape(delta)[0]

    Mu = EqMu(gam)

    logTheta = np.zeros((N,1))
    log1_Theta = np.zeros((N,1))

    for i in xrange(N):
        logTheta[i] = EqlogTheta(delta[i,:])
        log1_Theta[i] = Eqlog1_Theta(delta[i,:])

    EqlogPtheta = N*EqlogGamma(gam,M)
    for i in xrange(N):
        EqlogPtheta += (M*Mu-1) * logTheta[i] + (M*(1-Mu)-1)*log1_Theta[i]

    return  EqlogPtheta 

def neg_ELBO_M_j(logM, delta, gam):
    return -ELBO_M_j(np.exp(logM), delta, gam)

def opt_M_j(args):

    (M, delta, gam) = args
    #bnds = np.array([[-1, 11]]) # limit delta to [0.0001, 10000]
    # limit delta to [0.0001, 10000], np.log(delta) is [-9.21, 9.21]
    bnds = np.array([[-10, 10]]) 
    M = np.array(M)
    args = (delta, gam)

    logM = opt_par(neg_ELBO_M_j, np.log(M), args, bnds, 'M')
    M = np.exp(logM)

    return M

def opt_M(M, delta, gam, pool = None):
    # M = opt_M(M, delta, gam, pool = pool)
    logging.debug("Optimizing M")

    J= np.shape(M)[0]
    # pdb.set_trace()
    # M = np.array(M)

    if pool is not None:
        args = zip(M, np.transpose(delta, axes =(1,0,2)), gam)
        M = pool.map(opt_M_j, args)

    else:
        for j in xrange(J):
            args = (M[j],delta[:,j,:],gam[j,:])
            M[j] = opt_M_j(args)

    return M

def opt_par(func, x, args, bnds, parlabel):

    # often the fastest method to minimize functions of many 
    # variables uses the Newton-Conjugate Gradient algorithm. 
    # A function which computes the Hessian must be provided. 

    # res = so.minimize(func, x, 
    #     args=args, bounds=bnds, 
    #     method='Newton-CG') 


    # logging.debug("Inside of optimize function. got res")
    # Nelder-Mead is the simplest way to minimize a fairly well-behaved function. 
    # Good for simple minimization problems.
    # Does not use any gradient evaluations, might take longer to find the minimum.
    # if res.success == False: # There is no bounds for Nelder-Mead method
    #     logging.debug(2)  
    #     res = so.minimize(func, x, 
    #         args=args, bounds=bnds, method='Nelder-Mead')
    # pdb.set_trace()

    '''res = so.minimize(func, x, 
        args=args, bounds=bnds, 
        method='L-BFGS-B' ) # limited memory BFGS method'''

    #if res.success == False:
    #    logging.debug(3)  
    res = so.minimize(func, x, 
            args=args, bounds=bnds, method='SLSQP') #Sequential Least SQuares Programming to minimize 
        # a function of several variables with any combination of bounds, equality and inequality constraints

    if res.success == False and parlabel != 'M': 
        logging.debug(1)  
        res = so.minimize(func, x, 
        args=args, bounds=bnds, method='TNC') 
            # truncated Newton algorithm to minimize a function with variables subject to bounds.

    if res.success == False:
        logging.debug(2)
        res = so.minimize(func, x, args=args, bounds=bnds, method='L-BFGS-B' ) # limited memory BFGS method

    # use the gradient of the objective function, which can be given by the user. 
    if res.success == False:
        logging.debug(3)  
        pdb.set_trace()
        res = so.minimize(func, x, bounds=bnds,\
            args=args, method='BFGS') # quasi-Newton method of Broyden, Fletcher, Goldfarb, and Shanno 

    if res.success == False or np.any ( np.isnan(res.x) ) or np.any(np.isinf(res.x)):
        logging.warning("Could not optimize %s or %s is NaN." %(parlabel, parlabel))
        x = np.random.uniform(low=np.amin(bnds), high=np.amax(bnds), size = np.shape(x))      
        return x
        
    return res.x

def beta_mean(p):
    return p[0]*1.0/np.sum(p)

def ELBO_opt(r, n, phi = None, q = None, seed = None, pool = None, vaf = None):

    if pool is not None:
        pool = mp.Pool(processes=pool)
    # t = str(datetime.now)
    f = open('ELBO%s.txt' % str(vaf).replace(".", "_", 1),'w')
    t = time.time()

    # print("ELBO optimization trace starting from %s: \n" %t, file=f)

    if np.ndim(r) == 1: N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: N, J = np.shape(r)
    elif np.ndim(r) == 3: 
        r = np.sum(r, 2) 
        (N, J) = r.shape# sum over non-reference bases
    # r = r.T
    # n = n.T

    if seed is not None: np.random.seed(seed = seed)

    h5file = tempfile.NamedTemporaryFile(suffix='.hdf5')
    logging.info('Storing model updates in %s' % h5file.name)   
    #temp = "tmp.hdf5"
    #logging.info('Storing model updates in %s' % temp)  

    ## Define optimization stopping criterion
    MAXITER = 80
    ELBOTOLPCT = 0.001 *100
    MAXVARITER = 80     
    NORMTOL = 0.1

    ## Initialize model parameters
    if phi is None:
        phi, mu, theta = estimate_mom(r, n)
    else:
        _, mu, theta = estimate_mom(r, n)
    mu0 = phi['mu0']
    M0 = phi['M0']
    M = phi['M']

    ## Initialize the variational parameters
   
    if q is None:
        #delta = np.random.uniform(low = 0.1, high = 100, size = (N,J,2))
        #gam = np.random.uniform(low=0.1, high=100, size = (J,2))
        delta = np.random.uniform(low = 0.0001, high = 10000, size = (N,J,2))
        gam = np.random.uniform(low=0.0001, high=10000, size = (J,2))
        
    else:
        delta = q['delta']
        gam = q['gam']

    phi = {'mu0':mu0, 'M0':M0, 'M':M}
    q = {'delta':delta, 'gam':gam}
    #save_model('initial_value.hdf5', r, n, phi, q)  

    '''# Look at the initial random value of \mu_j
    logging.info("Initial gam: %s" % gam[344,:])
    logging.info("Initial $\mu$: %s" % beta_mean(gam[344,:]))'''
    ## Initialize ELBO    

    elbo = [ELBO(r, n, M, mu0, M0, delta, gam)]
    logging.info("Initial ELBO: %0.2f" % elbo[-1])

    print("M-iteration\tE-iteration\tELBO\tIncrease Percentage\tdelta-deltaprev\tgam-gamprev\tt-gam\tt-delta\tt-mu0\tt-M0\tt-M", file=f)

    print("%d\t%d\t%0.2f\t%0.3f%%\t\t\t\t\t\t\t" %(0, 0, elbo[-1], 0), file=f)

    # print("Initial \tELBO: \t%0.2f" % elbo[-1], file = f)


    ## Optimization
    moditer = 0
    delta_elbo_pct = np.inf

    while moditer < MAXITER and np.abs(delta_elbo_pct) > ELBOTOLPCT:
        # E-step: Update the variational distribution
        variter = 0
        var_elbo = [ elbo[-1] ]
        (norm_delta_delta, norm_delta_gam) = (np.inf, np.inf)
        delta_varelbo_pct = np.inf
        logging.info("E-step")
        while variter < MAXVARITER \
            and delta_varelbo_pct > ELBOTOLPCT \
            and (norm_delta_delta > NORMTOL or norm_delta_gam > NORMTOL):

            #Store the previous parameter values
            (delta_prev, gam_prev) = (np.copy(delta), np.copy(gam))

            #Update the variational distribution
            # pdb.set_trace()
            t0=time.time()
            gam = opt_gam( M, mu0, M0, delta, gam, pool = pool) # mu~Beta(gam)
            t1=time.time()           
            delta = opt_delta(r, n, M, delta, gam, pool = pool)  # theta~Beta(delta)
            t2=time.time()                        

            #Test for convergence
            var_elbo.append(ELBO(r, n, M, mu0, M0, delta, gam))
            delta_varelbo_pct = 100.0*(var_elbo[-1] - var_elbo[-2])/abs(var_elbo[-2])
            logging.info("********Variational Iteration %d of %d********" % (variter+1, MAXVARITER))
            logging.info("ELBO: %0.2f; Percent Change: %0.3f%%" % (var_elbo[-1], delta_varelbo_pct))
            # print("Variational \tELBO: \t%0.2f \tPercent Change: \t%0.3f%%" % (var_elbo[-1], delta_varelbo_pct), file = f)
           

            norm_delta_delta = linalg.norm(delta - delta_prev)
            norm_delta_gam = linalg.norm(gam - gam_prev)
            logging.debug("||delta - delta_prev|| = %0.2f; ||gam - gam_prev|| = %0.2f" 
                % (norm_delta_delta, norm_delta_gam))

            print("%d\t%d\t%0.2f\t%0.3f%%\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t\t\t" %(moditer, variter+1, var_elbo[-1],\
             delta_varelbo_pct,norm_delta_delta, norm_delta_gam, t1-t0,t2-t1), file=f)
            variter += 1

        logging.info("M-step")
        # M-step: Update model parameters
        t0=time.time()
        mu0 = opt_mu0(mu0, M0, gam)
        t1=time.time()
        M0 = opt_M0(mu0, M0, gam)
        t2=time.time()
        M = opt_M(M, delta, gam, pool = pool)
        t3=time.time()

        elbo.append(ELBO(r, n, M, mu0, M0, delta, gam))
        delta_elbo_pct = 100*(elbo[-1] - elbo[-2])/abs(elbo[-2])
        moditer += 1


        # ibic

        # Display results for debugging
        logging.info("----------Iteration %d of %d.----------" % (moditer, MAXITER))
        logging.info("ELBO: %0.2f; Percent Change: %0.3f%%" \
                    % (elbo[-1], delta_elbo_pct))
        
        # print("M-iteration\tE-iteration\tELBO\tIncrease Percentage\tdelta-deltaprev\tgam-gamprev\tt-gam\tt-beta\tt-mu0\tt-M0\tt-M")        
        # print("Iteration %d of %d.\tELBO: \t%0.2f \tPercent Change:\t %0.3f%%" % (moditer, MAXITER, elbo[-1], delta_elbo_pct), file = f)
        print("%d\t%d\t%0.2f\t%0.3f%%\t\t\t\t\t%0.2f\t%0.2f\t%0.2f" %(moditer,0, elbo[-1],delta_elbo_pct, t1-t0,t2-t1,t3-t2), file=f)

        logging.info("M0 = %0.2e" % M0)
        logging.info("mu0 = %0.2f" % mu0)

        '''# Store the model for viewing
        phi = {'mu0':mu0, 'M0':M0, 'M':M}
        q = {'delta':delta, 'gam':gam}
        save_model(h5file.name, r, n, phi, q)'''
        
    print("Total elapsed time is %0.3f seconds." %(time.time()-t), file=f)
    
    f.close()
    return(phi, q)

def estimate_mom(r, n):
    """ Return model parameter estimates using method-of-moments.
    """
    theta = r/(n + np.finfo(np.float).eps) # make sure this is non-truncating division
    if np.ndim(r) == 1: mu = theta
    elif np.ndim(r) > 1: mu = np.mean(theta, 0)
    
    mu0 = np.mean(mu)
    M0 = (mu0*(1-mu0))/(np.var(mu) + np.finfo(np.float).eps) + np.finfo(np.float).eps

    # estimate M. If there is only one replicate, set M as 10 times of M0.
    # If there is multiple replicates, set M according to the moments of beta distribution

    if np.shape(theta)[0] is 1:
        M = 10*M0*np.ones_like(mu)
    else:
        M = (mu*(1-mu))/(np.var(theta, 0) + np.finfo(np.float).eps ) 

    J = len(M)
    for i in xrange(J):
        if M[i] < 1:
            M[i] = 1
   
    phi = {'mu0':mu0, 'M0':M0, 'M':M}
    return phi, mu, theta

def save_model(h5Filename, r, n, phi, q, loc=None, refb=None):

    f = h5py.File(h5Filename, 'w')
    
    f.create_dataset('r', data=r)
    f.create_dataset('n', data=n)

    f.create_group('phi')
    f['phi'].create_dataset('mu0', data=phi['mu0'])
    f['phi'].create_dataset('M0', data=phi['M0'])
    f['phi'].create_dataset('M', data=phi['M'])
    
    f.create_group('q')
    f['q'].create_dataset('delta', data=q['delta'])
    f['q'].create_dataset('gam', data=q['gam'])

    # Save the reference data
    if loc is not None:
        f.create_dataset('loc', data=loc, 
                              chunks=True, fletcher32=True, compression='gzip')
    if refb is not None:
        f.create_dataset('refb', data=refb)

    f.close()

def load_model(h5Filename):
    
    f = h5py.File(h5Filename, 'r')

    out = []
    
    # pdb.set_trace()
    r = f['r'][...]
    out.append(r)

    n = f['n'][...]
    out.append(n)

    phi = {}    
    phi['mu0'] = f['phi/mu0'][...]
    phi['M0'] = f['phi/M0'][...]
    phi['M'] = f['phi/M'][...]
    out.append(phi)

    q = {}
    q['delta'] = f['q/delta'][...]
    q['gam'] = f['q/gam'][...]
    out.append(q)

    if u"loc" in f.keys():
        loc = f['loc'][...]
        out.append(loc)

    if u"refb" in f.keys():
        refb = f['refb'][...]
        out.append(refb)
    
    f.close()
    # pdb.set_trace()

    return tuple(out)

def load_depth(dcFileNameList):
    """ Return (r, n, location, reference base) for a list of depth charts. The
        variable r is the error read depth and n is the total read depth.
    """
    r=[]; n=[]
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    
    loc = []
    refb = {}
    cd = []

    # pdb.set_trace()
    
    for dcFileName in dcFileNameList:
        with open(dcFileName, 'r') as dcFile:
            header = dcFile.readline().strip()
            dc = dcFile.readlines()
            dc = [x.strip().split("\t") for x in dc]
            loc1 = [x[1]+':'+str(x[2]).strip('\000') for x in dc if x[4] in acgt.keys()]
            loc.append( loc1 )            
            refb1 = dict(zip(loc1, [x[4] for x in dc if x[4] in acgt.keys()]))
            refb.update(refb1)
            cd.append( dict(zip(loc1, [map(int, x[5:9]) for x in dc if x[4] in acgt.keys()])) )
            
    loc = list(reduce(set.intersection, map(set, loc)))
    
    def stringSplitByNumbers(x):
        r = re.compile('(\d+)')
        l = r.split(x)
        return [int(y) if y.isdigit() else y for y in l]

    loc = sorted(loc,key = stringSplitByNumbers)
    logging.debug(loc)
    refb = [refb[k] for k in loc]
    
    J = len(loc)
    N = len(dcFileNameList)
    for i in xrange(0, N):
        logging.debug("Processing %s" % dcFileNameList[i])
        c = np.array( [cd[i][k] for k in loc] )
        n1 = np.sum(c, 1)
        #r1 = np.zeros(J)
        refIdx=np.zeros(J)

        for j in xrange(0,J):
            #r1[j] = n1[j] - c[j, acgt[refb[j]]]
            refIdx[j] = 4*j+acgt[refb[j]]
        c = np.delete(c, refIdx, None)
        c = np.reshape(c, (J, 3) )
        #r.append(r1)
        n.append(n1)
        r.append(c)
    r = np.array(r)
    n = np.array(n)
    
    return (r, n, loc, refb)

if __name__ == "__main__":
    main()
