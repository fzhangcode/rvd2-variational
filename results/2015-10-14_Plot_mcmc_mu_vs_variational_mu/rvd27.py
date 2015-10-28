#!/usr/bin/env python

"""rvd27.py: Compute MAP estimates for RVD2.7 model."""

from __future__ import print_function
from __future__ import division

import numpy as np

import scipy.stats as ss
from scipy.special import gammaln

import logging
import multiprocessing as mp
from itertools import repeat
import h5py
import tempfile

import sys
import os
import subprocess
from datetime import date

import re
import pdb
import time

#import pandas as pd
import random

def main():
    import argparse

    # Populate our options, -h/--help is already there.
    description = """
                    RVD is a hierarchical bayesian model for identifying
                    rare variants from short-read sequence data. """
                    
    # create the top-level parser
    argp = argparse.ArgumentParser(prog='rvd', description=description)
    argp.add_argument('--version', action='version', version='%(prog)s 2.7')
    argp.add_argument('-v', '--verbose', dest='verbose', action='count',
                help="increase verbosity (specify multiple times for more)")
    
    # argp.add_argument('cmd', action='store', nargs='*',
    #             choices=['gen', 'gibbs'])
    subparsers = argp.add_subparsers(help='sub-command help')

    ## create subparser to sample the model
    argpGen = subparsers.add_parser('gen', 
        help='Demo: generate simulation sample data from the RVD model')
    argpGen.add_argument('-N', type = int, default=3,
        help='Number of replicates in computer simulation data')
    argpGen.add_argument('-J', type = int, default=10,
        help='Number of positions in computer simulation data')
    argpGen.add_argument('-s', dest='seedint', type = int, default=None,
        help='random process seed.')
    argpGen.set_defaults(func=gen)
    
    ## create subparser for gibbs fitting
    argpGibbs = subparsers.add_parser('gibbs', 
                        help='fit the RVD model using Gibbs sampling')
    argpGibbs.add_argument('dcfile', nargs='+',
                        help='depth chart file name')
    argpGibbs.add_argument('-o', dest='outputFile', 
                default='output',
                help='output HDF5 file name, default (output)')
    argpGibbs.add_argument('-p', '--pool', type=int, default=None,
                help='number of workers in multithread pool, default None')
    argpGibbs.add_argument('-g', '--ngibbs', type=int, default=4000,
                help='gibbs sampling size, default 4000')
    argpGibbs.add_argument('-m', '--nmh', type=int, default=10,
                help='Metropolis-Hastings sampling size, default 10')
    argpGibbs.add_argument('-b', '--burnin', type=float, default=0.2,
                help='burnin, default 0.2')
    argpGibbs.add_argument('-t', '--thin', type=int, default=2,
                help='thin, default 2')
    argpGibbs.add_argument('-s', dest='seedint', type = int, default=None,
        help='random process seed.')
    argpGibbs.set_defaults(func=gibbs_main)

    ## create subparse for one sample one sided test
    argpOneTest = subparsers.add_parser('one_sample_test', 
        help='One side Bayesian posterior density test of one single sample')
    argpOneTest.add_argument('HDF5Name', default=None, help ='HDF5 sample file')
    argpOneTest.add_argument('-i', '--intvl', nargs =2, type=float, default=(0.05, np.inf), 
        help='interval of interest in in posterior distribution.')
    argpOneTest.add_argument('-a','--alpha', type=float, default = 0.05, 
        help='hypothesis test credible level')
    argpOneTest.add_argument('-o', dest='outputFile', 
                default='variants_one_sample',
                help='output HDF5 file name, default (output)')
    argpOneTest.set_defaults(func=one_sample_test_main)

    ## create subparse for germline test
    argpGermTest = subparsers.add_parser('germline_test',
        help='Germline test on a single sample, which includes a one side Bayesian density test and an optional chi square test.')
    argpGermTest.add_argument('HDF5Name', default=None, help ='HDF5 sample file')
    argpGermTest.add_argument('-i', '--intvl', nargs =2, type=float, default=(0.05, np.inf), 
        help='interval of interest in in posterior distribution.')
    argpGermTest.add_argument('-a','--alpha', type=float, default = 0.05, 
        help='hypothesis test credible level')
    argpGermTest.add_argument('-o', dest='outputFile', 
                default='variants_germline',
                help='output HDF5 file name, default (variants_germline)')
    argpGermTest.add_argument('-c','--chi2',
        action='store_false',default=True,
        help='Whether to include chi square test in the germline test, default True (Include)')
    argpGermTest.set_defaults(func=germline_test_main)

    ## creat subparse for two sided paired difference test
    argpDiffTest = subparsers.add_parser('paired_difference_test',
        help='One sided posterior density difference test on control-case paired sample, with an optional chi square test.')
    argpDiffTest.add_argument('controlHDF5Name', default=None, help ='HDF5 control sample file')
    argpDiffTest.add_argument('caseHDF5Name', default=None, help ='HDF5 case sample file')
    argpDiffTest.add_argument('-i', '--intvl', nargs =2, type=float, default=(0.0, np.inf), 
        help='interval of interest in in posterior distribution.')
    argpDiffTest.add_argument('-a','--alpha', type=float, default = 0.05, 
        help='hypothesis test credible level')
    argpDiffTest.add_argument('-o', dest='outputFile', 
                default='variants_paired_difference',
                help='output HDF5 file name, default (variants_paired_difference)')
    argpDiffTest.add_argument('-c','--chi2',
        action='store_false',default=True,
        help='Whether to include chi square test in the paired difference test, default True (Include)')
    argpDiffTest.add_argument('-s',dest='seedint', type = int, default=None,
        help='random process seed.')
    argpDiffTest.add_argument('-n',dest='N', type = int, default=1000,
        help='Posterior difference distribution sampling size.')
    argpDiffTest.set_defaults(func=paired_difference_test_main)

    ## creat subparse for somatic test
    argpSomTest = subparsers.add_parser('somatic_test',
        help='Somatic test, which includes a two sided posterior density difference test and chi square test on the control-case paired sample.')
    argpSomTest.add_argument('controlHDF5Name', default=None, help ='HDF5 control sample file')
    argpSomTest.add_argument('caseHDF5Name', default=None, help ='HDF5 case sample file')
    argpSomTest.add_argument('-i', '--intvl', nargs =2, type=float, default=(0.0, np.inf), 
        help='interval of interest in in posterior distribution.')
    argpSomTest.add_argument('-a','--alpha', type=float, default = 0.05, 
        help='hypothesis test credible level')
    argpSomTest.add_argument('-o', dest='outputFile', 
                default='variants_somatic',
                help='output HDF5 file name, default (variants_somatic)')
    argpSomTest.add_argument('-c','--chi2',
        action='store_false',default=True,
        help='Whether to include chi square test in the germline test, default True (Include)')
    argpSomTest.add_argument('-n',dest='N', type = int, default=1000,
        help='Posterior difference distribution sampling size.')
    argpSomTest.add_argument('-s',dest='seedint', type = int, default=None,
        help='random process seed.')
    argpSomTest.set_defaults(func=somatic_test_main)
           
    ## Parse the arguments (defaults to parsing sys.argv)
    args = argp.parse_args()

    # TODO check what came in on the command line and call optp.error("Useful message") to exit if all is not well

    log_level = logging.WARNING  # default
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    
    # Set up basic configuration, out to stderr with a reasonable format
    logging.basicConfig(level=log_level,
                        format='%(levelname)s:%(module)s:%(message)s')
                        
    # Do actual work here
    args.func(args)
    

def gen(args):
    sample_run(N=args.N, J=args.J, seedint = args.seedint)

def sample_run( N=3, J=10, seedint = None):

    n = np.empty([N,J], dtype=np.int64) # read depth (rep, location)
    n.fill(2000)

    if seedint is not None: 
        np.random.seed(seedint)
    loc = range(J)
    loc = ['chr0:'+str(j+1) for j in loc]
    refb = ['A']*J

    # simulate the control sample
    phi = {'M0':100, 'mu0':0.01, 'M':[1000]*J}

    r, theta, mu= generate_sample(phi, n=n)

    # recontruct the sum of non reference reads r to reads across three non reference reads
    # replicate r in the control sample for uniformity across three non reference reads
    r = np.tile(r, (3,1,1)) # replicate
    r = np.transpose(r, (1,2,0)) # reorder the dimension


    (phi, theta_s, mu_s) = mh_sample(r, n, ngibbs=1000,nmh=5, burnin=0.2, thin=2, pool=None)
    save_model('control_simulation.hdf5', phi, mu=mu_s, theta=theta_s, r=r, n=n, loc=loc,refb=refb)

    # simulate the case sample
    phi = {'M0':100, 'mu0':0.1, 'M':[1000]*J}

    r, theta, mu= generate_sample(phi, n=n)

    # recontruct the sum of non reference reads r to reads across three non reference reads
    # fill the non reference reads for the other two non reference as zero for non-uniformity
    r = [r, np.zeros_like(r), np.zeros_like(r)]
    r = np.transpose(r, (1,2,0)) # reorder the dimension


    (phi, theta_s, mu_s) = mh_sample(r, n, ngibbs=1000,nmh=5, burnin=0.2, thin=2, pool=None)
    save_model('case_simulation.hdf5', phi, mu=mu_s, theta=theta_s, r=r, n=n, loc=loc,refb=refb)


def generate_sample(phi, n):
    """Returns a sample with n reads, N replicates, and
    J locations. The parameters of the model are in the structure phi.
    """
    [N,J] = np.shape(n)
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

def gibbs_main(args):
    """ Top-level function to use gibbs sampling on a set of depth chart files
    """
    gibbs(args.dcfile,ngibbs=args.ngibbs,nmh=args.nmh,burnin=args.burnin, thin=args.thin, pool=args.pool, outputFile=args.outputFile, seedint = args.seedint )

def gibbs(dcfile,ngibbs=4000,nmh=10,burnin=0.2, thin=2, pool=None, outputFile='output', seedint = None ):

    if seedint is not None: 
        np.random.seed(seedint)

    pool = pool
    if pool is not None:
        pool = mp.Pool(processes=pool)

    (r, n, loc, refb) = load_depth(dcfile)
    (phi, theta_s, mu_s) = mh_sample(r, n, ngibbs=ngibbs,nmh=nmh, burnin=burnin, thin=thin, pool=pool)
    save_model(outputFile+'.hdf5', phi, mu=mu_s, theta=theta_s, r=r, n=n, loc=loc,
               refb=refb)

def one_sample_test_main(args):
    intvl=(min(args.intvl),max(args.intvl))
    one_sample_test(HDF5Name=args.HDF5Name, intvl=intvl, alpha=args.alpha, outputFile = args.outputFile)

def one_sample_test(HDF5Name, intvl=[0.05, np.inf], alpha=0.05, chi2 = False, outputFile = 'variants_one_sample'):
    # Bayesian posterior density test of one single sample, one sided; 
    # If 1-alpha percent of samples are within a close intvl in the posterior distribution of a position, 
    # the position will be reported in a vcf file and hdf5 file

    logging.debug('Running one-sample-test on posterior distribution of sample %s' %os.path.basename(HDF5Name))
    (Phi, Theta, Mu, loc, R, N, refb) = load_model(HDF5Name)

    if len(np.shape(intvl))==1:
        intvl = [intvl]
    postP = bayes_test(Mu, intvl,'close')
    bayescall = postP > 1-alpha

    if chi2:
        chi2call , chi2P = chi2combinetest(R, N, bayescall)
        # combine the calls from bayesian posterior test and chi2 test
        call = np.logical_and(bayescall, chi2call)
    else:
        call = bayescall

    if outputFile is not None:
        # output vcf file
        vcfFilename = outputFile+'.vcf'

        altb = np.copy(refb)
        altb = ['.' for b in altb]

        write_vcf(vcfFilename, loc, bayescall, refb, altb, Mu, R, N)

        # output hdf5 file
        h5Filename = outputFile +'.hdf5'
        h5file = h5py.File(h5Filename, 'w')

        h5file.create_dataset('call', data=bayescall)
        h5file.create_dataset('refb', data=refb)
        h5file.create_dataset('loc', data=loc, 
                              chunks=True, fletcher32=True, compression='gzip')        
        h5file.create_dataset('Mu', data=Mu, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('N', data=N, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('R', data=R, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.close()
    return call, loc, Mu

def germline_test_main(args):
    intvl=(min(args.intvl),max(args.intvl))
    germline_test(HDF5Name=args.HDF5Name, intvl=intvl, alpha=args.alpha, chi2=args.chi2, outputFile=args.outputFile )

def germline_test(HDF5Name, intvl=[0.05, np.inf], alpha=0.05, chi2=True, outputFile='variants_germline'):
    # Germline test: Bayesian posterior density test of one single sample, one sided together with optional chi square test; 
    [call, loc, Mu] = one_sample_test(HDF5Name, intvl, alpha, chi2, outputFile)
    return call, loc, Mu

def paired_difference_test_main(args):
    intvl=(min(args.intvl),max(args.intvl))
    paired_difference_test(controlHDF5Name=args.controlHDF5Name, caseHDF5Name=args.caseHDF5Name, 
        N =args.N, intvl=intvl, alpha=args.alpha, chi2=args.chi2, outputFile=args.outputFile, seedint = args.seedint )

def paired_difference_test(controlHDF5Name, caseHDF5Name, N = 1000, intvl=[0, np.inf], alpha=0.05,
 chi2 = False, outputFile = 'variants_paired_difference', seedint = None):
    # Bayesian posterior difference test of case-control paired samples, one sided; 
    # If 1-alpha percent of samples are within a open intvl in the posterior distribution of a position, 
    # the position will be reported in a vcf file and hdf5 file

    if seedint is not None: 
        np.random.seed(seedint)

    logging.debug('Running Bayesian posterior difference test on paired sample %s and %s)' \
    %(os.path.basename(caseHDF5Name), os.path.basename(controlHDF5Name)))
    
    loc, refb, altb, controlMu, controlMu0, controlR, controlN, \
             caseMu, caseMu0, caseR, caseN = load_dualmodel(controlHDF5Name, caseHDF5Name)

    (Z, _, _) = sample_post_diff(caseMu-caseMu0, controlMu-controlMu0, N)

    if len(np.shape(intvl))==1:
        intvl = [intvl]
    postP = bayes_test(Z, intvl,'open')

    bayescall = postP > 1-alpha

    if chi2:
        chi2call , chi2P = chi2combinetest(caseR,caseN,bayescall)
        call = np.logical_and(bayescall, chi2call)
    else:
        call = bayescall

    if outputFile is not None:
        # output vcf file
        vcfFilename = outputFile+'.vcf'

        write_dualvcf(vcfFilename, loc, bayescall, refb, altb, controlMu, controlR, controlN,
                  caseMu, caseR, caseN)
        # output hdf5 file
        h5Filename = outputFile +'.hdf5'
        h5file = h5py.File(h5Filename, 'w')

        h5file.create_dataset('call', data=bayescall)
        h5file.create_dataset('refb', data=refb)
        h5file.create_dataset('loc', data=loc, 
                              chunks=True, fletcher32=True, compression='gzip')        
        h5file.create_dataset('controlMu', data=controlMu, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('caseMu', data=caseMu, 
                              chunks=True, fletcher32=True, compression='gzip')        
        h5file.create_dataset('controlN', data=controlN, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('caseN', data=caseN, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('controlR', data=controlR, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('caseR', data=caseR, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.close()

    return  call, loc, refb, altb, controlMu, controlMu0, controlR, controlN, \
             caseMu, caseMu0, caseR, caseN

def somatic_test_main(args):
    intvl=(min(args.intvl),max(args.intvl))
    somatic_test(controlHDF5Name=args.controlHDF5Name, caseHDF5Name=args.caseHDF5Name, 
        N =args.N, intvl=intvl, alpha=args.alpha, chi2=args.chi2, outputFile=args.outputFile, seedint = args.seedint )

def somatic_test(controlHDF5Name, caseHDF5Name, N = 1000, intvl=[0, np.inf], 
    alpha=0.05, chi2 = True, outputFile='variants_somatic', seedint = None):
    logging.debug('Running somatic test on paired sample %s and %s)' 
    %(os.path.basename(caseHDF5Name), os.path.basename(controlHDF5Name)))
    # Somatic test: Bayesian posterior difference test of case-control paired samples, one sided, together with optional chi square test; 
    bayescall_pos, loc, refb, altb, controlMu, controlMu0, controlR, controlN, \
             caseMu, caseMu0, caseR, caseN = paired_difference_test(controlHDF5Name, caseHDF5Name, N, intvl, alpha, chi2, outputFile = None, seedint = seedint)
    bayescall_neg, _, _, _, _, _, _, _, \
             _, _, _, _ = paired_difference_test(caseHDF5Name, controlHDF5Name, N, intvl, alpha, chi2,  outputFile = None, seedint = seedint)

    # Combine the calls from two sided paired difference test
    call = [bayescall_pos,bayescall_neg]
    call = np.sum(call,0) > 0

    # output result
    if outputFile is not None:
        # output vcf file
        vcfFilename = outputFile+'.vcf'

        write_dualvcf(vcfFilename, loc, call, refb, altb, controlMu, controlR, controlN,
                  caseMu, caseR, caseN)
        # output hdf5 file
        h5Filename = outputFile +'.hdf5'
        h5file = h5py.File(h5Filename, 'w')

        h5file.create_dataset('call', data=call)
        h5file.create_dataset('refb', data=refb)
        h5file.create_dataset('loc', data=loc, 
                              chunks=True, fletcher32=True, compression='gzip')        
        h5file.create_dataset('controlMu', data=controlMu, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('caseMu', data=caseMu, 
                              chunks=True, fletcher32=True, compression='gzip')        
        h5file.create_dataset('controlN', data=controlN, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('caseN', data=caseN, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('controlR', data=controlR, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.create_dataset('caseR', data=caseR, 
                              chunks=True, fletcher32=True, compression='gzip')
        h5file.close()

    return call, loc, controlMu, caseMu


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

def load_model(h5Filename):
    """ Returns the RVD2.7 model samples and parameters.
    Takes an hdf5 filename and returns phi and other parameters
    """

    out = []

    with h5py.File(h5Filename, 'r') as h5file:
        # Load phi - it always exists
        phi = {'mu0': h5file['phi/mu0'][()],
               'M0': h5file['phi/M0'][()],
               'M': h5file['phi/M'][...]}
        out.append(phi)
        
        # Load theta if it exists
        if u"theta" in h5file.keys():
            theta = h5file['theta'][...]
            out.append(theta)
            
        # Load mu if it exists
        if u"mu" in h5file.keys():
            mu = h5file['mu'][...]
            out.append(mu)
            
        # Load loc if it exists
        if u"loc" in h5file.keys():
            loc = h5file['loc'][...]
            out.append(loc)
    
        # Load r if it exists
        if u"r" in h5file.keys():
            r = h5file['r'][...]
            out.append(r)

        if u"n" in h5file.keys():
            n = h5file['n'][...]
            out.append(n)

        if u"refb" in h5file.keys():
            refb = h5file['refb'][...]
            out.append(refb)

        h5file.close()

            
    return tuple(out)
    
def save_model(h5Filename, phi, mu=None, theta=None, r=None, n=None, loc=None, refb=None):
    """ Save the RVD2.7 model samples and parameters """
    
    # TODO add attributes to hdf5 file
    h5file = h5py.File(h5Filename, 'w')
    
    # Save the model parameters (phi)
    h5file.create_group('phi')
    h5file['phi'].create_dataset('mu0', data=phi['mu0'])
    h5file['phi'].create_dataset('M0', data=phi['M0'])
    h5file['phi'].create_dataset('M', data=phi['M'], 
                                      chunks=True, 
                                      fletcher32=True, 
                                      compression='gzip')
    
    # Save the latent variables if available.
    if mu is not None:
        h5file.create_dataset('mu', data=mu, 
                              chunks=True, fletcher32=True, compression='gzip')
    if theta is not None:
        h5file.create_dataset('theta', data=theta, 
                              chunks=True, fletcher32=True, compression='gzip')
    
    # Save the data used for fitting the model if available
    if r is not None:
        h5file.create_dataset('r', data=r, 
                              chunks=True, fletcher32=True, compression='gzip')
    if n is not None:
        h5file.create_dataset('n', data=n, 
                              chunks=True, fletcher32=True, compression='gzip')

    # Save the reference data
    if loc is not None:
        h5file.create_dataset('loc', data=loc, 
                              chunks=True, fletcher32=True, compression='gzip')
    if refb is not None:
        h5file.create_dataset('refb', data=refb)

    h5file.close()

def load_dualmodel(controlHDF5Name, caseHDF5Name):
    '''
        Load and synchonize the Case and Control Model files 
    '''
    # Load the Case and Control Model files
    (controlPhi, controlTheta, controlMu, controlloc, controlR, controlN, controlrefb) = load_model(controlHDF5Name)
    (casePhi, caseTheta, caseMu, caseloc, caseR, caseN,_) = load_model(caseHDF5Name)

    # Extract the common locations in case and control
    caselocIdx = [i for i in xrange(len(caseloc)) if caseloc[i] in controlloc]
    controllocIdx = [i for i in xrange(len(controlloc)) if controlloc[i] in caseloc]

    caseMu = caseMu[caselocIdx,:]
    controlMu = controlMu[controllocIdx,:]
    caseR = caseR[:,caselocIdx,:]
    controlR = controlR[:,controllocIdx,:]

    loc = caseloc[caselocIdx]

    J = len(loc)
    
    refb = controlrefb[controllocIdx]

    altb = []
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    for i in xrange(J):
        r = np.squeeze(caseR[:,i,:]) # replicates x bases
        
        # Make a list of the alternate bases for each replicate
        acgt_r = ['A','C','G','T']
        del acgt_r[ acgt[refb[i]] ]

        if not np.iterable(np.argmax(r, axis=-1)):
            altb_r = acgt_r[np.argmax(r, axis=-1)]
        else:
            altb_r = [acgt_r[x] for x in np.argmax(r, axis=-1)]
        altb.append(altb_r[0])

    return loc, refb, altb, controlMu, controlPhi['mu0'], controlR, controlN,\
           caseMu, casePhi['mu0'], caseR, caseN

def write_dualvcf(outputFile, loc, call, refb, altb, controlMu=None, controlR=None, controlN=None,\
                  caseMu=None, caseR=None, caseN=None):
    controlMu = np.mean(controlMu, axis=1);
    caseMu = np.mean(caseMu, axis=1)

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

    print('##Posterior test in cancer-normal-paired sample.', file=vcfF)

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

            print ("chr%s\t%d\t.\t%s\t%s\t.\tPASS\tCOAF=%0.3f;CAAF=%0.3f\tAU:CU:GU:TU\t%d:%d:%d:%d\t%d:%d:%d:%d" \
                   % (chrom[i], pos[i], refb[i], altb[i],controlMu[i]*100.0, caseMu[i]*100.0,\
                      controlR4[0], controlR4[1], controlR4[2], controlR4[3],\
                      caseR4[0], caseR4[1], caseR4[2], caseR4[3]), file=vcfF)
                
    vcfF.close()

def write_vcf(outputFile, loc, call, refb, altb, Mu, R, N):
    '''
        Write high confidence variant calls from one sample test to VCF 4.2 file.
    '''
    Mu = np.mean(Mu, axis=1);

    R = np.median(R,0)


    J = len(loc)
    
    today=date.today()
    
    chrom = [x.split(':')[0][3:] for x in loc]
    pos = [int(x.split(':')[1]) for x in loc]
    
    vcfF = open(outputFile,'w')
    
    print("##fileformat=VCFv4.1", file=vcfF)
    print("##fileDate=%0.4d%0.2d%0.2d" % (today.year, today.month, today.day), file=vcfF)

    print("##source=rvd2", file=vcfF)

    print('##Posterior test in a single sample.', file=vcfF)

    # print("##Posterior difference threshold = %0.2f" %tau, file=vcfF)
    # print("##Probability threshold alpha = %0.2f" %alpha, file=vcfF)
    
    # print("##Chi square test is included", file=vcfF)

    uniquechrom = set(chrom)
    uniquechrom = list(uniquechrom)

    for i in xrange(len(uniquechrom)):
        seq = [idx for idx, name in enumerate(chrom) if name==uniquechrom[i]]
        seqlen = len(seq)
        print("##contig=<ID=%(chr)s,length=%(seqlen)d>" %{'chr': uniquechrom[i],'seqlen': seqlen}, file=vcfF)

    
    print("##INFO=<ID=COAF,Number=1,Type=Float,Description=\"Allele Frequency\">", file=vcfF)

    print("##FORMAT=<ID=AU,Number=1,Type=Integer,Description=\"Number of 'A' alleles\">", file=vcfF)
    print("##FORMAT=<ID=CU,Number=1,Type=Integer,Description=\"Number of 'C' alleles\">", file=vcfF)
    print("##FORMAT=<ID=GU,Number=1,Type=Integer,Description=\"Number of 'G' alleles\">", file=vcfF)
    print("##FORMAT=<ID=TU,Number=1,Type=Integer,Description=\"Number of 'T' alleles\">", file=vcfF)
    
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample", file=vcfF)

    for i in xrange(J):
        if call[i]:           
            # restore R
            actg = ['A','C','G','T']

            idx = actg.index(refb[i])
            R4 = np.zeros(4)
            R4[idx] = np.median(N[:,i])-np.sum(R[i,:])
            
            for d in xrange(idx):
                R4[d] = R[i,d]
            for d in xrange(3-idx):
                R4[d+idx+1] = R[i,d+idx]

            print ("chr%s\t%d\t.\t%s\t%s\t.\tPASS\tAF=%0.3f\tAU:CU:GU:TU\t%d:%d:%d:%d" \
                   % (chrom[i], pos[i], refb[i], altb[i],Mu[i]*100.0, R4[0], R4[1], R4[2], R4[3]), file=vcfF)
    vcfF.close()

def complete_ll(phi, r, n, theta, mu):
    """ Return the complete data log-likelihood.
    """
    alpha0 = phi['M0']*phi['mu0'] + np.finfo(np.float).eps
    beta0 = phi['M0']*(1-phi['mu0']) + np.finfo(np.float).eps
    
    alpha = phi['M']*mu + np.finfo(np.float).eps
    beta = phi['M']*(1 - mu) + np.finfo(np.float).eps
    
    # Bound theta away from 0 or 1
    theta[theta < np.finfo(np.float).eps] = np.finfo(np.float).eps
    theta[theta > 1-np.finfo(np.float).eps] = 1 - np.finfo(np.float).eps
    
    logPmu = beta_log_pdf(mu, alpha0, beta0)
    logPtheta = beta_log_pdf(theta, alpha, beta)
    logPr = ss.binom.logpmf(r, n, theta)
    
    return np.sum(logPmu + logPtheta + logPr)

def estimate_mom(r, n):
    """ Return model parameter estimates using method-of-moments.
    """

    theta = r/(n + np.finfo(np.float).eps) # make sure this is non-truncating division
    if np.ndim(r) == 1: mu = theta
    elif np.ndim(r) > 1: mu = np.mean(theta, 0)
    
    mu0 = np.mean(mu)
    M0 = (mu0*(1-mu0))/(np.var(mu) + np.finfo(np.float).eps)

    # estimate M. If there is only one replicate, set M as 10 times of M0.
    # If there is multiple replicates, set M according to the moments of beta distribution
    if np.shape(theta)[0] is 1:
        M = 10*M0*np.ones_like(mu)
    else:
        M = (mu*(1-mu))/(np.var(theta, 0) + np.finfo(np.float).eps )    

    phi = {'mu0':mu0, 'M0':M0, 'M':M}
    return phi, mu, theta
    
def samplelocMuMH(args):
    # Sample from the proposal distribution at a particular location
    mu, Qsd, theta, M, alpha0, beta0 = args
    
    # TODO put an escape in here to avoid an infinite loop
    while True:
        mu_p = ss.norm.rvs(mu, Qsd)
        if (0 < mu_p < 1): break
        
    # Log-likelihood for the proposal mu
    alpha_p = mu_p*M + np.finfo(np.float).eps
    beta_p = (1-mu_p)*M + np.finfo(np.float).eps
    logPmu_p = beta_log_pdf(mu_p, alpha0, beta0) \
                + np.sum(beta_log_pdf(theta, alpha_p, beta_p))
                    
    # Log-likelihood for the current mu
    alpha = mu*M + np.finfo(np.float).eps
    beta = (1-mu)*M + np.finfo(np.float).eps
    logPmu = beta_log_pdf(mu, alpha0, beta0) \
                + np.sum(beta_log_pdf(theta, alpha, beta))
                            
    # Accept new mu if it increases posterior pdf or by probability
    loga = logPmu_p - logPmu
    if (loga > 0 or np.log(np.random.random()) < loga): 
        mu = mu_p
    return mu
    
def sampleMuMH(theta, mu0, M0, M, mu=ss.beta.rvs(1, 1), Qsd=0.1, burnin=0, nmh=1, thin=0, pool=None):
    """ Return a sample of mu with parameters mu0 and M0.
    """
    if np.ndim(theta) == 1: (N, J) = (1, np.shape(theta)[0])
    elif np.ndim(theta) > 1: (N, J) = np.shape(theta)
        
    alpha0 = mu0*M0 + np.finfo(np.float).eps
    beta0 = (1-mu0)*M0 + np.finfo(np.float).eps

    mu_s = np.zeros( (nmh, J) ) 
    for ns in xrange(0, nmh):
        if pool is not None:
            args = zip(mu, Qsd, theta.T, M, repeat(alpha0, J), repeat(beta0, J))
            mu = pool.map(samplelocMuMH, args)
        else:
            for j in xrange(0, J):
                args = (mu[j], Qsd[j], theta[:,j], M[j], alpha0, beta0)
                mu[j] = samplelocMuMH(args)
        # Save the new sample
        mu_s[ns, :] = np.copy(mu)

    if burnin > 0.0:
        mu_s = np.delete(mu_s, np.s_[0:np.int(burnin*nmh):], 0)
    if thin > 0:
        mu_s = np.delete(mu_s, np.s_[::thin], 0)
    
    return mu_s
    
def mh_sample(r, n, ngibbs=4000,nmh=10, burnin=0.2, thin=2, pool=None):
    """ Return MAP parameter and latent variable estimates obtained by 

    Metropolis-Hastings sampling.
    By default, sample 4000 M-H with a 20% burn-in and thinning factor of 2. 
    Stop when the change in complete data log-likelihood is less than 0.01%.
    """

    if np.ndim(r) == 1: N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: N, J = np.shape(r)
    elif np.ndim(r) == 3: 
        r = np.sum(r, 2) # sum over non-reference bases
        N, J = np.shape(r)
    
    # Initialize a hdf5 file for logging model progress
    h5Filename = tempfile.NamedTemporaryFile(suffix='.hdf5').name
    logging.debug("Storing temp data in %s" % h5Filename)
    h5file = h5py.File(h5Filename, 'w')
    h5file.create_group('phi')
    h5file['phi'].create_dataset('mu0', (1,), dtype='f')
    h5file['phi'].create_dataset('M0', (1,), dtype='f')
    h5file['phi'].create_dataset('M', (J,), dtype='f')
    h5file.create_dataset('theta_s', (N, J, ngibbs), dtype='f')
    h5file.create_dataset('mu_s', (J, ngibbs), dtype='f')
    # Initialize estimates using MoM
    phi, mu, theta = estimate_mom(r, n)
    logging.debug("MoM: mu0 = %0.3e; M0 = %0.3e." % (phi['mu0'], phi['M0']) )

    # Correct MoM estimates to be non-trivial
    mu[mu < np.finfo(np.float).eps*1e4] = phi['mu0']
    theta[theta < np.finfo(np.float).eps*1e4] = phi['mu0']
    phi['M'][phi['M'] < np.finfo(np.float).eps *1e4] = 1

    # Set Qsd for proposal distribution according to mom of mu
    
    def boundfn(mu):       
        bound = 0.001
        if bound < mu < 1-bound:
            return mu*(1-mu)/10
        else:
            return bound/10
        

##    def boundfn(mu):       
##        bound = 0.001
##        if bound < mu < 1-bound:
##            return mu/10
##        else:
##            return bound/10
        
    Qsd = map(boundfn, mu)       

    # Sample theta, mu by gibbs sampling
    theta_s = np.zeros( (N, J, ngibbs) )
    mu_s = np.zeros( (J, ngibbs) )
    for i in xrange(ngibbs):
        if i % 100 == 0 and i > 0: logging.debug("Gibbs Iteration %d" % i)
            
        # Draw samples from p(theta | r, mu, M) by Gibbs
        alpha = r + mu*phi['M'] + np.finfo(np.float).eps
        beta = (n - r) + (1-mu)*phi['M'] + np.finfo(np.float).eps
        theta = ss.beta.rvs(alpha, beta)
        
        # Draw samples from p(mu | theta, mu0, M0) by Metropolis-Hastings
        mu_mh = sampleMuMH(theta, phi['mu0'], phi['M0'], phi['M'], mu=mu, Qsd=Qsd, nmh=nmh, pool=pool)
        mu = np.median(mu_mh, axis=0)
        # Store the sample
        theta_s[:,:,i] = np.copy(theta)
        mu_s[:,i] = np.copy(mu)
        # Update parameter estimates
        # phi['mu0'] = np.mean(mu)
        # phi['M0'] = (phi['mu0']*(1-phi['mu0']))/(np.var(mu) + np.finfo(np.float).eps)
        # TODO update for M
        
        # Store the current model
        h5file['phi']['mu0'][0] = phi['mu0']
        h5file['phi']['M0'][0] = phi['M0']
        h5file['phi']['M'][...] = phi['M']
        h5file['theta_s'][:,:,i] = theta
        h5file['mu_s'][:,i] = mu
        h5file.flush()
    # Apply the burn-in and thinning
    if burnin > 0.0:
        mu_s = np.delete(mu_s, np.s_[0:np.int(burnin*ngibbs):], 1)
        theta_s = np.delete(theta_s, np.s_[0:np.int(burnin*ngibbs):], 2)
    if thin > 1:
        mu_s = np.delete(mu_s, np.s_[::thin], 1)
        theta_s = np.delete(theta_s, np.s_[::thin], 2)
    
    h5file.close()
    return (phi, theta_s, mu_s)

def beta_log_pdf(x, a, b):
    return gammaln(a+b) - gammaln(a) - gammaln(b) \
            + (a-1)*np.log(x+np.finfo(np.float).eps) \
            + (b-1)*np.log(1-x+np.finfo(np.float).eps)


def ll(phi, r):
    """ Return the log-likelihood of the data, r, under the model, phi.
    """
    pass

def make_pileup(bamFileName, fastaFileName, region):
    """ Creates a pileup file using samtools mpileup in /pileup directory.
    """
    
    # Check that BAM file exists
    assert os.path.isfile(bamFileName), "BAM file does not exist: %s" % bamFileName
    
    # Check that FASTA reference file exists
    assert os.path.isfile(fastaFileName), "FASTA file does not exist: %s" % fastaFileName
    
    # Create pileup directory if it doesn't exist
    if not os.path.isdir("pileup"):
        os.makedirs("pileup")
    
    # Format the samtools call
    callString = ["samtools", "mpileup", "-d", "1000000", "-r", "%s" % region,
                  "-f", "%s" % fastaFileName, "%s" % bamFileName]
                  
    # Remove the extension from the bam filename and replace with .pileup
    pileupFileName = bamFileName.split("/")[-1]
    pileupFileName = os.path.join("pileup", 
                                  "%s.pileup" % pileupFileName.split(".", 1)[0])
    
    # Run samtools pileup only if the file doesn't already exist.
    #try:
    #    with open(pileupFileName, 'r'):
    #        logging.debug("Pileup file exists: %s" % pileupFileName)
    #except IOError:
    logging.debug("[call] %s", " ".join(callString))
    with open(pileupFileName, 'w') as fout:
        subprocess.call(callString, stdout=fout)
    return pileupFileName

def make_depth(pileupFileName):
    """ Generates a depth chart file for each pileup file and stores it in the
        /depth_chart directory. The folder will be created if it doesn't exist.
    """
    
    if not os.path.isdir("depth_chart"):
        os.makedirs("depth_chart")
    
    # TODO replace this with a python version.
    callString = ["../../bin/pileup2dc", "%s" % pileupFileName]

    dcFileName = pileupFileName.split("/")[-1]
    dcFileName = os.path.join("depth_chart", 
                                  "%s.dc" % dcFileName.split(".", 1)[0])
    #try:
    #    with open(dcFileName, 'r'):
#		logging.debug("Depth chart file exists: %s" % dcFileName) 
    #except IOError:
    logging.debug("Converting %s to depth chart." % pileupFileName)
    with open(dcFileName, 'w') as fout:
        subprocess.call(callString, stdout=fout)
    return dcFileName

def load_depth(dcFileNameList):
    """ Return (r, n, location, reference base) for a list of depth charts. The
        variable r is the error read depth and n is the total read depth.
    """
    r=[]; n=[]
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    
    loc = []
    refb = {}
    cd = []

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
 
def sample_post_diff(muCaseG, muControlG, N):
    """ Return N samples from the posterior distribution for 
         u_j|r_case - u_j|r_control. """
    
    nCase = muCaseG.shape[1]
    nControl = muControlG.shape[1]
    
    caseSample = np.random.choice(nCase, size=N, replace=True)
    controlSample = np.random.choice(nControl, size=N, replace=True)
    
    muCaseS = muCaseG[:, caseSample]
    muControlS = muControlG[:, controlSample]

    Z = muCaseS - muControlS

    return (Z, muCaseS, muControlS)
 
def bayes_test(Z, roi, type = 'close'):
    """ Return posterior probabilities in regions defined in list of tuples (roi)
        from samples in columns of Z. """
    (J,N)=np.shape(Z)
    
    nTest = len(roi) # get the number of regions to compute probabilities 

    p = np.zeros((J,nTest))
    for i in xrange(nTest):
        for j in xrange(J):
                if type == 'close':
                    # somatic test
                    p[j,i] = np.float( np.sum( np.logical_and( (Z[j,:] >= roi[i][0]), (Z[j,:] <= roi[i][1]) ) ) ) / N
                elif type == 'open':
                    # diff test
                    p[j,i] = np.float( np.sum( np.logical_and( (Z[j,:] > roi[i][0]), (Z[j,:] < roi[i][1]) ) ) ) / N
    
    p = np.sum(p,1) # add up all the regions in each position.
    
    return p # combine probabilities from all regions.
    
if __name__ == '__main__':
    main()
    
