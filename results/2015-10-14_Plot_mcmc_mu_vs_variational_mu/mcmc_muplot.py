# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>
# <codecell>
import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import h5py
import pdb
import rvd27

# Insert the src/python/directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)

path = "../2015-10-09_Run_rvd2_with_without_chi2_for_comparison_with_rvd3_synthetic_data/hdf5"
controlFile='%s/100/Control.hdf5' %path
caseFile='%s/100/Case1_0.hdf5' %path
#fname = 'synthetic_partial.hdf5'
position = [344]

def main():
    ## Import the data
    '''if os.path.isfile(fname):
            with h5py.File(fname, 'r') as f:
                        muControl1 = f['muControl'][...]
                        muCase1 = f['muCase'][...]
                        muZ = f['Z'][...]
    else:'''
    with h5py.File(controlFile, 'r') as f:
            muControl = f['mu'][...]
            locControl = f['loc'][...]
    with h5py.File(caseFile, 'r') as f:
            muCase = f['mu'][...]
            locCase = f['loc'][...]

    idx = []
    for pos in position:             
            idx.append(pos)   
            muControl1 = muControl[idx]
            muCase1 = muCase[idx]
            N = 2000

            (muZ,_,_) =rvd27.sample_post_diff(muCase1, muControl1, N) # sample Z

            '''# Save the model
            h5file = h5py.File(fname, 'w')
            if muControl1 is not None:
                        h5file.create_dataset('muControl', data=muControl1,
                                chunks=True, fletcher32=True, compression='gzip')
            if muCase1 is not None:
                        h5file.create_dataset('muCase', data=muCase1,
                                chunks=True, fletcher32=True, compression='gzip')
            if muZ is not None:
                        h5file.create_dataset('Z', data=muZ,
                                chunks=True, fletcher32=True, compression='gzip')
            if position is not None:
                        h5file.create_dataset('loc', data=position,
                                chunks=True, fletcher32=True, compression='gzip')
            h5file.close()'''

    ## plot histogram
    num_bins = 25
    for i in xrange(len(position)):
            fig = plt.figure(figsize=(15,8))
            plt.hist(muCase1[i,:].T, num_bins, facecolor='blue', alpha=0.8)
            plt.hist(muControl1[i,:].T, num_bins, facecolor='green', alpha=0.8)
            plt.legend( ['Case','Control'], frameon=False)
            plt.xlabel('$\mu-\mu_0$', fontsize = 15)
            plt.xticks(rotation=25) 

            fig.suptitle('Position%s' %(position[i]+1), fontsize = 15)
            plt.xticks(rotation=25) 
            plt.savefig('position_%s.png' %(position[i]+1))
            plt.tight_layout()

if __name__ == "__main__":
	main()