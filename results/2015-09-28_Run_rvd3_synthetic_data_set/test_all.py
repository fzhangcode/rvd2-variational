#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss

# Insert the src/python/rvd27 directory at front of the path
# rvddir = os.path.join('../../../src/python/rvd27')
# sys.path.insert(0, rvddir)
import rvd3

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    
    dilutionList = (0.1,0.3,1.0,10.0,100.0)

    folderList = ('hdf5/10',\
                  'hdf5/100',\
                  'hdf5/1000',\
                  'hdf5/10000') 

    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
    	for f in folderList:
            path='vcf/%s' % str(10**(folderList.index(f)+1))

            if not os.path.exists(path):
                os.makedirs(path)
                
            controlFile = "%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
            
            outputFile='%(path)s/vcf%(dilution)s' %{'path':path,'dilution':str(d).replace('.','_')}
            
            rvd3.test(caseFile, controlFile, alpha=0.05, tau=0, chi2=False, outputFile=outputFile)     

if __name__ == '__main__':
    main()

