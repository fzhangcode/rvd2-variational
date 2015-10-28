# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>
# <codecell>

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import multiprocessing as mp
import logging
import pdb
# <codecell>

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(module)s:%(message)s')

#rvddir = os.path.join('../../src/python/rvd27')
#sys.path.insert(0, rvddir)
import rvd3

##pool =None
pool = mp.Pool(processes=60)
tocfilename = "./synthetic_toc_p1_test_time.txt"
toc = pd.read_table(tocfilename)

# Estimate the model for one case to generate timing table with different regions of interest.
#for f in folderlist:
#    path = 'ELBO/%s' %str(folderlist.index(f))

for dilution in np.unique(toc[toc.isRef=='N'].Dilution):
        logging.debug("Processing dilution: %0.1f" % dilution)
        h5FileName = "Case%s.hdf5" % str(dilution).replace(".", "_", 1)
        try:
            with h5py.File(h5FileName, 'r') as f:
                pass
        except IOError as e:
            dsample = 10000
            region = "top_400_positions"
            #caseFileList = ['./depth_chart/%s/%s/20100916_c3_p1.07_CGT.dc' % (dsample, region)] # For single file
            caseFileList = ["./depth_chart/10000/top_400_positions/%s" % filename for filename in toc.Filename[toc.Dilution==dilution] ]
            logging.info('Estimate %s' % caseFileList)   
            (r, n, loc, refb) = rvd3.load_depth(caseFileList)
            casephi, caseq = rvd3.ELBO_opt(r, n, seed = 19860522, pool=60, dsample=dsample, region=region)
            logging.debug("Saving model in %s" % h5FileName)
            rvd3.save_model(h5FileName, r, n, casephi, caseq, loc, refb) 

