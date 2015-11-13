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
import time
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(module)s:%(message)s')

rvddir = os.path.join('../../../rvd2/src/python/rvd27')
sys.path.insert(0, rvddir)

##pool =None
pool = mp.Pool(processes=60)
tocfilename = "synthetic_toc_p1_test_time.txt"
toc = pd.read_table(tocfilename)

ngibbs = 4000
nmh = 5
logging.debug("Gibbs step size=%d" % ngibbs)
logging.debug("mh sample size=%d" % nmh)

# Estimate the model for the cases
for dilution in np.unique(toc[toc.isRef=='N'].Dilution):
    logging.debug("Processing dilution: %0.1f" % dilution)
    h5FileName = "Case%s.hdf5" % str(dilution).replace(".", "_", 1)
    try:
        with h5py.File(h5FileName, 'r') as f:
            pass
    except IOError as e:
        dsample = 100
        region = "top_100_positions"
        #caseFileList = ['../2015-10-15_Plot_time_vs_region_length_rvd3_synthetic_data/depth_chart/%s/%s/20100916_c3_p1.07_CGT.dc' % (dsample, region)]
    	#caseFileList = ["../2015-09-28_Run_synthetic_data_set_threshold_delta_gamma_ELBO_update/depth_chart/10000/%s" % filename for filename in toc.Filename[toc.Dilution==dilution]]
        caseFileList = ["../2015-10-15_Plot_time_vs_region_length_rvd3_synthetic_data/depth_chart/100/top_100_positions/%s" % filename for filename in toc.Filename[toc.Dilution==dilution]]
        (r, n, loc, refb) = rvd27.load_depth(caseFileList)
        t = time.time()
        phi, theta_s, mu_s = rvd27.mh_sample(r, n, ngibbs=ngibbs, nmh=nmh, burnin=0.2, pool=pool)
        print("Total elapsed time is %0.3f seconds." %(time.time()-t))
        logging.debug("Saving model in %s" % h5FileName)
        rvd27.save_model(h5FileName, phi, mu=mu_s, theta=theta_s, r=r, n=n,loc=loc, refb=refb)

