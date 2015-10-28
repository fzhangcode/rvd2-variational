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
from time import time
import time
# <codecell>

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(module)s:%(message)s')

import rvd3

##pool =None
pool = mp.Pool(processes=60)
tocfilename = "../2015-10-15_Plot_time_vs_region_length_rvd3_synthetic_data/synthetic_toc_p1_test_time.txt"
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

            time_all = []

            caseFileList = ["../2015-10-15_Plot_time_vs_region_length_rvd3_synthetic_data/depth_chart/100/top_400_positions/%s"\
                           % filename for filename in toc.Filename[toc.Dilution==dilution] ]
            logging.info('Estimate %s' % caseFileList)   
            
            t0 = time.time()
            (r, n, loc, refb) = rvd3.load_depth(caseFileList)
            t1 = time.time()
            time_load_depth = t1 - t0
            time_all.extend([time_load_depth])

            casephi, caseq, time_ini_model_para, time_ini_var_para, time_ini_ELBO, time_opt_gam, time_opt_delta, time_conv, \
            time_opt_mu0, time_opt_M0, time_opt_M, time_update_ELBO = rvd3.ELBO_opt(r, n, seed = 19860522, pool=60)
            
            time_all.extend([time_ini_model_para, time_ini_var_para, time_ini_ELBO, time_opt_gam, time_opt_delta, time_conv, \
                            time_opt_mu0, time_opt_M0, time_opt_M, time_update_ELBO])

            logging.debug("Saving model in %s" % h5FileName)
            t2 = time.time()
            rvd3.save_model(h5FileName, r, n, casephi, caseq, loc, refb) 
            t3 = time.time()
            time_save = t3 - t2
            time_all.extend([time_save])
            
            print time_all

            items = "load_depth, ini_model_para, ini_var_para, ini_ELBO, opt_gamma, opt_delta, conv, opt_mu0, opt_M0, opt_M, update_ELBO, save_model"
            np.savetxt("time_400.csv", time_all, fmt="%.3f", delimiter=",", header = str(items))

