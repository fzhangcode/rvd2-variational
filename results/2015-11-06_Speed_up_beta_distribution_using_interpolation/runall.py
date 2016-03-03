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
tocfilename = "./synthetic_toc_p1.txt"
toc = pd.read_table(tocfilename)

'''# Estimate the model for the control
logging.debug("Processing control data.")
h5FileName = "Control.hdf5"
try:
    with h5py.File(h5FileName, 'r') as f:
        pass
except IOError as e:
    controlFileList = ["../2015-09-28_Run_rvd3_synthetic_data_set/depth_chart/10/%s" % filename for filename in toc.Filename[toc.isRef=='Y']]
    (r, n, loc, refb) = rvd3.load_depth(controlFileList)
    controlphi, controlq = rvd3.ELBO_opt(r, n, seed = 1986, pool=60)
    logging.debug("Saving model in %s" % h5FileName)
    rvd3.save_model(h5FileName, r, n, controlphi, controlq, loc, refb)'''


# Estimate the model for the cases
for dilution in np.unique(toc[toc.isRef=='N'].Dilution):
    logging.debug("Processing dilution: %0.1f" % dilution)
    
    h5FileName = "Case%s.hdf5" % str(dilution).replace(".", "_", 1)

    try:
        with h5py.File(h5FileName, 'r') as f:
            pass
    except IOError as e:
        caseFileList = ["../2015-09-28_Run_rvd3_synthetic_data_set/depth_chart/10/%s" % filename for filename in toc.Filename[toc.Dilution==dilution]]
        (r, n, loc, refb) = rvd3.load_depth(caseFileList)
        casephi, caseq = rvd3.ELBO_opt(r, n, seed = 1986, pool=60, vaf = dilution)
        logging.debug("Saving model in %s" % h5FileName)
        rvd3.save_model(h5FileName, r, n, casephi, caseq, loc, refb) 

