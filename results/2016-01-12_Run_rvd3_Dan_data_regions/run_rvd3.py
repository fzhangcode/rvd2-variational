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

rvddir = os.path.join('./')
sys.path.insert(0, rvddir)
import rvd3

##pool =None
pool = mp.Pool(processes=62)

# Estimate the model for the control
logging.debug("Processing control data.")
h5FileName = "c4-8_ADE16_200286.hdf5"
try:
    with h5py.File(h5FileName, 'r') as f:
        pass
except IOError as e:
    filename = "c4-8_ADE16_200286.dc"
    controlFileList = ["./dc/E1/%s" % filename]
    (r, n, loc, refb) = rvd3.load_depth(controlFileList)
    controlphi, controlq = rvd3.ELBO_opt(r, n, seed = 20160210, pool=62)
    logging.debug("Saving model in %s" % h5FileName)
    rvd3.save_model(h5FileName, r, n, controlphi, controlq, loc, refb)
