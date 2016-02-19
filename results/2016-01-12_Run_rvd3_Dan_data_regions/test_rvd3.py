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


import rvd3

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():

    controlFile = "./hdf5/E1/c4-1_ADE16_200286.hdf5" #c4-1_MTH1_top400_1986.hdf5"
    caseFile = "./hdf5/E1/c4-8_ADE16_200286.hdf5" #c4-34_MTH1_top400_20160201.hdf5"

    outputFile="./vcf/E1/c4-1_ADE16_200286_c4-8_ADE16_200286" 
	
    #controlFile = "./hdf5/E1/c4-1_MTH1_100.hdf5" #c4-1_MTH1_top400_1986.hdf5"
    #caseFile = "./hdf5/E1/c4-8_MTH1_20160114.hdf5" #c4-34_MTH1_top400_20160201.hdf5"
    ##caseFile = "./hdf5/E1/c4-15_MTH1_20160114.hdf5"
    ##caseFile = "./hdf5/E1/c4-22_MTH1_198605220525.hdf5"
    ##caseFile = "./hdf5/E1/c4-29_MTH1_198605220525.hdf5"
    ##caseFile = "./hdf5/E1/c4-34_MTH1_198605220525.hdf5"
    ##caseFile = "./hdf5/E1/c4-42_MTH1_20160114.hdf5"
    ##caseFile = "./hdf5/E1/c4-47_MTH1_198605220525.hdf5"
    #outputFile="./vcf/E1/c4-1_MTH1_100_c4-8_MTH1_20160114" 
	
    rvd3.test(caseFile, controlFile, alpha=0.05, tau=0, chi2=True, outputFile=outputFile)     

if __name__ == '__main__':
    main()

