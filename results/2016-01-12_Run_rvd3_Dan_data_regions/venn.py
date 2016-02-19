# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import pdb

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3, venn2

def main():

    #plt.figure(figsize=(4,2))
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ############################# E1: plot the minimum intersection of the concomitant mutations in gene MTH1 and ADE16
    ADE16 = 98.886
    MTH1 = 99.603
    inter = MTH1 - (100 - ADE16) 
    region_1 = ADE16 - inter      
    region_2 = MTH1 - inter      
    v = venn2(subsets={'10': region_1, '01': region_2, '11': inter}, set_labels = ('ADE16', 'MTH1'), ax = ax)  
    # set the color, size for the label
    v.get_patch_by_id('10').set_color('r')
    v.get_patch_by_id('01').set_color('b')
    v.get_patch_by_id('10').set_alpha(1.0)
    v.get_patch_by_id('01').set_alpha(1.0)
    '''
    
    ############################# E2: plot the minimum intersection of the concomitant mutations in gene MTH1 and DCD55
    CDC55 = 95.21
    MTH1 = 72.53
    inter = MTH1 - (100 - CDC55) # 67.74
    region_1 = CDC55- inter      # 27.47
    region_2 = MTH1 - inter      # 4.79
    v = venn2(subsets={'10': region_1, '01': region_2, '11': inter}, set_labels = ('CDC55', 'MTH1'), ax = ax)  
    v.get_patch_by_id('10').set_color('blue')
    v.get_patch_by_id('01').set_color('red')
    v.get_patch_by_id('10').set_alpha(0.8)
    v.get_patch_by_id('01').set_alpha(0.8)    
    
    '''   
    '''                                       
    ############################## E3: plot the minimum intersection of the concomitant mutations in gene MUK1 and DAL81
    MUK1 = 82.48
    DAL81 = 56.98
    inter = DAL81 - (100 - MUK1) # 39.46
    region_1 = MUK1- inter      # 43.02
    region_2 = DAL81 - inter      # 17.52
    v = venn2(subsets={'10': region_1, '01': region_2, '11': inter}, set_labels = ('MUK1', 'DAL81'), ax=ax)    # set_labels = ('MUK1(82.48%)', 'DAL81(56.98%)') 
    v.get_patch_by_id('10').set_color('m')
    v.get_patch_by_id('01').set_color('g')
    v.get_patch_by_id('10').set_alpha(0.6)
    v.get_patch_by_id('01').set_alpha(0.7)    
    '''
    
    ##
    v.get_label_by_id('10').set_size(35)
    v.get_label_by_id('01').set_size(35)
    v.get_label_by_id('11').set_size(35) 
    
    v.get_label_by_id('10').set_text('%0.2f%%' %region_1 )
    v.get_label_by_id('01').set_text('%0.2f%%' %region_2 )
    v.get_label_by_id('11').set_text('%0.2f%%' %inter )
        
    #plt.show()
    plt.savefig('venn.png')

if __name__ == "__main__":
    main()