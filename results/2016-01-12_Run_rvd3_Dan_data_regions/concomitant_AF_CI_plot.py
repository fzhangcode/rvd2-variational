# -*- coding: utf-8 -*-

import numpy as np
import pdb
import matplotlib.pyplot as plt
import h5py
import rvd3

def main():
    # mutant positions
    loc = [339, 44] #[1014740, 200286]
    
    # time point, generations
    ind = np.array([7, 70, 133, 196, 266, 322, 385, 448])
    N = len(ind)
    
    # record the allele frequency of one position in all the time points
    AF_all = np.zeros(N)
    err_all = np.zeros(N)

    # plot the allele frequency (mu) and the credible interval of mu
    fig = plt.figure()
    ax = fig.add_subplot(111)  

    #################################### E1: MTH1 & ADE16 ################
    print '=======MTH1======='
    AF_all[0], err_all[0] = read('./hdf5/E1/c4-1_MTH1_100.hdf5', loc[0])  
    AF_all[1], err_all[1] = read('./hdf5/E1/c4-8_MTH1_20160114.hdf5', loc[0])  
    AF_all[2], err_all[2] = read('./hdf5/E1/c4-15_MTH1_20160114.hdf5', loc[0])  
    AF_all[3], err_all[3] = read('./hdf5/E1/c4-22_MTH1_198605220525.hdf5', loc[0])  
    AF_all[4], err_all[4] = read('./hdf5/E1/c4-29_MTH1_198605220525.hdf5', loc[0])  
    AF_all[5], err_all[5] = read('./hdf5/E1/c4-34_MTH1_198605220525.hdf5', loc[0])  
    AF_all[6], err_all[6] = read('./hdf5/E1/c4-42_MTH1_20160114.hdf5', loc[0])  
    AF_all[7], err_all[7] = read('./hdf5/E1/c4-47_MTH1_198605220525.hdf5', loc[0]) 	
    
    #rects1 = ax.bar(ind, AF_all, width, color='r', yerr=err_all, ecolor='k')
    plt.plot(ind, AF_all, linestyle='-', linewidth=12, marker='o', markersize=16, color='b', label='MTH1')
    plt.errorbar(ind, AF_all, yerr=err_all, linestyle='None', linewidth=4, marker = 'None', elinewidth=3, capsize=18, color='b')
    

    print '=======ADE16======='
    AF_all[0], err_all[0] = read('./hdf5/E1/c4-1_ADE16_200286.hdf5', loc[1])  
    AF_all[1], err_all[1] = read('./hdf5/E1/c4-8_ADE16_200286.hdf5', loc[1])  
    AF_all[2], err_all[2] = read('./hdf5/E1/c4-15_ADE16_200286.hdf5', loc[1])  
    AF_all[3], err_all[3] = read('./hdf5/E1/c4-22_ADE16_200286.hdf5', loc[1])  
    AF_all[4], err_all[4] = read('./hdf5/E1/c4-29_ADE16_200286.hdf5', loc[1])  
    AF_all[5], err_all[5] = read('./hdf5/E1/c4-34_ADE16_200286.hdf5', loc[1])  
    AF_all[6], err_all[6] = read('./hdf5/E1/c4-42_ADE16_200286.hdf5', loc[1])  
    AF_all[7], err_all[7] = read('./hdf5/E1/c4-47_ADE16_200286.hdf5', loc[1]) 	

    #rects2 = ax.bar(ind+width, AF_all, width, color='b', yerr=err_all, ecolor='g')
    plt.plot(ind, AF_all, linestyle='-', linewidth=12, marker='o', markersize=16, color='r', label='ADE16')
    plt.errorbar(ind, AF_all, yerr=err_all, linestyle='None', linewidth=4, marker = 'None',  elinewidth=3, capsize=18, color='r', alpha = 0.5)
    
    
    ########################################### plot figure
    
    #ax.legend( (rects1[0], rects2[0], rects3[0], rects4[0]), ('MTH1', 'CDC55', 'DAL81', 'MUK1') , loc = 'upper left')
    #autolabel(rects1)
    #autolabel(rects2)
    #autolabel(rects3)
    #autolabel(rects4)	
    
    # legend 
    legend = ax.legend(loc='upper left')
    for label in legend.get_texts():
        label.set_fontsize(35)

    for label in legend.get_lines():
        label.set_linewidth(3) 
    
    width = 5
    ax.set_ylim(0, 105)
    ax.set_xlim( -width, ind[N-1] + width)
    ax.set_xticks(ind)
    #ax.set_xticks( ind + 2*width )
    ax.set_xticklabels( ('G7', 'G70', 'G133', 'G196', 'G266', 'G322', 'G385', 'G448'))
    
    #ax.text(5, 80, r'A pair of concomitant mutations', fontsize=30)

    plt.setp(plt.gca().get_xticklabels(), fontsize=40)
    plt.setp(plt.gca().get_yticklabels(), fontsize=42)
    plt.xlabel('Time (generation)', fontsize=45)
    plt.ylabel('Variant Allele Frequency (%)', fontsize=50)
    plt.show()
    #plt.savefig('E1_venn.pdf')
	
'''def autolabel(rects):
    # attach some text labels
    fig, ax = plt.subplots()
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')'''


def read(filename, pos):

    def beta_mean(p):
        return p[0]*1.0/np.sum(p)    

    caseR, caseN, casephi, caseq, loc, refb = rvd3.load_model(filename)
    casegam = caseq['gam']
    caseMu = beta_mean(casegam[pos,:])-casephi['mu0']

    # calculate the lower and upper credible value
    alpha = 0.05
    cred = caseMu*alpha/2
    conf_l = caseMu - cred
    conf_u = caseMu + cred 

    # calculate the error bar value
    err = np.array(cred, cred)
    print 100*caseMu, conf_l, conf_u
	
    return 100*caseMu, 100*err

if __name__ == "__main__":
    main()
