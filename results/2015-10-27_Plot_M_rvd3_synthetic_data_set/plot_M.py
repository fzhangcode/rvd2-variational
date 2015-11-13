import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb

rvddir = os.path.join('./')
sys.path.insert(0, rvddir)
import rvd3

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)
    
    folder = '2015-09-28_Run_rvd3_synthetic_data_set/hdf5/10'

    fig=plt.figure(figsize=(12,20))
    #plt.suptitle('Read depth/M across position')

    controlFile = "../%s/Control.hdf5" %folder
    controlR, controlN, controlPhi, controlq, controlLoc, _ = rvd3.load_model(controlFile)
	
    sub0=len(dilutionList)+1
    ax=fig.add_subplot(sub0,2,1)
    #TODO: use index of controlN rather than directly controlLOC
    controlLoc = [int(x.split(':')[1]) for x in controlLoc]
    ax.plot(controlLoc,controlN.T)
    ax.set_title('Control')
    ax.set_ylabel('Coverage')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    ax=fig.add_subplot(sub0,2,2)
    ax.semilogy(controlLoc,np.mean(controlPhi['M'],axis=1))
    ax.set_title('Control')
    ax.set_ylabel('M')
    #ax.set_ylim([1e-4,1e5])
    ax.semilogy([controlLoc[0],controlLoc[-1]],[controlPhi['M0'],controlPhi['M0']],color='r',ls='--')
    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        caseFile = "Case%s.hdf5" % str(d).replace(".","_")
        caseFile = "../%(folder)s/%(file)s" %{'folder':folder,'file':caseFile}
        caseR, caseN, casePhi, caseq, caseLoc, _ = rvd3.load_model(caseFile)
        ax=fig.add_subplot(sub0,2,2*dilutionList.index(d)+3)
        caseLoc = [int(x.split(':')[1]) for x in caseLoc]
        ax.plot(caseLoc,caseN.T)
        ax.set_title("Dilution %0.1f%%" % d)
        if dilutionList.index(d)==len(dilutionList)-1:
            ax.set_xlabel('Position')
        ax.set_ylabel('Coverage')
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax=fig.add_subplot(sub0,2,2*dilutionList.index(d)+4)
        ax.semilogy(caseLoc,np.mean(casePhi['M'],axis=1))
        ax.set_title("Dilution %0.1f%%" % d)
        if dilutionList.index(d)==len(dilutionList)-1:
            ax.set_xlabel('Position')
        ax.set_ylabel('M')
        #ax.set_ylim([1e-4,1e5])
        ax.semilogy([caseLoc[0],caseLoc[-1]],[casePhi['M0'],casePhi['M0']],color='r',ls='--')
    plt.savefig('M_dsample=10.png')

if __name__ == '__main__':
    main()

