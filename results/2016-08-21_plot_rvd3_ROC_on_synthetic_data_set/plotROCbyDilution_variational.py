import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import logging
import scipy.stats as ss
import rvd3


# Insert the rvd29 directory at front of the path
rvddir = os.path.join('./')
sys.path.insert(0, rvddir)


logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0)

    folderList = ('../2015-09-28_Run_rvd3_synthetic_data_set/hdf5/10',\
                   '../2015-09-28_Run_rvd3_synthetic_data_set/hdf5/100',\
                   '../2015-09-28_Run_rvd3_synthetic_data_set/hdf5/1000',\
                   '../2015-09-28_Run_rvd3_synthetic_data_set/hdf5/10000')
    
    N=1000 # Z sampling size  
    fig=plt.figure(figsize=(15, 14))
    chi2=False
    lstyle=('--','-.','-',':')
    lcolor=('c','r','g','b')
    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        ax = fig.add_subplot(2,2,dilutionList.index(d)+1)
        label=[]
        nscrewc = -1
        for f in folderList:
            controlFile = "%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
            # ROC
            [fpr,tpr,cov] = ROCpoints(controlFile,caseFile, d, N, P=0.95, chi2=chi2)
            ax.plot(fpr+0.01*nscrewc,tpr-0.01*nscrewc,linestyle=lstyle[folderList.index(f)], color=lcolor[folderList.index(f)], linewidth=3, label='%d' % cov)
            #ax.plot(fpr,tpr,marker='o', markersize = 5, markerfacecolor=lcolor[folderList.index(f)])

        l = ax.legend(loc=4,prop={'size':20},title='Read depth')
        l.get_title().set_fontsize(25)
        ax.plot([0,1],[0,1],color='k',linestyle='dashed')
        ax.set_title('%0.1f%% Mutant Mixture' % d, fontsize=25)
        ax.set_xlim((-0.03,1.03))
        ax.set_ylim((-0.03,1.03))
                 
        ax.set_xlabel('1-Specificity (FPR)',fontsize=25)
        ax.set_ylabel('Sensitivity (TPR)',fontsize=25)

        plt.setp(plt.gca().get_xticklabels(), fontsize=25)
        plt.setp(plt.gca().get_yticklabels(), fontsize=25)

    plt.tight_layout()  
	
    if chi2:
        title='ROC_with_chi2'
    else:
        title='ROC_without_chi2'
    figformat='.pdf'
    plt.savefig(title+figformat)


def ROCpoints(controlFile,caseFile, d, N, P, chi2):

    # Load the model samples
    controlR, controlN, controlphi, controlq, controlLoc, _ = rvd3.load_model(controlFile)
    controlgam = controlq['gam']

    caseR, caseN, casephi, caseq, caseLoc, refb = rvd3.load_model(caseFile)
    casegam = caseq['gam'] 


    #(N,J) = np.shape(caseR)[0:2]
    J = len(controlLoc)
    def beta_mean(p):
        return p[0]*1.0/np.sum(p)    

    def beta_var(p):
        s = np.sum(p)
        return p[0]*p[1]/(s**2*(s+1))
    
    # Draw random samples from Beta distribution    
    controlMu = np.zeros(shape=(J, 4000))
    caseMu = np.zeros(shape=(J, 4000))
    for j in xrange(J):
        controlMu[j] = np.random.beta(controlgam[j,:][0], controlgam[j,:][1], 4000)   
        caseMu[j] = np.random.beta(casegam[j,:][0], casegam[j,:][1], 4000) 		


    # Extract the common locations in case and control
    caseLocIdx = [i for i in xrange(len(caseLoc)) if caseLoc[i] in controlLoc]
    controlLocIdx = [i for i in xrange(len(controlLoc)) if controlLoc[i] in caseLoc]
	

    caseMu = caseMu[caseLocIdx,:]
    controlMu = controlMu[controlLocIdx,:]
    # caseR = caseR[:,caseLocIdx,:]
    # controlR = controlR[:,controlLocIdx,:]
    # caseN = caseN[:,caseLocIdx]
    # controlN = controlN[:,controlLocIdx]

    loc = caseLoc[caseLocIdx]
    J = len(loc)
    pos = np.arange(85,346,20)
    posidx = [i for i in xrange(J) if int(loc[i][8:]) in pos]
    
    # Sample from the posterior Z = muCase - muControl        
    (Z, caseMuS, controlMuS) = sample_post_diff(caseMu-casephi['mu0'], controlMu-controlphi['mu0'], N)
	
    # Compute cumulative posterior probability for regions (Threshold,np.inf)
    T = np.linspace(np.min(np.min(Z)), np.max(np.max(Z)), num=300)

    pList = [bayes_test(Z, [(t, np.inf)]) for t in T]

    # mutation classification
    clsList = np.array((np.array(pList)>P).astype(int))
    clsList = clsList.reshape((clsList.shape[0],clsList.shape[1]))# category list
    
    # chi2 test for goodness-of-fit to a uniform distribution for non-ref bases
    if chi2:
        nRep = caseR.shape[0]
        chi2Prep = np.zeros((J,nRep))
        chi2P = np.zeros(J)
        for j in xrange(J):
                chi2Prep[j,:] = np.array([rvd3.chi2test( caseR[i,j,:] ) for i in xrange(nRep)] )
                if np.any(np.isnan(chi2Prep[j,:])):
                    chi2P[j] = 1
                else:
                   chi2P[j] = 1-ss.chi2.cdf(-2*np.sum(np.log(chi2Prep[j,:] + np.finfo(float).eps)), 2*nRep) # combine p-values using Fisher's Method
        
        clsList2 = np.array((np.array(chi2P)<0.05/J).astype(int))        
        clsList2 = np.tile(clsList2,(clsList.shape[0],1))
        clsList = np.array(((clsList+clsList2)==2).astype(int))

    # false postive rate
    fpr = np.array([float(sum(clsList[i])-sum(clsList[i,np.array(posidx)]))/(clsList.shape[1]-len(posidx)) for i in xrange(clsList.shape[0])])

    # true positive rate
    tpr = np.array([float(sum(clsList[i,np.array(posidx)]))/len(posidx) for i in xrange(clsList.shape[0])])

    cov = np.median(caseN)

    # # return information for mu bar plot at called positions under optimal threshold.

    # # using EL distance.
# ##    distance=np.sum(np.power([fpr,tpr-1],2),0) 
# ##    Tidx=distance.argmin()
# ##    print Tidx


     # # Using L1 distance 
    # distance = 1+tpr-fpr
    # Tidx=distance.argmax()
    

    # outputFile=os.path.join(path,'vcf%s.vcf' %str(d).replace('.','_'))
    
    # with h5py.File(controlFile, 'r') as f:
        # refb = f['/refb'][...]
        # f.close()
    # refb = refb[controlLocIdx]
    
    # altb = []
    # call=[]
    # acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    # for i in xrange(J):
        # r = np.squeeze(caseR[:,i,:]) # replicates x bases
        
        # # Make a list of the alternate bases for each replicate
        # acgt_r = ['A','C','G','T']
        # del acgt_r[ acgt[refb[i]] ]

        # altb_r = [acgt_r[x] for x in np.argmax(r, axis=1)]

        # if clsList[Tidx,i]==1:
            # call.append(True)
            # altb.append(altb_r[0])
        # else:
            # altb.append(None)
            # call.append(False)
            
   # rvd30.write_vcf(outputFile, loc, call, refb, altb, np.mean(caseMu, axis=1), np.mean(controlMu, axis=1))
    return fpr,tpr, cov

	
	
	
def sample_post_diff(muCaseG, muControlG, N):
    """ Return N samples from the posterior distribution for 
         u_j|r_case - u_j|r_control. """
    
    nCase = muCaseG.shape[1]
    nControl = muControlG.shape[1]
    
    caseSample = np.random.choice(nCase, size=N, replace=True)
    controlSample = np.random.choice(nControl, size=N, replace=True)
    
    muCaseS = muCaseG[:, caseSample]
    muControlS = muControlG[:, controlSample]

    Z = muCaseS - muControlS

    return (Z, muCaseS, muControlS)
 
def bayes_test(Z, roi, type = 'close'):
    """ Return posterior probabilities in regions defined in list of tuples (roi)
        from samples in columns of Z. """
    (J,N)=np.shape(Z)
    
    nTest = len(roi) # get the number of regions to compute probabilities 

    p = np.zeros((J,nTest))
    for i in xrange(nTest):
        for j in xrange(J):
                if type == 'close':
                    # somatic test
                    p[j,i] = np.float( np.sum( np.logical_and( (Z[j,:] >= roi[i][0]), (Z[j,:] <= roi[i][1]) ) ) ) / N
                elif type == 'open':
                    # diff test
                    p[j,i] = np.float( np.sum( np.logical_and( (Z[j,:] > roi[i][0]), (Z[j,:] < roi[i][1]) ) ) ) / N
    
    p = np.sum(p,1) # add up all the regions in each position.
    
    return p # combine probabilities from all regions
	
	
	
if __name__ == '__main__':
    main()

