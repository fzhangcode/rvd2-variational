import sys
import os
import numpy as np

import logging
import pdb

import vcf
import xlwt

##os.chdir('S:/yhe2/Research/rvd2/results/2013-09-19_operating_characteristics')

# Insert the src/python/rvd27 directory at front of the path
#rvddir = os.path.join('../../src/python/rvd27')
#sys.path.insert(0, rvddir)
import rvd3

def main():
    book=xlwt.Workbook(encoding="utf-8")
    sheet1=book.add_sheet("TPR_TNR")
    sheet1.write(0, 0, "VAF")
    sheet1.write(0, 1, "Median Depth")

    sheet2=book.add_sheet("Multi-measures")
    sheet2.write(1, 0, "VAF")
    sheet2.write(1, 1, "Median Depth")

    sheet3=book.add_sheet("FDR")
    sheet3.write(0, 0, "VAF")
    sheet3.write(0, 1, "Median Depth")

    sheet4=book.add_sheet("MCC")
    sheet4.write(0, 0, "VAF")
    sheet4.write(0, 1, "Median Depth")

    # method = {'RVD2(T*)(R=6)':'./../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/six_replicates_synthetic_optT/vcf/MCC',
    #           'RVD2(T*)(R=1)':'./../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/one_replicate_synthetic_optT/vcf/MCC',
    #           'RVD2(T=0)(R=6)':'./../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/six_replicates_synthetic_T0/vcf',
    #           'RVD2(T=0)(R=1)':'./../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/one_replicate_synthetic_T0/vcf',
    #           'VarScan2 somatic':'./../2013-09-23_SNP_calling_using_varscan2_somatic/vcf',
    #           'SAMtools':'./../2013-09-10_SNP_calling_using_samtools/vcf',
    #           'GATK':'./../2013-09-13_SNP_calling_using_GATK/vcf',
    #           'MuTect':'./../2013-10-02_SNP_calling_using_MuTect/work',
    #           'Strelka':'./../2013-10-01_SNP_calling_using_strelka/vcf',
    #           'VarScan2 mpileup':'./../2013-09-20_SNP_calling_using_varscan2/vcf'}

    # method = {'RVD2_MCMC(T=0,R=6)':'./../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/six_replicates_synthetic_T0/vcf',
    #           'MuTect':'./../2013-10-02_SNP_calling_using_MuTect/work',
    #           'RVD2_Var(T=0,R=6)':'./vcf'
    #           }

    method = {'RVD3(T=0,R=6)':'./vcf',}
    
    DilutionList = (0.1, 0.3, 1.0, 10.0,100.0)
    DepthList = (10000, 1000, 100, 10)
    i=0
    
    for k, v in method.iteritems():
        i=i+1
        print 'Method %(number)d: %(method)s' %{'number':i, 'method': k}
        sheet1.write(0, i+1, k)
        sheet2.write(0, 9*(i-1)+6, k)
        sheet3.write(0, i+1, k)
        sheet4.write(0, i+1, k)
        character=('Sensitiviy', 'Specificity', 'FPR', 'FNR', 'PPV', 'NPV', 'FDR', 'ACC', 'MCC')
        for j in xrange(9):
            sheet2.write(1,9*(i-1)+j+2,character[j])
        for d in DilutionList:
            if i==1:
                sheet1.write(DilutionList.index(d)*len(DepthList)+1,0,"%0.1f%%" %d)
                sheet2.write(DilutionList.index(d)*len(DepthList)+2,0,"%0.1f%%" %d)
                sheet3.write(DilutionList.index(d)*len(DepthList)+1,0,"%0.1f%%" %d)
                sheet4.write(DilutionList.index(d)*len(DepthList)+1,0,"%0.1f%%" %d)
                
            for r in DepthList:
                # read in the median coverage
                #hdf5Dir='../2013-08-14_Compute_ROC_Synthetic_avg%s' %str(r)
                hdf5Dir = './hdf5/%s' %str(r)
                caseFile = 'Case%s.hdf5' %str(d).replace('.','_')
                caseFile = "%(dir)s/%(file)s" %{'dir':hdf5Dir,'file':caseFile}
               
                # pdb.set_trace()
                #(_, _, _, _, _, caseN,_) = rvd27.load_model(caseFile)
                (_, caseN, _, _, _, _) = rvd3.load_model(caseFile)
                cov = int(np.median(caseN))

                # print the median coverage
                if i==1:
                    sheet1.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+1, 1, "%s" % str(cov))
                    sheet2.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+2, 1, "%s" % str(cov))
                    sheet3.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+1, 1, "%s" % str(cov))
                    sheet4.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+1, 1, "%s" % str(cov))

                # read in called positions from vcf files
                # pdb.set_trace()
                vcfFile=os.path.join(v,"%s" %r,
                                     "vcf%s.vcf" %str(d).replace('.','_'))
                
                logging.debug(vcfFile)

                vcf_reader = vcf.Reader(open(vcfFile, 'r'))
                # pdb.set_trace()
                callpos=np.array([record.POS for record in vcf_reader])

                # prediction classification
                PredictClass = np.zeros(400)
                if len(callpos) != 0:
                    PredictClass[callpos-1] = np.ones_like(callpos)
                    
                # actual classification
                RefClass = np.zeros(400)
                pos = np.arange(85,346,20)
                RefClass[pos-1] = np.ones_like(pos)

                # characteristics computation
                [TPR, TNR, FPR, FNR, PPV, NPV, FDR, ACC, MCC]=characteristics(RefClass, PredictClass)
                ncharacter=(TPR, TNR, FPR, FNR, PPV, NPV, FDR, ACC, MCC)

                # print characteristics
                sheet1.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+1,i+1,"%(TPR)0.2f/%(TNR)0.2f" %{'TPR':TPR,'TNR':TNR})               
                for j in xrange(9):
                    sheet2.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+2,9*(i-1)+j+2,'%0.2f' %ncharacter[j])
                if not np.isnan(FDR):
                    sheet3.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+1,i+1,"%0.2f" %FDR)
                if not np.isnan(FDR):
                    sheet4.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+1,i+1,"%0.2f" %MCC)
    book.save('statistics_no_chi2.xls')  

def characteristics(RefClass = None, PredictClass = None):

    if RefClass is None:
        RefClass = np.zeros(400)
        pos = np.arange(85,346,20)
        RefClass[pos-1] = np.ones_like(pos)
    if PredictClass is None:
        PredictClass = np.copy(RefClass)
        pos = np.arange(85,346,20)
        PredictClass[pos-1] = np.zeros_like(pos)

    #True Positive
    TP = len([i for i in range(len(RefClass)) if RefClass[i]==1 and PredictClass[i]==1])
    #True Negative
    TN = len([i for i in range(len(RefClass)) if RefClass[i]==0 and PredictClass[i]==0])
    #False Positive
    FP = len([i for i in range(len(RefClass)) if RefClass[i]==0 and PredictClass[i]==1])
    #False Negative
    FN = len([i for i in range(len(RefClass)) if RefClass[i]==1 and PredictClass[i]==0])
  
    #RefClassence Postive
    P1 = sum(RefClass)
    #RefClassence Negative
    N1 = len(RefClass) - P1

    #PredictClassion Postive
    P2 = sum(PredictClass)
    #PredictClassion Negative
    N2 = len(PredictClass) - P2

    #Sensitivity (TPR,true positive rate)
    if TP+FN != 0:
        TPR = float(TP)/(TP+FN)
    else:
        TPR=np.nan 
    #Specificity (TNR, true negative rate)
    if FP+TN != 0:
        TNR = float(TN)/(FP+TN)
    else:
        TNR = np.nan

    #FPR (1-Specificity, false negative rate)
    if FP+TN != 0:
        FPR = float(FP)/(FP+TN)
    else:
        FPR = np.nan
    #FNR
    if TP+FN !=0:
        FNR = float(FN)/(TP+FN)
    else:
        FNR = np.nan

    #PPV (Positive PredictClassive value)
    #pdb.set_trace()
    if TP+FP !=0:
        PPV = float(TP)/(TP+FP)
    else:
        PPV = np.nan
    #NPV (Negative PredictClassive value)
    if TN+FN !=0:
        NPV = float(TN)/(TN+FN)
    else:
        NPV = np.nan

    #FDR (false discovery rate)
    if FP+TP!=0:
        FDR = float(FP)/(FP+TP)
    else:
        FDR = np.nan

    #ACC (Accuracy)
    ACC = (TP+TN)/(P1+N1)
    #MCC (Matthews correlation coefficient)
    if P1*N1*P2*N2 != 0:
        MCC = float(TP*TN-FP*FN)/np.sqrt(P1*N1*P2*N2)
    else:
        MCC= np.nan
    # AUC

    return  TPR, TNR, FPR, FNR, PPV, NPV, FDR, ACC, MCC

if __name__ == '__main__':
    main()
