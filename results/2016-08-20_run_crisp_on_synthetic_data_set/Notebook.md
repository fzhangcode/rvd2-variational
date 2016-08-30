2016-08-20 run crisp on synthetic data set
==============================

Purpose
------------
Run CRISP on the synthetic data set.

Conclusions
-----------------
CRISP works very well in the case of VAF = 100%.  
It failed to call true positives when the VAF is 1.0%, 03%, and 0.1%.
   
Background
-----------------


Materials and Equipment
------------------------------


Experimental Protocol
---------------------------
1. download crisp from: https://sites.google.com/site/vibansal/software/crisp 
2. run crisp

100

	[fzhang@redwood CRISP-071812]$ ./CRISP --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.12_GTT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.12_GTT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.07_CGT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.07_CGT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.14_TCT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.14_TCT.reheadered.sorted.bam  --ref ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa --VCF variantcalls_downsample=10_vaf=100.VCF --poolsize 6 > variantcalls_downsample=10_vaf=100.log
	
	CRISP options: QVoffset 33 min_base_quality 10 min_mapping_score 20 max_permutations 20000 poolsize 6 CT-pvalue-thresh -3.5 QV-pvalue-thresh -5.0
	
	reading fasta index file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa.fai 1
	fasta file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa has 1 chromosomes/contigs
	
	processing 6 bamfiles: ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.12_GTT.reheadered.sorted.bam ..... ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.14_TCT.reheadered.sorted.bam 
	reading next chromosome 1_1_400 
	processed 1000000 reads QSIZE 244865 1_1_400:164 203
	processed 2000000 reads QSIZE 271972 1_1_400:302 341
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.12_GTT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.07_CGT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.07_CGT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.12_GTT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.14_TCT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.14_TCT.reheadered.sorted.bam 

10

	[fzhang@redwood CRISP-071812]$ ./CRISP --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.02_ACT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.02_ACT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.04_ATT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.04_ATT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.05_CAT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.05_CAT.reheadered.sorted.bam  --ref ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa --VCF variantcalls_downsample=10_vaf=10.VCF --poolsize 6 > variantcalls_downsample=10_vaf=10.log
	
	CRISP options: QVoffset 33 min_base_quality 10 min_mapping_score 20 max_permutations 20000 poolsize 6 CT-pvalue-thresh -3.5 QV-pvalue-thresh -5.0
	
	reading fasta index file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa.fai 1
	fasta file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa has 1 chromosomes/contigs
	
	processing 6 bamfiles: ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.02_ACT.reheadered.sorted.bam ..... ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.05_CAT.reheadered.sorted.bam 
	reading next chromosome 1_1_400 
	processed 1000000 reads QSIZE 208632 1_1_400:181 219
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.02_ACT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.04_ATT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.02_ACT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.05_CAT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p1.04_ATT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c3_p2.05_CAT.reheadered.sorted.bam

1

	[fzhang@redwood CRISP-071812]$ ./CRISP --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.07_CGT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.07_CGT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.12_GTT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.12_GTT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.14_TCT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.14_TCT.reheadered.sorted.bam  --ref ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa --VCF variantcalls_downsample=10_vaf=1.VCF --poolsize 6 > variantcalls_downsample=10_vaf=1.log
	
	CRISP options: QVoffset 33 min_base_quality 10 min_mapping_score 20 max_permutations 20000 poolsize 6 CT-pvalue-thresh -3.5 QV-pvalue-thresh -5.0
	
	reading fasta index file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa.fai 1
	fasta file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa has 1 chromosomes/contigs
	
	processing 6 bamfiles: ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.07_CGT.reheadered.sorted.bam ..... ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.14_TCT.reheadered.sorted.bam 
	reading next chromosome 1_1_400 
	processed 1000000 reads QSIZE 413143 1_1_400:76 115
	processed 2000000 reads QSIZE 411406 1_1_400:172 211
	processed 3000000 reads QSIZE 482389 1_1_400:253 292
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.14_TCT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.12_GTT.reheadered.sorted.bam 
	processed 4000000 reads QSIZE 380748 1_1_400:325 364
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.12_GTT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.07_CGT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.07_CGT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.14_TCT.reheadered.sorted.bam

0.3

	[fzhang@redwood CRISP-071812]$ ./CRISP --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.02_ACT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.02_ACT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.04_ATT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.04_ATT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.05_CAT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.05_CAT.reheadered.sorted.bam  --ref ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa --VCF variantcalls_downsample=10_vaf=03.VCF --poolsize 6 > variantcalls_downsample=10_vaf=03.log
	
	CRISP options: QVoffset 33 min_base_quality 10 min_mapping_score 20 max_permutations 20000 poolsize 6 CT-pvalue-thresh -3.5 QV-pvalue-thresh -5.0
	
	reading fasta index file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa.fai 1
	fasta file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa has 1 chromosomes/contigs
	
	processing 6 bamfiles: ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.02_ACT.reheadered.sorted.bam ..... ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.05_CAT.reheadered.sorted.bam 
	reading next chromosome 1_1_400 
	processed 1000000 reads QSIZE 292093 1_1_400:128 167
	processed 2000000 reads QSIZE 398921 1_1_400:229 269
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.02_ACT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.04_ATT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.02_ACT.reheadered.sorted.bam 
	processed 3000000 reads QSIZE 272943 1_1_400:324 363
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p2.05_CAT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.04_ATT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c2_p1.05_CAT.reheadered.sorted.bam 

0.1

	[fzhang@redwood CRISP-071812]$ ./CRISP --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p1.07_CGT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p2.07_CGT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p1.12_GTT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p2.12_GTT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p1.14_TCT.reheadered.sorted.bam --bam ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p2.14_TCT.reheadered.sorted.bam  --ref ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa --VCF variantcalls_downsample=10_vaf=01.VCF --poolsize 6 > variantcalls_downsample=10_vaf=01.log
	
	CRISP options: QVoffset 33 min_base_quality 10 min_mapping_score 20 max_permutations 20000 poolsize 6 CT-pvalue-thresh -3.5 QV-pvalue-thresh -5.0
	
	reading fasta index file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa.fai 1
	fasta file ../../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa has 1 chromosomes/contigs
	
	processing 6 bamfiles: ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p1.07_CGT.reheadered.sorted.bam ..... ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p2.14_TCT.reheadered.sorted.bam 
	reading next chromosome 1_1_400 
	processed 1000000 reads QSIZE 338457 1_1_400:96 135
	processed 2000000 reads QSIZE 356062 1_1_400:217 256
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p2.07_CGT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p1.12_GTT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p2.14_TCT.reheadered.sorted.bam 
	processed 3000000 reads QSIZE 287470 1_1_400:325 364
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p2.12_GTT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p1.07_CGT.reheadered.sorted.bam 
	finished reading bam file ../../../../rvd/results/2013-10-02_problematic_header_removal/bam/10/20100916_c1_p1.14_TCT.reheadered.sorted.bam 

Change the name of the vcf files to variantcalls_vcf0_1.VCF, variantcalls_vcf0_3.VCF variantcalls_vcf1_0.VCF variantcalls_vcf10_0.VCF variantcalls_vcf100_0.VCF.


Results
----------- 


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _______Fan Zhang_______     Date: ___________2016/08/20_________


Witnessed by: ________________________
